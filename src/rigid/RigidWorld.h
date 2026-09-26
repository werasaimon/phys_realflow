#pragma once
// Rigid body dynamics.
//  * shapes: sphere, box, convex polyhedron (ConvexShape hierarchy, support mappings)
//  * broad phase: BVH over body AABBs
//  * narrow phase: analytic / SAT / GJK+EPA with perturbation manifolds (NarrowPhase)
//  * static environment: domain walls (planes) and a triangle mesh (BVH -> triangles -> GJK/EPA)
//  * solver: sequential impulses with Coulomb friction on persistent manifolds with warm starting

#include "spatial/BVH.h"
#include "math/Math.h"
#include "rigid/BroadPhase.h"
#include "rigid/Ccd.h"
#include "rigid/Joints.h"
#include "rigid/RigidBody.h"
#include "rigid/Narrowphase.h"
#include "rigid/Shapes.h"

#include <memory>
#include <unordered_map>
#include <vector>

namespace rf {

enum class RigidSolver {
    XPBD,              // Mueller et al. 2020, "Detailed Rigid Body Simulation with Extended PBD" (experimental)
    SequentialImpulse  // Catto-style impulses: warm starting, split impulse, shock propagation (default)
};

struct RigidParams {
    RigidSolver solver = RigidSolver::SequentialImpulse;
    Vector3 gravity{0, -9.81f, 0};
    // Many small steps with few iterations (TGS / "small steps", as PhysX and Jolt do) keep tall
    // stacks stable where many iterations with one large step do not.
    int iterations = 6;         // sequential-impulse iterations (XPBD uses one per substep)
    int substeps = 10;          // per frame (XPBD wants 20-40)
    int collisionInterval = 4;  // XPBD: run collision detection every n substeps (anchors persist)
    int positionIterations = 4; // XPBD: contact position passes per substep
    int manifoldIterations = 4; // local relaxation of each manifold's normal impulses (block solve)
    float baumgarte = 0.2f;
    float slop = 0.004f;
    float contactMargin = 0.01f; // speculative contact distance [m]
    float linearDamping = 0.02f;
    float angularDamping = 0.05f;
    bool collideWithDomain = true;
    bool warmStarting = true;
    bool splitImpulse = true;   // penetration recovery on pseudo velocities (no energy gain)
    // Shock propagation (Guendelman, Bridson, Fedkiw 2003): final passes solve contacts level by
    // level from the ground up, treating the lower body as infinitely heavy, so impulses travel
    // in one direction only and tall stacks behave like short ones.
    bool shockPropagation = true;
    float rollingResistance = 0.02f; // fraction of the normal load resisting relative rolling
    bool rotationalLock = true;      // SO(3) angular constraint for resting face-face contacts
    bool blockSolver = true;         // exact LCP for each manifold's normal impulses (<= 4 points)
    float blockCfm = 1e-3f;          // relative regularisation making 4 coplanar points well posed
    // Sleeping islands: an island whose bodies all stay below the thresholds for sleepTime sleeps.
    bool sleeping = true;
    float sleepLinear = 0.05f;       // [m/s]
    float sleepAngular = 0.1f;       // [rad/s]
    float sleepTime = 0.5f;          // [s]
    float stickVelocity = 0.01f;     // below this sliding speed static friction applies [m/s]
    // Resting-contact damping: residual motion of touching bodies below these thresholds is
    // removed by an opposing velocity change of at most restDamping * dt.
    // Off by default: removing velocity after the solve fights the support of stacked bodies.
    float restDamping = 0.0f;           // [m/s^2]
    float restLinearThreshold = 0.05f;  // [m/s]
    float restAngularThreshold = 0.2f;  // [rad/s]
    int shockIterations = 2;
    bool shockFriction = true;
    int jointPositionIterations = 3; // nonlinear Gauss-Seidel passes of the joint position stage
    // Continuous collision detection (conservative advancement on GJK) for fast bodies.
    bool ccd = true;
    float ccdThreshold = 0.5f;   // CCD when the motion in a step exceeds half of the body's smallest extent (as Bullet)
    float ccdTolerance = 0.002f; // [m] distance at which the time of impact is accepted       // friction in the shock pass: supports hold what rests on them
};

class RigidWorld {
public:
    RigidParams params;

    void clear() { bodies_.clear(); cache_.clear(); manifolds_.clear(); joints_.clear(); grab_ = GrabJoint(); }

    // Joints (b = -1: attached to the static world). Anchors/axes are given in world space at the
    // current poses; the joint stores them in the bodies' local frames.
    BallJoint& addBallJoint(int a, int b, const Vector3& worldAnchor);
    HingeJoint& addHingeJoint(int a, int b, const Vector3& worldAnchor, const Vector3& worldAxis);
    SliderJoint& addSliderJoint(int a, int b, const Vector3& worldAxis);
    FixedJoint& addFixedJoint(int a, int b);
    DistanceJoint& addDistanceJoint(int a, int b, const Vector3& worldAnchorA, const Vector3& worldAnchorB);
    const std::vector<std::unique_ptr<Joint>>& joints() const { return joints_; }
    int addSphere(const Vector3& pos, float radius, float density, const Vector3& color);
    int addBox(const Vector3& pos, const Vector3& halfExtents, const Quaternion& rot, float density, const Vector3& color);
    // Convex polyhedron from a closed convex mesh; `pos` is where the mesh's centre of mass goes.
    int addConvex(const TriMesh& convexMesh, const Vector3& pos, const Quaternion& rot, float density, const Vector3& color);
    // Non-convex body from a prepared compound shape (`pos` = where its centre of mass goes, `rot`
    // = orientation of the original model).
    int addCompound(std::shared_ptr<const CompoundShape> shape, const Vector3& pos, const Quaternion& rot, float density, const Vector3& color);
    int addBody(std::shared_ptr<const ConvexShape> shape, const Vector3& pos, const Quaternion& rot, float density, const Vector3& color);

    void setDomain(const AABB& d) { domain_ = d; }
    void setStaticMesh(const MeshBVH* bvh) { mesh_ = bvh; }
    void setBroadPhase(std::unique_ptr<BroadPhase> bp) { broadphase_ = std::move(bp); }
    const BroadPhase& broadPhase() const { return *broadphase_; }
    size_t pairCount() const { return pairs_.size(); }
    int bodyLevel(int i) const { return i < int(levels_.size()) ? levels_[i] : -1; }

    void step(float dt);
    void applyImpulse(int body, const Vector3& impulse, const Vector3& worldPoint);
    void applyBiasImpulse(int body, const Vector3& impulse, const Vector3& worldPoint);
    // Impulse from outside the solver (fluid coupling, user input): also wakes the body up.
    void applyExternalImpulse(int body, const Vector3& impulse, const Vector3& worldPoint);
    // Linear impulse through the centre of mass plus an angular impulse (e.g. gas pressure).
    void applyExternalWrench(int body, const Vector3& impulse, const Vector3& angularImpulse);
    // Velocity change from outside the solver (particle coupling): also wakes the body up.
    void applyVelocityChange(int body, const Vector3& dv, const Vector3& dw);
    void wake(int body);
    size_t sleepingCount() const;
    size_t ccdHits() const { return ccdHits_; } // bodies clamped by CCD in the last step

    // Ray cast against all bodies: nearest hit (body index, distance along the unit direction).
    bool raycast(const Vector3& origin, const Vector3& dir, float maxT, int& body, float& t, Vector3& normal) const;

    // Mouse joint (Box2D b2MouseJoint): a soft point constraint pulling a body point to a target.
    struct GrabJoint {
        bool active = false;
        int body = -1;
        Vector3 localAnchor;  // grab point in the body frame
        Vector3 target;       // world target
        float frequency = 5.0f, damping = 0.7f, maxForce = 1000.0f; // maxForce per unit weight
        Vector3 impulse;      // accumulated, warm started
    };
    void grab(int body, const Vector3& worldPoint);
    void setGrabTarget(const Vector3& target) { grab_.target = target; }
    void releaseGrab() { grab_ = GrabJoint(); }
    const GrabJoint& grabJoint() const { return grab_; }

    std::vector<RigidBody>& bodies() { return bodies_; }
    const std::vector<RigidBody>& bodies() const { return bodies_; }
    size_t contactCount() const { return contactCount_; }
    size_t manifoldCount() const { return manifolds_.size(); }
    float kineticEnergy() const;
    // Time spent per stage in the last step [ms].
    struct Timings { float broad = 0, narrow = 0, solve = 0, ccd = 0, islands = 0, total = 0; };
    const Timings& timings() const { return timings_; }

    // Contact points of the last step (for visualisation / debugging).
    struct DebugContact {
        int a, b;
        Vector3 position, normal;
        float depth, impulse;
    };
    std::vector<DebugContact> debugContacts() const;

private:
    struct SolverPoint {
        Vector3 position, normal, localA; // localA: anchor in A's frame, used to match points between steps
        uint64_t id = 0;               // hash of the quantised localA
        float depth = 0;
        Vector3 t1, t2;
        float massN = 0, massT1 = 0, massT2 = 0;
        float velocityBias = 0; // restitution / speculative gap (velocity level)
        float positionBias = 0; // penetration recovery (split impulse, position level)
        float jn = 0, jt1 = 0, jt2 = 0, jp = 0;
    };
    struct Manifold {
        int a = -1, b = -1; // b < 0: static (walls / mesh)
        float friction = 0.5f, staticFriction = 0.7f, restitution = 0.2f, rolling = 0.0f;
        std::vector<SolverPoint> points;
        // Manifold-level friction, as in ReactPhysics3D: two tangents and twist at the centre of
        // the contact patch, plus rolling resistance. Accumulated impulses are warm started.
        Vector3 center, normal, t1, t2;
        float massT1 = 0, massT2 = 0, massTwist = 0, patchRadius = 0;
        Matrix3x3 rollMass = Matrix3x3::zero();
        float jt1 = 0, jt2 = 0, jtwist = 0;
        Vector3 jroll;
        // Rotational lock of a resting face contact (angular part of a fixed joint on SO(3)).
        bool locked = false;
        Quaternion lockRef;        // reference relative orientation qB^-1 * qA
        Vector3 lockError;      // log(E) in world frame, releases the lock when large
        Vector3 jlock;
        float K[4][4] = {};  // normal effective-mass matrix of the points (block solver)
        int activeSet = -1;  // block solver: active points of the last solve (tried first)
    };
    // What survives from one step to the next for a body pair.
    struct CachedPair {
        std::vector<SolverPoint> points;
        Vector3 friction;     // world tangential impulse at the patch centre
        float twist = 0;   // impulse moment about the normal
        Vector3 roll;         // rolling resistance impulse moment
        bool locked = false;
        Quaternion lockRef;
        Vector3 lock;         // rotational lock impulse moment
    };
    static uint64_t key(int a, int b) { return (uint64_t(uint32_t(a)) << 32) | uint32_t(b + 64); }

    // Warm-start cache: each point carries a hash of its quantised position in A's frame (id), so
    // the same contact is found again next step without per-point searches.
    static float contactCell(const RigidBody& A);
    static uint64_t positionHash(const Vector3& localA, float cell);
    static const SolverPoint* findCached(const std::vector<SolverPoint>& old, const SolverPoint& p, float cell);

    void collide();
    void collideStatic(int i, std::vector<Manifold>& out) const;
    void addManifold(std::vector<Manifold>& out, int a, int b, ContactManifold& cm) const;
    void prepare(float dt);
    void solve();
    void buildColors();
    template <class F> void forEachManifold(F&& f);
    void prepareManifold(Manifold& m, float dt);
    void solveManifold(Manifold& m);
    void computeLevels();
    void blockNormalSolve(Manifold& m);
    // Sleeping: frozen bodies act as static during the step.
    void freezeSleepers();
    void unfreezeAll(bool onlyAwake = false);
    bool updateIslands(bool decideSleep, float dt);
    struct Frozen { int body; float invMass; Vector3 invInertiaLocal; };
    std::vector<Frozen> frozen_;
    GrabJoint grab_;
    std::vector<std::unique_ptr<Joint>> joints_;
    template <class J> J& attach(std::unique_ptr<J> j, const Vector3& anchorA, const Vector3& anchorB, const Vector3& axis);
    void solveJointPositions();
    void continuousCollision();
    size_t ccdHits_ = 0;
    std::vector<char> ccdClamped_; // bodies stopped by CCD in the last step (their next contact is an impact)
    Matrix3x3 grabMass_ = Matrix3x3::zero();
    Vector3 grabBias_;
    float grabGamma_ = 0;
    void prepareGrab(float dt);
    void solveGrab(float dt);
    std::vector<int> islandParent_;
    int nextIsland_ = 0;
    bool wakeIsland(int body); // wakes every body of the sleeping island; true if anything woke

    // --- XPBD solver (XpbdSolver.cpp) ---
    struct XContact {
        int a = -1, b = -1;
        Vector3 rA, rB;      // anchors: local to A / local to B (world point if B is static)
        Vector3 n;           // from B to A
        float lambdaN = 0, lambdaT = 0;
        float vnPrev = 0;
        float friction = 0.5f, staticFriction = 0.7f, restitution = 0.2f;
        bool active = false;
    };
    void stepXPBD(float h);
    void buildXContacts();
    void solveXContactPosition(XContact& c, float h);
    void solveXContactVelocity(XContact& c, float h);
    float applyPositionalCorrection(int a, int b, const Vector3& corr, const Vector3& pA, const Vector3& pB);
    void applyVelocityChange(int a, int b, const Vector3& dv, const Vector3& pA, const Vector3& pB);
    Vector3 anchorA(const XContact& c, bool prev = false) const;
    Vector3 anchorB(const XContact& c, bool prev = false) const;
    std::vector<XContact> xcontacts_;
    int substepCounter_ = 0;
    void solveManifoldShock(Manifold& m);

    std::vector<RigidBody> bodies_;
    std::vector<Manifold> manifolds_;               // current step
    std::unordered_map<uint64_t, CachedPair> cache_; // previous step (warm starting), by body pair
    std::vector<std::vector<int>> colors_; // manifold batches without shared bodies
    bool parallelColors_ = false;
    std::vector<int> levels_; // contact-graph distance from the static environment
    std::unique_ptr<BroadPhase> broadphase_ = std::make_unique<SweepAndPruneBroadPhase>(0.1f);
    std::vector<std::pair<int, int>> pairs_;
    NarrowPhase narrow_;
    AABB domain_{{-1, 0, -1}, {1, 2, 1}};
    const MeshBVH* mesh_ = nullptr;
    size_t contactCount_ = 0;
    float lastDt_ = 1.0f / 600.0f;
    bool shockFrictionPass_ = false;
    Timings timings_;
};

} // namespace rf
