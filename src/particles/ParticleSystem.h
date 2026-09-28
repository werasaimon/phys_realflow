#pragma once
// Unified particle solver (as NVIDIA FleX: Macklin, Mueller, Chentanez, Kim 2014, "Unified Particle
// Physics for Real-Time Applications"). All particles share one neighbour grid and one constraint
// loop; each has a phase:
//   fluid - incompressible liquid: Position Based Fluids (Macklin & Mueller 2013), an SPH method
//   soft  - soft body: shape matching on overlapping clusters (SoftBody.h)
//   cloth - cloth: XPBD stretch / shear / bending constraints (Cloth.h)
// Particles of different phases (and distant particles of the same cloth) collide with each other.
//  * poly6 density / spiky gradient kernels, uniform-grid neighbour search
//  * unilateral density constraint (no clumping at free surfaces)
//  * XSPH viscosity, vorticity confinement
//  * boundaries: domain box, static meshes (BVH signed distance), two-way coupling with rigid bodies

#include "spatial/BVH.h"
#include "particles/Cloth.h"
#include "particles/SoftBody.h"
#include "math/Math.h"
#include "rigid/RigidWorld.h"

#include <functional>
#include <vector>

namespace rf {

enum class ParticlePhase : uint8_t { Fluid, Soft, Cloth };

// The hot gas around a point as a burning cloth feels it: its temperature [K above ambient] and the
// thermal radiation of the flame arriving there [W/m^2].
struct GasHeat {
    float temperature = 0;
    float irradiance = 0;
};

// What a burning cloth particle gives to the gas around it (see ParticleSystem::burnCloths).
struct FireOutput {
    Vector3 position;
    float heat = 0; // [J]
    float fuel = 0; // [m^3] fuel gas at unit concentration
};

struct ParticleParams {
    float particleRadius = 0.015f;  // [m], spacing = 2r, kernel radius h = 4r   (*reset)
    float restDensity = 1000.0f;    // [kg/m^3]
    int solverIterations = 4;
    int solidIterations = 2;        // passes of the soft / cloth / contact constraints per iteration
    // Stiff stacks (Macklin et al. 2014, "Unified Particle Physics for Real-Time Applications",
    // sec. 5.2, eq. 21): in the contacts between particles a particle counts as lighter the higher
    // it sits, m* = m exp(-k h), with h its height along -gravity in particle spacings. A contact
    // then lifts the upper particle and leaves the lower one where it is, so a pile of soft bodies
    // is held up within a few passes instead of sinking into itself. 0 = off.
    float stackMassScaling = 1.0f;
    int clothSubsteps = 8;          // small steps of the cloth inside every substep (Macklin et al. 2019)
    float clothSpacing = 1.0f;      // distance between cloth particles, in particle radii: 1 = a
                                    // particle of radius r cannot slip through the sheet (as in FleX)
    int substeps = 3;               // per frame
    float viscosity = 0.02f;        // XSPH coefficient, 0..1
    float vorticity = 0.0f;         // confinement strength [m/s^2 scale]
    float relaxation = 0.1f;        // constraint-force-mixing term (scaled by 1/h^2)
    float tensileK = 0.03f;         // artificial pressure k (Macklin & Muller 2013: 0.1), dimensionless, x h^2 inside
    float wallFriction = 0.1f;      // 0 = free slip, 1 = no slip on obstacles
    Vector3 gravity{0, -9.81f, 0};
    int maxParticles = 250000;
};

struct ParticleEmitter {
    bool enabled = false;
    Vector3 position;             // nozzle center
    Vector3 direction{1, 0, 0};   // unit
    float radius = 0.05f;
    float speed = 2.0f;
    float accumulated = 0;     // distance travelled since last layer
};

class ParticleSystem {
public:
    ParticleParams params;
    ParticleEmitter emitter;

    void reset(const AABB& domain);
    // Liquid filling the box; returns its particle group.
    int addBlock(const AABB& box, const Vector3& velocity = Vector3(0.0f));
    // Soft body: particles on a lattice (spacing 2r) filling the closed mesh (world space);
    // stiffness 0..1 (1 = rigid). Returns the soft body index, -1 if nothing fitted.
    int addSoftBody(const TriMesh& shape, float density, float stiffness, const Vector3& color,
                    const Vector3& velocity = Vector3(0.0f));
    // Cloth: particle grid from `origin` along the edges u (warp) and v (weft), spacing
    // params.clothSpacing * r. pinMask pins corners: 1 = origin, 2 = origin + u, 4 = origin + v,
    // 8 = origin + u + v; or whole edges: 16 = the edge along u at the origin (e.g. a curtain rod),
    // 32 = the opposite edge, 64 = the edge along v at the origin, 128 = the opposite one.
    // Returns the cloth index.
    int addCloth(const Vector3& origin, const Vector3& u, const Vector3& v, const ClothMaterial& material, int pinMask,
                 const Vector3& color);
    void setStaticMesh(const MeshBVH* bvh) { mesh_ = bvh; }
    void setRigidWorld(RigidWorld* w) { rigid_ = w; }

    // Groups: every add makes one - a block of liquid, a soft body, a cloth - and the emitter keeps
    // one of its own for everything it releases. removeGroup takes one group out without touching
    // the rest (the editor's meta-objects: a role changed on one object): its particles go, the
    // arrays close up in order, and the soft bodies and cloths left over are renumbered. The
    // simulation of everything else goes on as if the group had never been there.
    void removeGroup(int group);
    int groupOf(int particle) const { return group_[size_t(particle)]; }
    int softBodyGroup(int body) const { return softBodies_[size_t(body)].group; }
    int clothGroup(int cloth) const { return cloths_[size_t(cloth)].group; }
    int emitterGroup() const { return emitterGroup_; } // -1 until the emitter released something
    size_t groupSize(int group) const;

    // Advances one substep of length dt.
    void step(float dt);
    // Fire (cloths whose material burns): heating by the gas and the flame's radiation (gas(x)),
    // ignition, burning, charring. The fabric loses mass as it burns. What the cloths give to the
    // gas is appended to `out` (only particles with something to give).
    void burnCloths(float dt, const std::function<GasHeat(const Vector3&)>& gas, float ambientTemperature,
                    std::vector<FireOutput>& out);

    size_t size() const { return x_.size(); }
    const std::vector<Vector3>& positions() const { return x_; }
    const std::vector<Vector3>& velocities() const { return v_; }
    const std::vector<float>& densities() const { return rho_; }
    const std::vector<uint8_t>& phases() const { return phase_; }
    const std::vector<float>& invMasses() const { return invMass_; }
    const std::vector<SoftBody>& softBodies() const { return softBodies_; }
    const std::vector<Cloth>& cloths() const { return cloths_; }
    // Surface of a soft body now (its mesh skinned to the clusters).
    void softBodySurface(size_t body, std::vector<Vector3>& out) const { skinSurface(softBodies_[body], x_, out); }

    // Mouse grab: the particle nearest to `point` (within 3 spacings) and its neighbours of the same
    // object (cloth: that particle alone; soft body / liquid: a small ball) follow the target
    // kinematically; everything else is dragged along by the constraints - a hard yank tears cloth.
    bool grab(const Vector3& point);
    void setGrabTarget(const Vector3& target) { grab_.target = target; }
    void releaseGrab();
    bool grabbing() const { return !grab_.particles.empty(); }
    Vector3 grabAnchor() const { return grabbing() ? x_[grab_.particles[0]] : Vector3(0.0f); }
    Vector3 grabTarget() const { return grab_.target; }
    // Pins every solid particle (soft body, cloth) inside `region` (world point -> true): its inverse
    // mass becomes 0 and it stays put - a beam clamped in a wall, a sheet nailed to a board. Shape
    // matching treats pinned particles as an outside support. Returns how many were pinned.
    int pinParticles(const std::function<bool(const Vector3&)>& region);
    size_t fluidCount() const { return fluidCount_; }
    bool hasSolids() const { return fluidCount_ < x_.size(); } // soft bodies or cloth
    void addVelocity(int i, const Vector3& dv) { v_[i] += dv; } // external forces (e.g. gas drag)
    size_t particleContactCount() const { return contacts_.size(); } // candidate pairs of the last substep
    float kernelRadius() const { return h_; }
    float particleMass() const { return mass_; }
    float spacing() const { return 2.0f * params.particleRadius; }
    const AABB& domain() const { return domain_; }
    float averageDensityError() const { return avgDensityError_; }
    // Part of the kernel of a particle at p lying beyond the domain walls (0 .. 1/2 per wall) and
    // its gradient: the walls count in the density as liquid at rest (see computeLambda).
    float wallVolume(const Vector3& p, Vector3& gradient) const;
    float maxSpeed() const { return maxSpeed_; }

private:
    void emitParticles(float dt);
    void buildGrid(const std::vector<Vector3>& pts);
    void findNeighbors();
    void computeLambda();
    void computeDeltaP();
    // Walls, obstacle mesh, rigid bodies. `start`: where the particle was at the start of the step
    // of length dt that moved it to p (the friction acts on that motion).
    void collide(int i, Vector3& p, const Vector3& start, bool recordImpulse, float dt);
    void collideWallsAndMesh(Vector3& p, const Vector3& start) const; // the first two of collide's three
    void applyViscosityAndVorticity(float dt);
    // A soft body's signed distance field at a point inside it: depth under the surface, outward normal.
    void measureSurface(const MeshBVH& bvh, const Vector3& p, float& depth, Vector3& normal) const;
    // (Colours live on the objects: SoftBody::color, Cloth::color.)
    void addParticle(const Vector3& x, const Vector3& v, ParticlePhase phase, int object, int group, float invMass,
                     float volume = 1.0f);
    // The steps of removeGroup: which particles stay and where they go, then the arrays, the
    // soft bodies and the cloths follow the new numbering.
    std::vector<int> renumberWithout(int group) const;
    void compactParticles(const std::vector<int>& newIndex, size_t kept);
    void renumberSolids(int group, const std::vector<int>& newIndex);
    // Particles of different phases, of different soft bodies and non-adjacent particles of one
    // cloth keep d0 = 2r apart (ParticleContacts.cpp): the candidate pairs are collected once per
    // substep; pairs of bodies found inside each other are pulled apart first (pre-stabilization);
    // then every pair is projected Gauss-Seidel style in the solid passes.
    struct ParticleContact {
        int i, j;
        Vector3 normal;    // the way i is pushed (j the other way), fixed for the substep
        float lift;        // stack mass scaling: i's inverse mass counts x lift, j's / lift
        float target;      // how far apart the main solve keeps them: d0, less for an intersection
        bool intersecting; // its two bodies are inside each other: the pre-stabilization's to undo
    };
    void findParticleContacts();
    float stackLift(int i, int j) const; // FleX's stiff-stack mass scaling (eq. 21)
    // Pre-stabilization: bodies found inside each other at the start of the step pulled apart, x
    // and p alike (no speed), every soft body as a whole (no dent), the way out of the bodies.
    void preStabilizeContacts();
    void findMovers();
    long long moverPair(const ParticleContact& c) const;
    int findIntersections(); // how many contacts belong to an intersection
    void rememberUnresolved();
    Vector3 intersectionNormal(int i, int j) const; // FleX's signed-distance normal (eq. 17, 20)
    bool pushMoversApart();
    void pushMoversOutOfSupports();
    void pushOutOfSupports(Vector3& p);
    void moveByMovers();
    // The main solve: its normals and targets for the step, then one pass per solid iteration.
    void setMainSolveTargets();
    void solveParticleContacts();
    // Cloth dynamics in small steps (gravity, one constraint pass with every thread solved exactly,
    // tearing each) from the start of the substep: "small steps" converge far better than more
    // iterations, so the thread tensions - and with them the tearing - are physical, not solver lag.
    // clothSmallSteps: clothSubsteps of them, more while the threads are pulled hard.
    void stepClothsInSmallSteps(float dt);
    int clothSmallSteps(const Cloth& c, float dt) const;
    std::vector<ParticleContact> contacts_;
    struct ParticleGrab {
        std::vector<int> particles;       // [0] = the picked one
        std::vector<Vector3> offsets;     // from the picked particle at grab time
        std::vector<float> savedInvMass;  // restored on release
        Vector3 target;
    } grab_;
    bool isFluid(int i) const { return phase_[i] == uint8_t(ParticlePhase::Fluid); }

    inline float W(float r2) const {
        if (r2 >= h2_) return 0.0f;
        float d = h2_ - r2;
        return poly6_ * d * d * d;
    }
    inline Vector3 gradW(const Vector3& r) const {
        float l2 = length2(r);
        if (l2 >= h2_ || l2 < 1e-20f) return Vector3(0.0f);
        float l = std::sqrt(l2);
        float d = h_ - l;
        return r * (spikyGrad_ * d * d / l);
    }

    AABB domain_;
    float h_ = 0.06f, h2_ = 0, poly6_ = 0, spikyGrad_ = 0, mass_ = 1, deltaQW_ = 1;

    std::vector<Vector3> x_, v_, p_, dp_, omega_, vtmp_;
    std::vector<float> rho_, lambda_;
    // per particle: phase, object id (-1 = fluid), inverse mass (0 = pinned) and rest position
    // (cloth self-collision filter)
    std::vector<uint8_t> phase_;
    std::vector<int> object_;
    std::vector<int> group_;   // the group of every particle (removeGroup)
    int nextGroup_ = 0, emitterGroup_ = -1;
    std::vector<float> invMass_;
    std::vector<float> volume_; // volume relative to a fluid particle (cloth sheets are thinner)
    std::vector<Vector3> rest_;
    // The signed distance field of a soft body, sampled at its particles (Macklin et al. 2014,
    // sec. 5.1): how deep under its body's surface a particle sits at rest (-1: liquid and cloth
    // have no surface) and the outward normal there, at rest and turned with the clusters now.
    std::vector<float> surfaceDepth_;
    std::vector<Vector3> restSurfaceNormal_, surfaceNormal_;
    bool hasSurface(int i) const { return surfaceDepth_[size_t(i)] >= 0; }
    std::vector<SoftBody> softBodies_;
    std::vector<Cloth> cloths_;
    int nextObject_ = 0;
    size_t fluidCount_ = 0;

    // neighbours
    static constexpr int kMaxNeighbors = 80;
    std::vector<int> nbrCount_;
    std::vector<int> nbr_;
    // grid
    int gx_ = 1, gy_ = 1, gz_ = 1;
    std::vector<int> cellStart_, cellOf_, sorted_;
    // Pre-stabilization's movers: every soft body is one (numbered n + its index), every other
    // particle is its own (numbered as the particle). Per particle its mover, its soft body (-1:
    // none) and its push out of the supports in the current pass; per mover its inverse mass
    // (0: held), how far this step has shifted it, and the move of the current pass with the
    // number of pairs that asked for it. The pairs of movers inside each other (moverPair, sorted):
    // this step's, and those the pre-stabilization has not pulled apart yet (for the next step).
    std::vector<int> moverOf_, softBodyOf_, moverAsks_;
    std::vector<float> moverInvMass_;
    std::vector<Vector3> supportPush_, moverShift_, moverMove_;
    std::vector<long long> intersections_, unresolved_;

    // Two-way coupling with rigid bodies inside the iterations (XPBD contacts, Mueller et al. 2020,
    // "Detailed Rigid Body Simulation with Extended Position Based Dynamics"): every contact splits
    // its correction between particle and body by their generalized inverse masses; the bodies move
    // (bodyShift_, small rotation bodyTurn_) during the substep, so the particles meet them where
    // they are, and the motion becomes their velocity change at the end. The contacts of a body are
    // solved one after another (Gauss-Seidel): a light body hit by a lot of water at once is pushed
    // no farther than the water pushes it, and the pressure all around a floating body leaves
    // exactly its buoyancy.
    std::vector<int> contactBody_;                     // body touched in the current pass, -1 none
    std::vector<Vector3> contactNormal_, contactPoint_; // out of the body; on its surface
    std::vector<float> contactDepth_;                   // penetration when found
    std::vector<Vector3> bodyShift_, bodyTurn_;
    void solveBodyContacts(float dt); // after a collision pass
    // The steps of step() (ParticleSystem.cpp).
    void beginStep(int n);
    void predictPositions(float dt);
    void solveIteration(bool solids, float dt);
    void finishStep(float dt);
    // The research layers of the particles (ParticleDebugDraw.cpp): neighbours of the particle at
    // the probe point, density error, soft-body clusters, cloth tension. Only the layers that are on.
    void drawDebug(float dt) const;
    void drawNeighbours() const;
    void drawDensityError() const;
    void drawSoftClusters() const;
    void drawClothTension(float dt) const;

    const MeshBVH* mesh_ = nullptr;
    RigidWorld* rigid_ = nullptr;
    // Bodies near a particle come from the rigid world's tree, not from a loop over all bodies:
    // before a collision pass the tree is refreshed and the query reach is set to the particle
    // radius plus the farthest any body has been shifted or turned in this substep (so the
    // candidates are a superset of what the exact test accepts - the contacts are the same as
    // with the full loop). One candidate vector per worker thread: the passes allocate nothing.
    void prepareBodyQuery(bool coupled);
    float bodyReach_ = 0;
    std::vector<std::vector<int>> bodyCandidates_;
    float avgDensityError_ = 0, maxSpeed_ = 0;
};

} // namespace rf
