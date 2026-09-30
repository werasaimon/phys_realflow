// Optional Jolt v5.3.0 reference for sample 39; built by tools/rigid_reference/CMakeLists.txt.
// Usage: rf_jolt_reference NEW_DIRECTORY [iterations=6] [frames=600] [slop=0.004] [ccd=1]
// Public Jolt APIs only. No Jolt code or dependency is added to rfcore.
#include <Jolt/Jolt.h>
#include <Jolt/RegisterTypes.h>
#include <Jolt/Core/Factory.h>
#include <Jolt/Core/TempAllocator.h>
#include <Jolt/Core/JobSystemSingleThreaded.h>
#include <Jolt/Physics/PhysicsSystem.h>
#include <Jolt/Physics/Body/BodyCreationSettings.h>
#include <Jolt/Physics/Body/MotionProperties.h>
#include <Jolt/Physics/Collision/Shape/ConvexHullShape.h>
#include <Jolt/Physics/Collision/Shape/StaticCompoundShape.h>
#include <Jolt/Physics/Collision/Shape/OffsetCenterOfMassShape.h>
#include <Jolt/Physics/Collision/Shape/PlaneShape.h>
#include <Jolt/Physics/Collision/ObjectLayerPairFilterTable.h>
#include <Jolt/Physics/Collision/BroadPhase/BroadPhaseLayerInterfaceTable.h>
#include <Jolt/Physics/Collision/BroadPhase/ObjectVsBroadPhaseLayerFilterTable.h>
#include "samples/TorusChainsScene.h"
#include "samples/Samples.h"
#include "core/Format.h"
#include "rigid_reference/ReferenceTiming.h"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <unordered_map>

namespace {
JPH::Vec3 vector(const rf::Vector3& v) { return {v.x,v.y,v.z}; }
rf::Vector3 vector(JPH::Vec3Arg v) { return {v.GetX(),v.GetY(),v.GetZ()}; }
JPH::Quat quaternion(const rf::Quaternion& q) { return {q.x,q.y,q.z,q.w}; }

struct Registration {
    Registration() {
        JPH::RegisterDefaultAllocator();
        JPH::Factory::sInstance = new JPH::Factory();
        JPH::RegisterTypes();
    }
    ~Registration() { JPH::UnregisterTypes(); delete JPH::Factory::sInstance; JPH::Factory::sInstance = nullptr; }
};

struct Contacts : JPH::ContactListener {
    float deepest = 0;
    void OnContactAdded(const JPH::Body&, const JPH::Body&, const JPH::ContactManifold& m, JPH::ContactSettings&) override {
        deepest = std::max(deepest, m.mPenetrationDepth);
    }
    void OnContactPersisted(const JPH::Body&, const JPH::Body&, const JPH::ContactManifold& m, JPH::ContactSettings&) override {
        deepest = std::max(deepest, m.mPenetrationDepth);
    }
};

struct Reference {
    JPH::BroadPhaseLayerInterfaceTable layers{2,2};
    JPH::ObjectLayerPairFilterTable pairs{2};
    std::unique_ptr<JPH::ObjectVsBroadPhaseLayerFilterTable> broad;
    Contacts contacts;
    JPH::PhysicsSystem world;
    JPH::TempAllocatorImpl temporary{32*1024*1024};
    JPH::JobSystemSingleThreaded jobs{2048};
    std::vector<JPH::BodyID> bodies;
    std::unordered_map<const rf::CompoundShape*, JPH::RefConst<JPH::Shape>> imported;
    float supportError = 0;
    bool ccd = true;

    Reference() {
        layers.MapObjectToBroadPhaseLayer(1,JPH::BroadPhaseLayer(1));
        pairs.EnableCollision(0,1);
        pairs.EnableCollision(1,1);
    }

    void validateHull(const rf::ConvexHullShape& source, const JPH::ConvexHullShape& target) {
        const auto& mesh = *source.mesh();
        for (size_t i = 0; i < mesh.triangles.size(); ++i) {
            const rf::Vector3 n = mesh.faceNormal(i);
            float support = -rf::kInf;
            for (JPH::uint k = 0; k < target.GetNumPoints(); ++k)
                support = std::max(support, rf::dot(n, vector(target.GetPoint(k) + target.GetCenterOfMass())));
            supportError = std::max(supportError, std::fabs(support - rf::dot(n, source.support(n))));
        }
        if (supportError > 5e-6f) throw std::runtime_error("Jolt hull support differs by more than 5 micrometers");
    }

    JPH::RefConst<JPH::Shape> compound(const rf::CompoundShape& source) {
        const auto found = imported.find(&source);
        if (found != imported.end()) return found->second;
        JPH::StaticCompoundShapeSettings compound;
        JPH::uint32 childIndex = 0;
        for (const auto& child : source.children()) {
            const auto& hull = static_cast<const rf::ConvexHullShape&>(*child.shape);
            JPH::ConvexHullShapeSettings settings;
            settings.mMaxConvexRadius = 0;
            settings.mHullTolerance = 1e-5f;
            for (const auto& v : hull.vertices()) settings.mPoints.push_back(vector(v));
            auto result = settings.Create();
            if (result.HasError()) throw std::runtime_error(result.GetError().c_str());
            validateHull(hull, static_cast<const JPH::ConvexHullShape&>(*result.Get()));
            compound.AddShape(vector(child.t), quaternion(rf::Quaternion::fromMatrix3x3(child.R)), result.Get(),childIndex++);
        }
        auto result = compound.Create();
        if (result.HasError()) throw std::runtime_error(result.GetError().c_str());
        const auto& built = static_cast<const JPH::CompoundShape&>(*result.Get());
        for (const auto& sub : built.GetSubShapes()) {
            const auto& child = source.children()[sub.mUserData];
            const auto& expected = static_cast<const rf::ConvexHullShape&>(*child.shape);
            const auto& actual = static_cast<const JPH::ConvexHullShape&>(*sub.mShape);
            for (size_t t = 0; t < expected.mesh()->triangles.size(); ++t) {
                const auto local = expected.mesh()->faceNormal(t);
                const auto normal = child.R * local;
                float support = -rf::kInf;
                for (JPH::uint k = 0; k < actual.GetNumPoints(); ++k) {
                    const auto point = built.GetCenterOfMass() + sub.GetPositionCOM() + sub.GetRotation()*actual.GetPoint(k);
                    support = std::max(support,rf::dot(normal,vector(point)));
                }
                supportError = std::max(supportError,std::fabs(support-rf::dot(normal,child.t+child.R*expected.support(local))));
            }
        }
        if (supportError > 5e-6f) throw std::runtime_error("Jolt compound transform import exceeds 5 micrometers");
        JPH::RefConst<JPH::Shape> shape = new JPH::OffsetCenterOfMassShape(result.Get(), -result.Get()->GetCenterOfMass());
        imported.emplace(&source, shape);
        return shape;
    }

    void load(const rf::RigidWorld& source, int iterations, float slop) {
        // Construct the broad-phase filter after its source tables have been configured.
        broad = std::make_unique<JPH::ObjectVsBroadPhaseLayerFilterTable>(layers,2,pairs,2);
        world.Init(1024,0,8192,32768,layers,*broad,pairs);
        world.SetGravity(vector(source.params.gravity));
        auto physics = world.GetPhysicsSettings();
        physics.mNumVelocitySteps = JPH::uint(iterations);
        physics.mPenetrationSlop = slop;
        physics.mSpeculativeContactDistance = source.params.contactMargin;
        world.SetPhysicsSettings(physics);
        world.SetContactListener(&contacts);
        auto& api = world.GetBodyInterface();
        for (const auto& body : source.bodies()) {
            auto shape = compound(static_cast<const rf::CompoundShape&>(*body.shape));
            JPH::BodyCreationSettings settings(shape, vector(body.pos), quaternion(body.rot),
                body.mass > 0 ? JPH::EMotionType::Dynamic : JPH::EMotionType::Static, body.mass > 0 ? 1 : 0);
            settings.mAllowSleeping = false;
            settings.mMotionQuality = ccd ? JPH::EMotionQuality::LinearCast : JPH::EMotionQuality::Discrete;
            settings.mLinearVelocity = vector(body.vel);
            settings.mAngularVelocity = vector(body.angVel);
            settings.mLinearDamping = source.params.linearDamping;
            settings.mAngularDamping = source.params.angularDamping;
            settings.mFriction = body.friction;
            settings.mRestitution = body.restitution;
            settings.mApplyGyroscopicForce = true;
            if (body.mass > 0) {
                settings.mOverrideMassProperties = JPH::EOverrideMassProperties::MassAndInertiaProvided;
                settings.mMassPropertiesOverride.mMass = body.mass;
                settings.mMassPropertiesOverride.mInertia = JPH::Mat44::sScale(
                    {1/body.invInertiaLocal.x,1/body.invInertiaLocal.y,1/body.invInertiaLocal.z});
            }
            const auto id = api.CreateAndAddBody(settings, JPH::EActivation::Activate);
            if (id.IsInvalid()) throw std::runtime_error("Jolt body capacity exceeded");
            if (body.mass > 0) {
                JPH::BodyLockRead lock(world.GetBodyLockInterface(),id);
                if (!lock.Succeeded() || std::fabs(lock.GetBody().GetMotionProperties()->GetInverseMass()-body.invMass)
                    > 1e-5f*std::max(body.invMass,1.0f)) throw std::runtime_error("Jolt mass import mismatch");
            }
            if (body.mass > 0) for (const rf::Vector3 axis : {rf::Vector3(1,0,0),rf::Vector3(0,1,0),rf::Vector3(0,0,1)}) {
                const auto expected = body.applyInvInertiaWorld(axis);
                const auto actual = vector(api.GetInverseInertia(id).Multiply3x3(vector(axis)));
                if (rf::length(expected-actual) > 1e-4f*rf::length(expected))
                    throw std::runtime_error("Jolt inertia import mismatch");
            }
            bodies.push_back(id);
        }
        const auto box = source.domain();
        const rf::Vector3 normals[] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
        const rf::Vector3 points[] = {box.lo,box.hi,box.lo,box.hi,box.lo,box.hi};
        for (int k = 0; k < 6; ++k) {
            JPH::RefConst<JPH::Shape> plane = new JPH::PlaneShape(JPH::Plane(vector(normals[k]), -rf::dot(normals[k], points[k])));
            JPH::BodyCreationSettings settings(plane,JPH::RVec3::sZero(),JPH::Quat::sIdentity(),JPH::EMotionType::Static,0);
            settings.mFriction = 0.6f;
            settings.mRestitution = 0;
            const auto id = api.CreateAndAddBody(settings,JPH::EActivation::DontActivate);
            if (id.IsInvalid()) throw std::runtime_error("Jolt wall creation failed");
        }
    }

    void copyPoses(rf::RigidWorld& target) const {
        const auto& api = world.GetBodyInterface();
        for (size_t i = 0; i < bodies.size(); ++i) {
            auto& body = target.bodies()[i];
            body.pos = vector(api.GetPosition(bodies[i]));
            const auto q = api.GetRotation(bodies[i]);
            body.rot = {q.GetW(),q.GetX(),q.GetY(),q.GetZ()};
            body.vel = vector(api.GetLinearVelocity(bodies[i]));
            body.angVel = vector(api.GetAngularVelocity(bodies[i]));
            if (!std::isfinite(rf::length2(body.pos)) || !std::isfinite(rf::length2(body.vel)))
                throw std::runtime_error("non-finite Jolt state");
        }
    }
};
}

int main(int argc, char** argv) {
    try {
        if (argc < 2 || argc > 6) throw std::runtime_error("NEW_DIRECTORY [iterations] [frames] [slop] [ccd]");
        const int iterations = argc > 2 ? std::stoi(argv[2]) : 6;
        const int frames = argc > 3 ? std::stoi(argv[3]) : 600;
        float slop = 0.004f;
        if (argc > 4 && !rf::parseNumber(argv[4],slop)) throw std::runtime_error("invalid slop");
        if (iterations < 1 || frames < 1 || !std::isfinite(slop) || slop < 0) throw std::runtime_error("invalid settings");
        const std::filesystem::path output(argv[1]);
        if (!std::filesystem::create_directory(output)) throw std::runtime_error("output must be a new directory");
        Registration registration;
        rf::Simulation sim;
        rf::loadSample(sim,rf::Preset::TorusChains);
        auto* scene = static_cast<rf::TorusChainsScene*>(sim.scene());
        Reference reference;
        reference.ccd = argc <= 5 || std::stoi(argv[5]) != 0;
        reference.load(sim.rigid,iterations,slop);
        const auto initial = sim.rigid.bodies();
        reference.copyPoses(sim.rigid);
        for (size_t i = 0; i < initial.size(); ++i) {
            const auto& a = initial[i];
            const auto& b = sim.rigid.bodies()[i];
            if (rf::length(a.pos-b.pos) > 1e-6f || rf::length(a.vel-b.vel) > 1e-6f ||
                rf::length(a.rot.rotate({1,2,3})-b.rot.rotate({1,2,3})) > 1e-5f)
                throw std::runtime_error("Jolt pose/velocity import roundtrip failed");
        }
        scene->afterStep(sim);
        if (scene->brokenPairs()) throw std::runtime_error("import lost linking");
        std::ofstream csv(output / "frames.csv"), summary(output / "summary.txt");
        if (!csv || !summary) throw std::runtime_error("cannot write report");
        csv << "frame,time_s,broken_pairs,deepest_jolt_contact_m,engine_ms,observer_ms\n";
        int first = -1, worst = 0;
        for (int frame = 1; frame <= frames; ++frame) {
            ReferenceTiming timing;
            for (int sub = 0; sub < 10; ++sub) {
                timing.start();
                if (reference.world.Update(1.0f/600,1,&reference.temporary,&reference.jobs) != JPH::EPhysicsUpdateError::None)
                    throw std::runtime_error("Jolt update exceeded capacity");
                timing.finishEngine();
                reference.copyPoses(sim.rigid);
                scene->afterStep(sim);
                worst = std::max(worst,scene->brokenPairs());
                if (scene->brokenPairs() && first < 0) first = frame;
                timing.finishObserver();
            }
            csv << frame << ',' << rf::numberText(double(frame)/60) << ',' << scene->brokenPairs()
                << ',' << rf::numberText(reference.contacts.deepest) << timing.csv() << '\n';
            csv.flush();
            if (frame % 30 == 0) std::cout << "frame " << frame << ", worst broken " << worst << std::endl;
        }
        summary << "Jolt v5.3.0\niterations " << iterations << "\nframes " << frames << "\nsubsteps 10\nslop " << rf::numberText(slop)
                << "\nlinear_cast " << reference.ccd << "\nfirst_break_frame " << first << "\nworst_broken " << worst
                << "\nmax_reported_contact_depth_m " << rf::numberText(reference.contacts.deepest)
                << "\nhull_face_support_error_m " << rf::numberText(reference.supportError)
                << "\nconvex radius 0; hull tolerance 1e-5; input masses/inertias provided; COM offset compensated\n"
                << "sleep disabled; gyro enabled; 2 position iterations; other settings Jolt defaults\n"
                << "single thread; manifold depth is not a swept-geometry certificate\n";
        std::cout << "first break " << first << ", worst broken " << worst << ", contact depth "
                  << rf::numberText(reference.contacts.deepest) << " m\n";
        return worst ? 1 : 0;
    } catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 2; }
}
