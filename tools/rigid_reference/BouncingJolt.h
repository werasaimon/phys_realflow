#pragma once
// Jolt adapter for the shared elastic sphere scenario, using discrete contacts and no stabilization.
#include "BouncingCase.h"
#include "JoltRuntime.h"
#include <Jolt/Core/TempAllocator.h>
#include <Jolt/Core/JobSystemSingleThreaded.h>
#include <Jolt/Physics/PhysicsSystem.h>
#include <Jolt/Physics/Body/BodyCreationSettings.h>
#include <Jolt/Physics/Body/MotionProperties.h>
#include <Jolt/Physics/Collision/Shape/BoxShape.h>
#include <Jolt/Physics/Collision/Shape/SphereShape.h>
#include <Jolt/Physics/Collision/ObjectLayerPairFilterTable.h>
#include <Jolt/Physics/Collision/BroadPhase/BroadPhaseLayerInterfaceTable.h>
#include <Jolt/Physics/Collision/BroadPhase/ObjectVsBroadPhaseLayerFilterTable.h>
#include <array>
#include <memory>
#include <stdexcept>

namespace rf::benchmark {
class BouncingJolt {
public:
    BouncingJolt(int iterations, bool) {
        layers_.MapObjectToBroadPhaseLayer(1, JPH::BroadPhaseLayer(1));
        pairs_.EnableCollision(0, 1);
        pairs_.EnableCollision(1, 1);
        filter_ = std::make_unique<JPH::ObjectVsBroadPhaseLayerFilterTable>(layers_, 2, pairs_, 2);
        world_.Init(128, 0, 1024, 1024, layers_, *filter_, pairs_);
        world_.SetGravity({0,-gravity,0});
        auto settings = world_.GetPhysicsSettings();
        settings.mNumVelocitySteps = JPH::uint(iterations);
        settings.mNumPositionSteps = 0;
        settings.mBaumgarte = settings.mPenetrationSlop = 0;
        settings.mSpeculativeContactDistance = 0.01f;
        settings.mMinVelocityForRestitution = 1e-6f; // API requires a positive threshold.
        settings.mAllowSleeping = false;
        world_.SetPhysicsSettings(settings);
        JPH::RefConst<JPH::Shape> floor = new JPH::BoxShape({25,5,25}, 0);
        add(floor, {0,-5,0}, false);
        JPH::RefConst<JPH::Shape> sphere = new JPH::SphereShape(radius);
        for (int i = 0; i < count; ++i) ids_[size_t(i)] = add(sphere, startPosition(i), true);
    }

    bool step(float dt, Work& work) {
        const auto result = world_.Update(dt, 1, &temporary_, &jobs_);
        if (result != JPH::EPhysicsUpdateError::None) return false;
        ++work.accepted;
        work.time += dt;
        return true;
    }

    std::array<State, count> states() const {
        std::array<State, count> out;
        const auto& api = world_.GetBodyInterface();
        for (int i = 0; i < count; ++i) {
            const auto id = ids_[size_t(i)];
            const auto p = api.GetPosition(id);
            const auto v = api.GetLinearVelocity(id), w = api.GetAngularVelocity(id);
            out[size_t(i)] = {{float(p.GetX()),float(p.GetY()),float(p.GetZ())},
                              {v.GetX(),v.GetY(),v.GetZ()}, {w.GetX(),w.GetY(),w.GetZ()}};
        }
        return out;
    }

private:
    JPH::BodyID add(JPH::RefConst<JPH::Shape> shape, const rf::Vector3& pos, bool dynamic) {
        JPH::BodyCreationSettings settings(shape, {pos.x,pos.y,pos.z}, JPH::Quat::sIdentity(),
            dynamic ? JPH::EMotionType::Dynamic : JPH::EMotionType::Static, dynamic ? 1 : 0);
        settings.mAllowSleeping = false;
        settings.mMotionQuality = JPH::EMotionQuality::Discrete;
        settings.mLinearDamping = settings.mAngularDamping = settings.mFriction = 0;
        settings.mRestitution = 1;
        if (dynamic) {
            settings.mOverrideMassProperties = JPH::EOverrideMassProperties::MassAndInertiaProvided;
            settings.mMassPropertiesOverride.mMass = mass;
            settings.mMassPropertiesOverride.mInertia = JPH::Mat44::sScale({0.04f,0.04f,0.04f});
        }
        const auto id = world_.GetBodyInterface().CreateAndAddBody(settings, JPH::EActivation::Activate);
        if (id.IsInvalid()) throw std::runtime_error("Jolt body creation failed");
        if (dynamic) {
            {
                JPH::BodyLockRead lock(world_.GetBodyLockInterface(), id);
                if (!lock.Succeeded() || std::fabs(lock.GetBody().GetMotionProperties()->GetInverseMass() - 0.1f) > 1e-7f)
                    throw std::runtime_error("Jolt inverse mass mismatch");
            }
            const auto inverse = world_.GetBodyInterface().GetInverseInertia(id);
            for (const auto axis : {JPH::Vec3::sAxisX(), JPH::Vec3::sAxisY(), JPH::Vec3::sAxisZ()})
                if ((inverse.Multiply3x3(axis) - 25 * axis).Length() > 1e-4f)
                    throw std::runtime_error("Jolt inverse inertia mismatch");
        }
        return id;
    }

    JoltRuntime runtime_; // First constructed, last destroyed.
    JPH::BroadPhaseLayerInterfaceTable layers_{2,2};
    JPH::ObjectLayerPairFilterTable pairs_{2};
    std::unique_ptr<JPH::ObjectVsBroadPhaseLayerFilterTable> filter_;
    JPH::PhysicsSystem world_;
    JPH::TempAllocatorImpl temporary_{16 * 1024 * 1024};
    JPH::JobSystemSingleThreaded jobs_{2048};
    std::array<JPH::BodyID, count> ids_;
};
} // namespace rf::benchmark
