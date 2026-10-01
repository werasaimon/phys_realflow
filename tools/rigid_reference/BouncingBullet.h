#pragma once
// Bullet discrete adapter for the same bodies; no hidden substeps or deactivation.
#include "BouncingCase.h"
#include <btBulletDynamicsCommon.h>
#include <array>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

namespace rf::benchmark {
rf::Vector3 vector(const btVector3& value) {
    return {float(value.x()), float(value.y()), float(value.z())};
}

class BouncingBullet {
public:
    explicit BouncingBullet(int iterations, bool) {
        world_.setGravity({0,-gravity,0});
        auto& info = world_.getSolverInfo();
        info.m_numIterations = iterations;
        info.m_erp = info.m_erp2 = info.m_globalCfm = info.m_linearSlop = 0;
        info.m_frictionERP = info.m_splitImpulseTurnErp = 0;
        info.m_restitutionVelocityThreshold = 0;
        shapes_.push_back(std::make_unique<btBoxShape>(btVector3(25,5,25)));
        shapes_.back()->setMargin(0);
        add(shapes_.back().get(), {0,-5,0}, 0);
        shapes_.push_back(std::make_unique<btSphereShape>(radius));
        for (int i = 0; i < count; ++i) add(shapes_.back().get(), startPosition(i), mass);
    }

    ~BouncingBullet() {
        for (const auto& body : bodies_) world_.removeRigidBody(body.get());
    }

    bool step(float dt, Work& work) {
        world_.stepSimulation(dt, 0); // One explicit step; no hidden time accumulator/substeps.
        ++work.accepted;
        work.time += dt;
        return true;
    }

    std::array<State, count> states() const {
        std::array<State, count> out;
        for (int i = 0; i < count; ++i) {
            const auto& body = *bodies_[size_t(i + 1)];
            out[size_t(i)] = {vector(body.getWorldTransform().getOrigin()), vector(body.getLinearVelocity()),
                              vector(body.getAngularVelocity())};
        }
        return out;
    }

private:
    void add(btCollisionShape* shape, const rf::Vector3& position, float bodyMass) {
        const btVector3 inertia = bodyMass > 0 ? btVector3(0.04f,0.04f,0.04f) : btVector3(0,0,0);
        btRigidBody::btRigidBodyConstructionInfo info(bodyMass, nullptr, shape, inertia);
        info.m_startWorldTransform.setIdentity();
        info.m_startWorldTransform.setOrigin({position.x,position.y,position.z});
        info.m_friction = info.m_linearDamping = info.m_angularDamping = 0;
        info.m_restitution = 1;
        auto body = std::make_unique<btRigidBody>(info);
        if (bodyMass > 0 && (std::fabs(float(body->getInvMass()) - 0.1f) > 1e-7f
            || (body->getInvInertiaDiagLocal() - btVector3(25,25,25)).length() > 1e-4f))
            throw std::runtime_error("Bullet mass/inertia mismatch");
        body->setActivationState(DISABLE_DEACTIVATION);
        world_.addRigidBody(body.get());
        bodies_.push_back(std::move(body));
    }

    btDefaultCollisionConfiguration configuration_;
    btCollisionDispatcher dispatcher_{&configuration_};
    btDbvtBroadphase broadphase_;
    btSequentialImpulseConstraintSolver solver_;
    btDiscreteDynamicsWorld world_{&dispatcher_, &broadphase_, &solver_, &configuration_};
    std::vector<std::unique_ptr<btCollisionShape>> shapes_;
    std::vector<std::unique_ptr<btRigidBody>> bodies_;
};

} // namespace rf::benchmark
