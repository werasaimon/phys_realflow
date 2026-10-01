#pragma once
// Box3D adapter: hertz=0 zeros overlap softness terms; the native relaxation pass still runs.
// Box3D exposes restitution iterations separately; its normal solver passes remain native.
#include "BouncingCase.h"
#include <box3d/box3d.h>
#include <box3d/collision.h>
#include <array>
#include <cmath>
#include <stdexcept>

namespace rf::benchmark {
class BouncingBox3d {
public:
    BouncingBox3d(int iterations, bool nativeContact) {
        auto settings = b3DefaultWorldDef();
        settings.gravity = {0,-gravity,0};
        settings.enableSleep = settings.enableContinuous = false;
        settings.workerCount = 1;
        settings.restitutionThreshold = 0;
        if (!nativeContact) settings.contactHertz = settings.contactSpeed = 0;
        settings.restitutionIterations = iterations;
        world_ = b3CreateWorld(&settings);
        if (B3_IS_NULL(world_)) throw std::runtime_error("Box3D world creation failed");
        try {
            const auto floor = add({0,-5,0}, false);
            auto shape = shapeSettings();
            const auto hull = b3MakeBoxHull(25,5,25);
            if (B3_IS_NULL(b3CreateHullShape(floor, &shape, &hull.base))) throw std::runtime_error("Box3D floor failed");
            for (int i = 0; i < count; ++i) {
                const auto id = add(startPosition(i), true);
                b3Sphere sphere{};
                sphere.radius = radius;
                if (B3_IS_NULL(b3CreateSphereShape(id, &shape, &sphere))) throw std::runtime_error("Box3D sphere failed");
                b3MassData properties{};
                properties.mass = mass;
                properties.inertia = {{0.04f,0,0},{0,0.04f,0},{0,0,0.04f}};
                b3Body_SetMassData(id, properties);
                const auto actual = b3Body_GetMassData(id);
                if (actual.mass != mass || std::fabs(actual.inertia.cx.x - 0.04f) > 1e-7f
                    || std::fabs(actual.inertia.cy.y - 0.04f) > 1e-7f || std::fabs(actual.inertia.cz.z - 0.04f) > 1e-7f)
                    throw std::runtime_error("Box3D mass/inertia mismatch");
                ids_[size_t(i)] = id;
            }
        } catch (...) { b3DestroyWorld(world_); throw; }
    }

    ~BouncingBox3d() { b3DestroyWorld(world_); }
    BouncingBox3d(const BouncingBox3d&) = delete;
    BouncingBox3d& operator=(const BouncingBox3d&) = delete;

    bool step(float dt, Work& work) {
        b3World_Step(world_, dt, 1);
        ++work.accepted;
        work.time += dt;
        return true;
    }

    std::array<State, count> states() const {
        std::array<State, count> out;
        for (int i = 0; i < count; ++i) {
            const auto id = ids_[size_t(i)];
            const auto p = b3Body_GetPosition(id), v = b3Body_GetLinearVelocity(id), w = b3Body_GetAngularVelocity(id);
            out[size_t(i)] = {{p.x,p.y,p.z}, {v.x,v.y,v.z}, {w.x,w.y,w.z}};
        }
        return out;
    }

private:
    static b3ShapeDef shapeSettings() {
        auto settings = b3DefaultShapeDef();
        settings.baseMaterial.friction = 0;
        settings.baseMaterial.restitution = 1;
        settings.updateBodyMass = false;
        return settings;
    }

    b3BodyId add(const rf::Vector3& position, bool dynamic) {
        auto settings = b3DefaultBodyDef();
        settings.type = dynamic ? b3_dynamicBody : b3_staticBody;
        settings.position = {position.x,position.y,position.z};
        settings.linearDamping = settings.angularDamping = 0;
        settings.enableSleep = false;
        const auto id = b3CreateBody(world_, &settings);
        if (B3_IS_NULL(id)) throw std::runtime_error("Box3D body creation failed");
        return id;
    }

    b3WorldId world_ = b3_nullWorldId;
    std::array<b3BodyId, count> ids_;
};
} // namespace rf::benchmark
