// Optional Bullet reference for the exact convex children, masses and poses of sample 39.
// Build separately; Bullet is never linked into rfcore:
// c++ -std=c++17 -O3 -Isrc -I. -I/usr/include/bullet tools/BulletTorusReference.cpp
//   build-core/librfsamples.a build-core/librfcore.a -lBulletDynamics -lBulletCollision
//   -lLinearMath -pthread -o /tmp/rf_bullet_reference
// Usage: /tmp/rf_bullet_reference NEW_DIRECTORY [iterations=6] [frames=600] [polyhedral=0]
#include "samples/TorusChainsScene.h"
#include "samples/Samples.h"
#include "core/Format.h"
#include <btBulletDynamicsCommon.h>
#include "rigid_reference/ReferenceTiming.h"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <unordered_map>

namespace {
btVector3 vector(const rf::Vector3& v) { return {v.x, v.y, v.z}; }
rf::Vector3 vector(const btVector3& v) { return {float(v.x()), float(v.y()), float(v.z())}; }
btQuaternion quaternion(const rf::Quaternion& q) { return {q.x, q.y, q.z, q.w}; }

struct Reference {
    btDefaultCollisionConfiguration config;
    btCollisionDispatcher dispatcher{&config};
    btDbvtBroadphase broadphase;
    btSequentialImpulseConstraintSolver solver;
    btDiscreteDynamicsWorld world{&dispatcher, &broadphase, &solver, &config};
    std::vector<std::unique_ptr<btCollisionShape>> shapes;
    std::vector<std::unique_ptr<btRigidBody>> bodies;
    std::unordered_map<const rf::CompoundShape*, btCollisionShape*> imported;
    bool polyhedral = false;

    ~Reference() { for (const auto& body : bodies) world.removeRigidBody(body.get()); }

    btCollisionShape* compound(const rf::CompoundShape& source) {
        const auto existing = imported.find(&source);
        if (existing != imported.end()) return existing->second;
        auto result = std::make_unique<btCompoundShape>();
        for (const auto& child : source.children()) {
            const auto* hull = dynamic_cast<const rf::ConvexHullShape*>(child.shape.get());
            if (!hull) throw std::runtime_error("reference requires convex hull children");
            auto shape = std::make_unique<btConvexHullShape>();
            shape->setMargin(0); // preserve the input hull; Bullet's default 4 cm changes this scene
            for (const auto& v : hull->vertices()) shape->addPoint(vector(v), false);
            shape->recalcLocalAabb();
            if (polyhedral && !shape->initializePolyhedralFeatures()) throw std::runtime_error("Bullet polyhedral initialization failed");
            const btTransform pose(quaternion(rf::Quaternion::fromMatrix3x3(child.R)), vector(child.t));
            result->addChildShape(pose, shape.get());
            shapes.push_back(std::move(shape));
        }
        auto* ptr = result.get();
        shapes.push_back(std::move(result));
        imported.emplace(&source, ptr);
        return ptr;
    }

    void add(btCollisionShape* shape, const btTransform& pose, float mass, const btVector3& inertia,
             float friction, float restitution) {
        btRigidBody::btRigidBodyConstructionInfo info(mass, nullptr, shape, inertia);
        info.m_startWorldTransform = pose;
        info.m_friction = friction;
        info.m_restitution = restitution;
        info.m_linearDamping = 0.02f;
        info.m_angularDamping = 0.05f;
        auto body = std::make_unique<btRigidBody>(info);
        body->setActivationState(DISABLE_DEACTIVATION);
        world.addRigidBody(body.get());
        bodies.push_back(std::move(body));
    }

    void load(const rf::RigidWorld& source) {
        world.setGravity(vector(source.params.gravity));
        for (const auto& body : source.bodies()) {
            const auto& shape = static_cast<const rf::CompoundShape&>(*body.shape);
            rf::Vector3 inertia(0);
            if (body.mass > 0) inertia = {1/body.invInertiaLocal.x, 1/body.invInertiaLocal.y, 1/body.invInertiaLocal.z};
            add(compound(shape), btTransform(quaternion(body.rot), vector(body.pos)), body.mass, vector(inertia),
                std::sqrt(body.friction), std::sqrt(body.restitution));
            bodies.back()->setLinearVelocity(vector(body.vel));
            bodies.back()->setAngularVelocity(vector(body.angVel));
        }
        const auto& box = source.domain();
        const rf::Vector3 normals[] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
        const rf::Vector3 points[] = {box.lo,box.hi,box.lo,box.hi,box.lo,box.hi};
        for (int k = 0; k < 6; ++k) {
            auto plane = std::make_unique<btStaticPlaneShape>(vector(normals[k]), rf::dot(normals[k], points[k]));
            add(plane.get(), btTransform::getIdentity(), 0, btVector3(0,0,0), std::sqrt(0.6f), std::sqrt(0.2f));
            shapes.push_back(std::move(plane));
        }
    }

    void copyPoses(rf::RigidWorld& target) const {
        for (size_t i = 0; i < target.bodies().size(); ++i) {
            auto& body = target.bodies()[i];
            const auto& pose = bodies[i]->getWorldTransform();
            const auto q = pose.getRotation();
            body.pos = vector(pose.getOrigin());
            body.rot = {float(q.w()),float(q.x()),float(q.y()),float(q.z())};
            body.vel = vector(bodies[i]->getLinearVelocity());
            body.angVel = vector(bodies[i]->getAngularVelocity());
        }
    }

    float validateGeometry(const rf::RigidWorld& source) const {
        float worst = 0;
        for (size_t i = 0; i < source.bodies().size(); ++i) {
            const auto& body = source.bodies()[i];
            if (std::fabs(float(bodies[i]->getInvMass())-body.invMass) > 1e-5f*std::max(body.invMass,1.0f))
                throw std::runtime_error("Bullet inverse mass import mismatch");
            if (body.mass > 0) for (const rf::Vector3 axis : {rf::Vector3(1,0,0),rf::Vector3(0,1,0),rf::Vector3(0,0,1)}) {
                const auto expected = body.applyInvInertiaWorld(axis);
                const auto actual = vector(bodies[i]->getInvInertiaTensorWorld()*vector(axis));
                if (rf::length(expected-actual) > 1e-4f*rf::length(expected))
                    throw std::runtime_error("Bullet inverse inertia import mismatch");
            }
            const auto& parts = static_cast<const rf::CompoundShape&>(*body.shape).children();
            const auto* shape = static_cast<const btCompoundShape*>(bodies[i]->getCollisionShape());
            if (shape->getNumChildShapes() != int(parts.size())) throw std::runtime_error("child count mismatch");
            for (size_t j = 0; j < parts.size(); ++j) {
                const auto& part = parts[j];
                const auto* hull = static_cast<const btConvexHullShape*>(shape->getChildShape(int(j)));
                const auto& vertices = static_cast<const rf::ConvexHullShape&>(*part.shape).vertices();
                if (hull->getNumPoints() != int(vertices.size()) || hull->getMargin() != 0)
                    throw std::runtime_error("hull geometry mismatch");
                const auto pose = bodies[i]->getWorldTransform() * shape->getChildTransform(int(j));
                for (size_t k = 0; k < vertices.size(); ++k) {
                    const auto expected = body.pos + body.rotation() * (part.t + part.R * vertices[k]);
                    worst = std::max(worst, rf::length(expected - vector(pose * hull->getUnscaledPoints()[k])));
                }
            }
        }
        if (worst > 5e-6f) throw std::runtime_error("world geometry import exceeds 5 micrometers");
        return worst;
    }
};
}

int main(int argc, char** argv) {
    try {
        if (argc < 2 || argc > 5) throw std::runtime_error("NEW_DIRECTORY [iterations] [frames] [polyhedral]");
        const int iterations = argc > 2 ? std::stoi(argv[2]) : 6;
        const int frames = argc > 3 ? std::stoi(argv[3]) : 600;
        if (iterations < 1 || frames < 1) throw std::runtime_error("positive iterations/frames required");
        const std::filesystem::path output(argv[1]);
        if (!std::filesystem::create_directory(output)) throw std::runtime_error("output must be a new directory");
        rf::Simulation sim;
        rf::loadSample(sim, rf::Preset::TorusChains);
        auto* scene = static_cast<rf::TorusChainsScene*>(sim.scene());
        Reference reference;
        reference.polyhedral = argc > 4 && std::stoi(argv[4]) != 0;
        reference.load(sim.rigid);
        const float importError = reference.validateGeometry(sim.rigid);
        const auto initial = sim.rigid.bodies();
        reference.copyPoses(sim.rigid);
        for (size_t i = 0; i < initial.size(); ++i) {
            const auto& a = initial[i];
            const auto& b = sim.rigid.bodies()[i];
            if (rf::length(a.pos - b.pos) > 1e-6f || rf::length(a.vel - b.vel) > 1e-6f ||
                rf::length(a.rot.rotate({1,2,3}) - b.rot.rotate({1,2,3})) > 1e-5f)
                throw std::runtime_error("Bullet pose/velocity import roundtrip failed");
        }
        scene->afterStep(sim);
        if (scene->brokenPairs() != 0) throw std::runtime_error("imported scene lost initial linking");
        reference.world.getSolverInfo().m_numIterations = iterations;
        std::ofstream csv(output / "frames.csv"), summary(output / "summary.txt");
        if (!csv || !summary) throw std::runtime_error("cannot write report");
        csv << "frame,time_s,broken_pairs,deepest_bullet_contact_m,engine_ms,observer_ms\n";
        int first = -1, worst = 0;
        float deepest = 0;
        for (int frame = 1; frame <= frames; ++frame) {
            ReferenceTiming timing;
            for (int sub = 0; sub < 10; ++sub) {
                timing.start();
                reference.world.stepSimulation(1.0f/600, 0); // one explicit step, no accumulated-time rounding
                timing.finishEngine();
                reference.copyPoses(sim.rigid);
                scene->afterStep(sim); // topology observer only; rfcore never advances this world
                worst = std::max(worst, scene->brokenPairs());
                if (scene->brokenPairs() && first < 0) first = frame;
                for (int k = 0; k < reference.dispatcher.getNumManifolds(); ++k) {
                    const auto* manifold = reference.dispatcher.getManifoldByIndexInternal(k);
                    for (int p = 0; p < manifold->getNumContacts(); ++p)
                        deepest = std::max(deepest, float(-manifold->getContactPoint(p).getDistance()));
                }
                timing.finishObserver();
            }
            csv << frame << ',' << rf::numberText(double(frame)/60) << ',' << scene->brokenPairs()
                << ',' << rf::numberText(deepest) << timing.csv() << '\n';
            csv.flush();
            if (frame % 30 == 0) std::cout << "frame " << frame << ", worst broken " << worst << std::endl;
        }
        summary << "Bullet " << BT_BULLET_VERSION << "\niterations " << iterations << "\nframes " << frames
                << "\npolyhedral_features " << reference.polyhedral << "\ngeometry_import_error_m " << rf::numberText(importError)
                << "\nsubsteps 10\nhull_margin_m 0\nfirst_break_frame " << first << "\nworst_broken " << worst
                << "\nmax_reported_contact_depth_m " << rf::numberText(deepest)
                << "\nCCD: Bullet defaults; compound bodies have no added swept proxy\n"
                << "sleep disabled; same input hulls, principal inertias, poses and velocities\n"
                << "Bullet solver defaults otherwise; material products matched; no separate static friction\n"
                << "contact depth is Bullet manifold data, not an independent swept-geometry certificate\n";
        summary.close();
        std::cout << "first break " << first << ", worst broken " << worst << ", contact depth "
                  << rf::numberText(deepest) << " m\n";
        return worst ? 1 : 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 2;
    }
}
