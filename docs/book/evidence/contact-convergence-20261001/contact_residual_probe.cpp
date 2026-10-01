// Shadow solve of a real late sandbox state. GCC -fno-access-control is used ONLY for this
// external diagnostic, without changing SDK visibility or the running reference simulation.
#include "tests/TestRunner.h"
#include "core/Format.h"
#include <iostream>

double normalResidual(const RigidWorld& w) {
    double worst = 0;
    for (const auto& m : w.manifolds_) {
        float cfm = 0;
        if (w.params.blockSolver && m.points.size() >= 2 && m.points.size() <= 4) {
            for (int i = 0; i < m.points.size(); ++i) cfm += m.K[i][i];
            cfm *= w.params.blockCfm / m.points.size();
        }
        for (const auto& p : m.points) {
            const float r = dot(w.relativeVelocity(m, p.position), p.normal) - p.velocityBias + cfm * p.jn;
            const double projected = std::fabs(double(p.jn) - std::max(0.0, double(p.jn) - p.massN * r));
            if (p.massN > 0) worst = std::max(worst, projected / p.massN);
        }
    }
    return worst;
}

double tangentResidual(const RigidWorld& w) {
    double worst = 0;
    for (const auto& m : w.manifolds_) {
        const Vector3 v = w.relativeVelocity(m, m.center);
        const float speed = length(v - m.normal * dot(v, m.normal));
        const float mu = speed < w.params.stickVelocity ? m.staticFriction : m.friction;
        float load = 0;
        for (const auto& p : m.points) load += p.jn;
        auto consider = [&](float lambda, float mass, const Vector3& axis) {
            if (mass <= 0) return;
            const float next = clampv(lambda - mass * dot(v, axis), -mu * load, mu * load);
            worst = std::max(worst, double(std::fabs(next - lambda) / mass));
        };
        consider(m.jt1, m.massT1, m.t1); consider(m.jt2, m.massT2, m.t2);
    }
    return worst;
}

void record(const RigidWorld& w, const char* stage, int iteration) {
    std::cout << "{\"stage\":\"" << stage << "\",\"iteration\":" << iteration
              << ",\"kinetic_J\":" << numberText(w.kineticEnergy())
              << ",\"normal_projected_residual_m_s\":" << numberText(normalResidual(w))
              << ",\"tangent_projected_residual_m_s\":" << numberText(tangentResidual(w))
              << ",\"cup_spin_rad_s\":" << numberText(length(w.bodies()[23].angVel)) << "}\n";
}

int main() {
    Simulation sim; loadSample(sim, Preset::RigidSandbox); sim.rigid.params.sleeping = false;
    for (int i = 0; i < 480; ++i) sim.stepFrame();
    auto& w = sim.rigid;
    const float h = sim.frameDt / float(w.params.substeps);
    w.lastDt_ = w.correctionDt_ = h;
    w.beginStep(); w.integrateVelocities(h); w.collide(); w.prepare(h);
    record(w, "warm-start", 0);
    for (int i = 1; i <= 96; ++i) {
        w.solve();
        if (i <= 12 || i == 24 || i == 48 || i == 96) record(w, "contact-iterations", i);
    }
    w.applyRestitution(); record(w, "restitution", 96);
    w.propagateShock(); record(w, "shock", 96);
    const auto& cup = w.bodies()[23];
    std::cout << "{\"cup_position\":[" << numberText(cup.pos.x) << ',' << numberText(cup.pos.y) << ',' << numberText(cup.pos.z)
              << "],\"cup_omega\":[" << numberText(cup.angVel.x) << ',' << numberText(cup.angVel.y) << ',' << numberText(cup.angVel.z) << "]}\n";
    for (const auto& m : w.manifolds_) if (m.a == 23 || m.b == 23) {
        std::cout << "{\"cup_contact\":true,\"a\":" << m.a << ",\"b\":" << m.b << ",\"locked\":" << m.locked
                  << ",\"points\":" << m.points.size() << ",\"normal\":[" << numberText(m.normal.x) << ',' << numberText(m.normal.y) << ',' << numberText(m.normal.z)
                  << "],\"center\":[" << numberText(m.center.x) << ',' << numberText(m.center.y) << ',' << numberText(m.center.z)
                  << "],\"twist\":" << numberText(m.jtwist) << ",\"roll\":[" << numberText(m.jroll.x) << ',' << numberText(m.jroll.y) << ',' << numberText(m.jroll.z) << "]}\n";
        for (const auto& p : m.points)
            std::cout << "{\"depth_m\":" << numberText(p.depth) << ",\"load_Ns\":" << numberText(p.jn)
                      << ",\"vn_m_s\":" << numberText(dot(w.relativeVelocity(m, p.position),p.normal))
                      << ",\"target_m_s\":" << numberText(p.velocityBias) << "}\n";
    }
}
