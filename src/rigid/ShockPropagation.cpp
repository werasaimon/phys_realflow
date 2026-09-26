// Shock propagation of RigidWorld (Guendelman, Bridson, Fedkiw 2003): every body's level in
// the contact graph - its distance from the static environment - and the one-sided solve from
// the ground up in which the lower body is infinitely heavy.
#include "rigid/RigidWorld.h"

#include <algorithm>

namespace rf {

void RigidWorld::computeLevels() {
    // BFS over the contact graph starting from bodies touching the static environment (level 0).
    const int n = int(bodies_.size());
    const int kFar = 1 << 28;
    levels_.assign(n, kFar);
    // Adjacency in compressed rows (one allocation instead of one vector per body).
    std::vector<int> start(n + 1, 0), adj, queue;
    auto dynamicPair = [&](const Manifold& m) {
        return m.b >= 0 && bodies_[m.a].invMass > 0 && bodies_[m.b].invMass > 0;
    };
    for (const Manifold& m : manifolds_) {
        bool bStatic = m.b < 0 || bodies_[m.b].invMass == 0;
        bool aStatic = bodies_[m.a].invMass == 0;
        if (bStatic && !aStatic && levels_[m.a] != 0) { levels_[m.a] = 0; queue.push_back(m.a); }
        if (aStatic && m.b >= 0 && !bStatic && levels_[m.b] != 0) { levels_[m.b] = 0; queue.push_back(m.b); }
        if (dynamicPair(m)) ++start[m.a + 1], ++start[m.b + 1];
    }
    for (int i = 0; i < n; ++i) start[i + 1] += start[i];
    adj.resize(start[n]);
    std::vector<int> fill(start.begin(), start.end() - 1);
    for (const Manifold& m : manifolds_)
        if (dynamicPair(m)) adj[fill[m.a]++] = m.b, adj[fill[m.b]++] = m.a;
    for (size_t h = 0; h < queue.size(); ++h) {
        int u = queue[h];
        for (int e = start[u]; e < start[u + 1]; ++e) {
            int v = adj[e];
            if (levels_[v] > levels_[u] + 1) { levels_[v] = levels_[u] + 1; queue.push_back(v); }
        }
    }
}

void RigidWorld::solveManifoldShock(Manifold& m) {
    // Level of each side: static environment counts as -1 (always "below").
    const int la = bodies_[m.a].invMass == 0 ? -1 : levels_[m.a];
    const int lb = (m.b < 0 || bodies_[m.b].invMass == 0) ? -1 : levels_[m.b];
    if (la == lb) { solveManifold(m); return; } // same level: ordinary two-sided solve
    const bool upperIsA = la > lb;
    const int ui = upperIsA ? m.a : m.b;
    RigidBody& U = bodies_[ui];
    const float s = upperIsA ? 1.0f : -1.0f; // impulses are expressed for A (normal from B to A)
    // Velocity correction of the upper body only. Uses its own accumulators (starting at zero) so the
    // warm-start impulses stay those of the two-sided solve; friction is left to the regular
    // iterations. The upper body against a frozen support is a tiny problem (<= 4 points, one
    // body): iterate it locally until converged.
    float acc[16] = {0};
    const int np = std::min<int>(int(m.points.size()), 16);
    for (int local = 0; local < 8; ++local)
        for (int k = 0; k < np; ++k) {
            SolverPoint& p = m.points[k];
            Vector3 r = p.position - U.pos;
            Vector3 vu = U.velocityAt(p.position);
            Vector3 vl(0.0f);
            if (upperIsA) { if (m.b >= 0) vl = bodies_[m.b].velocityAt(p.position); }
            else vl = bodies_[m.a].velocityAt(p.position);
            float vn = dot(upperIsA ? vu - vl : vl - vu, p.normal);
            float kk = U.invMass + dot(cross(U.applyInvInertiaWorld(cross(r, p.normal)), r), p.normal);
            if (kk < 1e-12f) continue;
            float target = std::min(p.velocityBias, 0.0f); // speculative gap may still close; never push apart
            float old = acc[k];
            acc[k] = std::max(old + (target - vn) / kk, 0.0f);
            applyImpulse(ui, p.normal * ((acc[k] - old) * s), p.position);
        }

    // Friction in the same one-directional pass: the upper body is dragged along by (or held on)
    // its frozen support up to the Coulomb limit, so a sideways pull reaches the top of a column in
    // one sweep instead of one level per iteration (which shears a stack into a "staircase").
    if (!shockFrictionPass_) return;
    float normalTotal = 0;
    Vector3 c(0.0f), nrm(0.0f);
    for (int k = 0; k < np; ++k) {
        normalTotal += m.points[k].jn + acc[k];
        c += m.points[k].position;
        nrm += m.points[k].normal;
    }
    if (normalTotal <= 0) return;
    c /= float(np);
    nrm = normalize(nrm);
    Vector3 vu = U.velocityAt(c);
    Vector3 vl(0.0f);
    if (upperIsA) { if (m.b >= 0) vl = bodies_[m.b].velocityAt(c); }
    else vl = bodies_[m.a].velocityAt(c);
    Vector3 vrel = vu - vl;                          // velocity of the upper body relative to its support
    Vector3 vt = vrel - nrm * dot(vrel, nrm);
    float vtl = length(vt);
    if (vtl < 1e-7f) return;
    // Pure translation (no torque): a correction pass must not spin bodies up; tipping torques are
    // left to the regular two-sided friction.
    Vector3 t = vt / vtl;
    // Two limits:
    //  * the support's friction has to move everything resting on it: in equilibrium the normal
    //    impulse equals the weight of the whole column above times h, so the Coulomb-limited
    //    velocity change of that column is mu * |g| * h per step, independent of its height (this
    //    also keeps the pure translation - it ignores the rotation - a small correction);
    //  * the Coulomb cone of this contact: mu times its normal impulse, minus the friction impulse
    //    the regular iterations have already applied. The pass only finishes what they could not
    //    carry up the column; it never adds friction beyond Coulomb (a sliding body keeps mu_k N).
    const float mu = vtl < params.stickVelocity ? m.staticFriction : m.friction;
    const float used = length(m.t1 * m.jt1 + m.t2 * m.jt2);
    const float budget = std::max(0.0f, mu * normalTotal - used);
    const float dv = std::min({vtl, mu * length(params.gravity) * lastDt_, budget * U.invMass});
    U.vel -= t * dv;
}

} // namespace rf
