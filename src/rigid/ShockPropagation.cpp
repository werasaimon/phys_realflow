// Shock propagation of RigidWorld (Guendelman, Bridson, Fedkiw 2003): every body's level in
// the contact graph - its distance from the static environment - and the one-sided solve from
// the ground up in which the lower body is infinitely heavy.
#include "rigid/RigidWorld.h"

#include "core/Parallel.h"

#include <algorithm>

namespace rf {

// Shock propagation (Guendelman, Bridson & Fedkiw 2003): the manifolds are solved once more from
// the ground up, each level against a frozen support, so the weight of a tall stack reaches the
// floor in one sweep instead of one level per iteration.
//
// In parallel within a level, where it is safe: a one-sided push writes only its upper body (the
// support below is frozen, its level done), so the pushes of one level whose upper bodies differ
// run at the same time - in a pile of a thousand cubes a level holds about a hundred of them. The
// two-sided solves of bodies on the same level do reach each other (a row of touching balls passes
// a hit along it); they keep their ground-up order, one by one, after the level's pushes - solved
// in coloured batches instead, a Newton's cradle of touching balls gained 4 % momentum. The batches
// (buildShockBatches) depend only on the contact graph, never on the number of threads, so the
// answer is the same bits on 1 or 24 threads.
void RigidWorld::propagateShock() {
    computeLevels();
    buildShockBatches();
    for (int pass = 0; pass < params.shockIterations; ++pass) {
        shockFrictionPass_ = params.shockFriction && pass == params.shockIterations - 1; // once per step
        for (int b = 0; b < shockBatchCount_; ++b) {
            const std::vector<int>& batch = shockBatches_[size_t(b)];
            if (!shockBatchSerial_[size_t(b)] && batch.size() >= kShockParallelBatch)
                parallelFor(int(batch.size()), [&](int k) { solveManifoldShock(manifolds_[size_t(batch[size_t(k)])]); }, 8);
            else
                for (int i : batch) solveManifoldShock(manifolds_[size_t(i)]);
        }
    }
}

// The batches of the shock pass, ground up:
//   1. every manifold's level - the lower of its two bodies' (the static world counts as -1);
//   2. the manifolds in that order, ties by the bottom-up solve order (solveRank_);
//   3. within each level, the one-sided pushes coloured greedily by their upper body (as the
//      regular iterations' buildColors): a push takes the first colour its upper body does not have
//      yet in this level; colour 63 collects the rest, solved one by one; the two-sided solves of
//      the level get colour 64 - the last batch, one by one, in ground-up order;
//   4. the batches: level by level, colour by colour, each in ground-up order.
void RigidWorld::buildShockBatches() {
    const int n = int(manifolds_.size());
    auto levelOf = [&](int body) { return body < 0 || bodies_[body].invMass == 0 ? -1 : levels_[size_t(body)]; };
    shockKeys_.resize(size_t(n));
    shockOrder_.resize(size_t(n));
    for (int i = 0; i < n; ++i) { // 1.
        const Manifold& m = manifolds_[size_t(i)];
        shockKeys_[size_t(i)] = std::min(levelOf(m.a), levelOf(m.b));
        shockOrder_[size_t(i)] = i;
    }
    std::sort(shockOrder_.begin(), shockOrder_.end(), [&](int x, int y) { // 2.
        const int kx = shockKeys_[size_t(x)], ky = shockKeys_[size_t(y)];
        return kx < ky || (kx == ky && solveRank_[size_t(x)] < solveRank_[size_t(y)]);
    });
    shockColorOf_.resize(size_t(n));
    bodyShockLevel_.assign(bodies_.size(), kNoShockLevel);
    bodyShockColors_.resize(bodies_.size());
    shockBatchCount_ = 0;
    for (int first = 0; first < n;) {
        const int level = shockKeys_[size_t(shockOrder_[size_t(first)])];
        int last = first;
        while (last < n && shockKeys_[size_t(shockOrder_[size_t(last)])] == level) ++last;
        bool twoSided = false;
        const int colors = colorShockLevel(first, last, level, twoSided); // 3.
        for (int c = 0; c < colors; ++c) addShockBatch(first, last, c); // 4.
        if (twoSided) addShockBatch(first, last, kShockTwoSided);
        first = last;
    }
}

// Step 3 of buildShockBatches for the manifolds shockOrder_[first, last) of one level: each one's
// colour into shockColorOf_ (kShockTwoSided: a two-sided solve). Returns how many colours the
// one-sided pushes used; `twoSided` says whether the level has two-sided solves too.
int RigidWorld::colorShockLevel(int first, int last, int level, bool& twoSided) {
    int used = 0;
    auto colorsOf = [&](int body) -> uint64_t& { // the colours `body` already has in this level
        if (bodyShockLevel_[size_t(body)] != level) bodyShockLevel_[size_t(body)] = level, bodyShockColors_[size_t(body)] = 0;
        return bodyShockColors_[size_t(body)];
    };
    auto levelOf = [&](int body) { return body < 0 || bodies_[body].invMass == 0 ? -1 : levels_[size_t(body)]; };
    for (int k = first; k < last; ++k) {
        const int i = shockOrder_[size_t(k)];
        const Manifold& m = manifolds_[size_t(i)];
        const int la = levelOf(m.a), lb = levelOf(m.b);
        if (la == lb) { // a two-sided solve: the level's last batch, in order
            shockColorOf_[size_t(i)] = kShockTwoSided;
            twoSided = true;
            continue;
        }
        const int upper = la > lb ? m.a : m.b;
        const uint64_t taken = colorsOf(upper);
        int c = 63;
        for (int bit = 0; bit < 63; ++bit)
            if (!(taken & (uint64_t(1) << bit))) { c = bit; break; }
        shockColorOf_[size_t(i)] = c;
        colorsOf(upper) |= uint64_t(1) << c;
        used = std::max(used, c + 1);
    }
    return used;
}

// Step 4 of buildShockBatches: the batch of colour c of the level shockOrder_[first, last), in
// ground-up order; the overflow colour and the two-sided solves share bodies, so they run one by one.
void RigidWorld::addShockBatch(int first, int last, int c) {
    std::vector<int>& batch = nextShockBatch();
    for (int k = first; k < last; ++k)
        if (shockColorOf_[size_t(shockOrder_[size_t(k)])] == c) batch.push_back(shockOrder_[size_t(k)]);
    if (batch.empty()) --shockBatchCount_; // a colour left unused: give the batch back
    else shockBatchSerial_[size_t(shockBatchCount_ - 1)] = c >= 63;
}

// A cleared batch at the end of the list (kept between steps: no allocation once grown).
std::vector<int>& RigidWorld::nextShockBatch() {
    if (shockBatchCount_ == int(shockBatches_.size())) shockBatches_.emplace_back(), shockBatchSerial_.push_back(0);
    std::vector<int>& batch = shockBatches_[size_t(shockBatchCount_++)];
    batch.clear();
    return batch;
}

void RigidWorld::computeLevels() {
    // BFS over the contact graph starting from bodies touching the static environment (level 0).
    const int n = int(bodies_.size());
    const int kFar = 1 << 28;
    levels_.assign(n, kFar);
    // Adjacency in compressed rows, in vectors kept between steps: no allocation at all.
    std::vector<int>&start = levelStart_, &adj = levelAdj_, &queue = levelQueue_, &fill = levelFill_;
    start.assign(n + 1, 0);
    adj.clear();
    queue.clear();
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
    fill.assign(start.begin(), start.end() - 1);
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

// A contact is solved one-sided only while its support does not move into the upper body faster
// than this, along the normal at the contact [m/s].
constexpr float kRestingSupport = 0.1f;

void RigidWorld::solveManifoldShock(Manifold& m) {
    // Level of each side: static environment counts as -1 (always "below").
    const int la = bodies_[m.a].invMass == 0 ? -1 : levels_[m.a];
    const int lb = (m.b < 0 || bodies_[m.b].invMass == 0) ? -1 : levels_[m.b];
    const bool upperIsA = la > lb;
    // Shock propagation is for resting support. A support frozen while it moves into the upper
    // body throws it off like a moving wall, with momentum the support never gives: a chain let go
    // stretched out sideways gained 171 J in the first second of its fall. A support that sinks
    // away (a stack under a cascade of landings) only takes the upper body's approach away. So a
    // support moving into the upper body gets the ordinary two-sided solve.
    const int support = upperIsA ? m.b : m.a;
    bool resting = true;
    if (support >= 0 && bodies_[size_t(support)].invMass > 0) {
        const RigidBody& s = bodies_[size_t(support)];
        const Vector3 toUpper = upperIsA ? m.normal : -m.normal; // m.normal points from b to a
        resting = dot(s.vel + cross(s.angVel, m.center - s.pos), toUpper) < kRestingSupport;
    }
    if (la == lb || !resting) { // two-sided - except for an impact, whose separation (applied by
        // applyRestitution just before) a two-sided re-solve would take back; shock propagation
        // is for resting support, not for a bullet meeting a box in the air.
        bool impact = false;
        for (const SolverPoint& p : m.points) impact |= p.bounce > 0;
        if (!impact) solveManifold(m);
        return;
    }
    RigidBody& upper = bodies_[size_t(upperIsA ? m.a : m.b)];
    float pushes[4] = {0, 0, 0, 0};
    const int np = pushUpperOffSupport(m, upperIsA, pushes);
    if (shockFrictionPass_) dragAlongSupport(m, upper, upperIsA, pushes, np);
}

// Velocity correction of the upper body only, with its own impulses (starting at zero) so the
// warm-start impulses stay those of the two-sided solve; friction is left to the regular
// iterations. The upper body against a frozen support is a tiny problem - at most four points
// pushing one body - solved by eight projected Gauss-Seidel sweeps, on a small matrix made once:
//   1. every point's push direction on the upper body d (the normal, turned to point away from the
//      support), its arm r, and b: the normal velocity the point has now minus the one it may keep
//      (a speculative gap may still close, never pull); the support is frozen;
//   2. nothing to do if no point needs a push (b >= 0 everywhere: the regular iterations already
//      hold it, or it separates) - the common case, known before the matrix is made;
//   3. the matrix K: how a unit push at point j changes the normal velocity at point i, with the
//      upper body's mass alone (the support does not give);
//   4. the sweeps: each point's push made just large enough (never negative) for its velocity to
//      reach its target, given the others' pushes; a sweep that changed nothing ends them;
//   5. the pushes applied to the upper body, once.
// The same sweeps once ran on the bodies themselves - a velocity and an impulse at every point of
// every sweep - and were the most expensive stage of the whole step; on the 4 x 4 matrix a sweep
// is 16 multiply-adds. (An exact solve of the same small LCP - ContactLcp.h - was tried: it let
// bullets push boxes 1 cm into each other in the CCD test, where the sweeps keep 1 mm.)
// Returns the number of points; `pushes` gets each point's push.
int RigidWorld::pushUpperOffSupport(Manifold& m, bool upperIsA, float* pushes) {
    constexpr int kSweeps = 8;
    const int ui = upperIsA ? m.a : m.b;
    const RigidBody& U = bodies_[ui];
    const float s = upperIsA ? 1.0f : -1.0f; // impulses are expressed for A (normal from B to A)
    const int np = std::min<int>(int(m.points.size()), 4);
    float b[4], K[4][4];
    Vector3 d[4], r[4];
    bool anyPush = false;
    for (int i = 0; i < np; ++i) { // 1.
        const SolverPoint& p = m.points[i];
        d[i] = p.normal * s;
        r[i] = p.position - U.pos;
        Vector3 support(0.0f);
        if (upperIsA) { if (m.b >= 0) support = bodies_[m.b].velocityAt(p.position); }
        else support = bodies_[m.a].velocityAt(p.position);
        b[i] = dot(U.velocityAt(p.position) - support, d[i]) - std::min(p.velocityBias, 0.0f);
        anyPush |= b[i] < 0;
    }
    if (!anyPush) return np; // 2.
    for (int j = 0; j < np; ++j) { // 3.
        const Vector3 turn = U.applyInvInertiaWorld(cross(r[j], d[j])); // the spin a unit push at j gives
        for (int i = 0; i < np; ++i) K[i][j] = dot(d[i], d[j] * U.invMass + cross(turn, r[i]));
    }
    for (int sweep = 0; sweep < kSweeps; ++sweep) { // 4.
        bool changed = false;
        for (int i = 0; i < np; ++i) {
            if (K[i][i] < 1e-12f) continue;
            float velocity = b[i];
            for (int j = 0; j < np; ++j) velocity += K[i][j] * pushes[j];
            const float pushed = std::max(0.0f, pushes[i] - velocity / K[i][i]);
            changed |= pushed != pushes[i];
            pushes[i] = pushed;
        }
        if (!changed) break;
    }
    for (int i = 0; i < np; ++i) // 5.
        if (pushes[i] > 0) applyImpulse(ui, d[i] * pushes[i], m.points[i].position);
    return np;
}

// Friction in the same one-directional pass: the upper body is dragged along by (or held on) its
// frozen support up to the Coulomb limit, so a sideways pull reaches the top of a column in one
// sweep instead of one level per iteration (which shears a stack into a "staircase"). `acc` are
// the normal impulses this pass just applied on top of the manifold's own.
void RigidWorld::dragAlongSupport(Manifold& m, RigidBody& U, bool upperIsA, const float* acc, int np) {
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
