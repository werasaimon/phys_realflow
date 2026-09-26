#include "particles/Cloth.h"

#include "core/Parallel.h"

#include <algorithm>
#include <queue>

namespace rf {

static uint64_t pairKey(int a, int b) {
    return (uint64_t(uint32_t(std::min(a, b))) << 32) | uint32_t(std::max(a, b));
}

static bool contains(const std::vector<int>& v, int x) { return std::find(v.begin(), v.end(), x) != v.end(); }

void buildClothConstraints(Cloth& cloth, const std::vector<Vector3>& rest) {
    const int w = cloth.width, h = cloth.height, first = cloth.firstParticle;
    const ClothMaterial& m = cloth.material;
    const float s = cloth.spacing;
    // A thread between neighbours stands for a strip of width s and length s: spring k = stiffness.
    const float threadCompliance = m.tensileStiffness > 0 ? 1.0f / m.tensileStiffness : 0.0f;
    const float shearCompliance = m.shearStiffness > 0 ? 1.0f / m.shearStiffness : 0.0f;
    // Colour of every constraint type / position so that equal colours never share a particle:
    // threads and shear alternate every other column (row), bending every fourth.
    std::vector<std::vector<DistanceConstraint>> byColor(16);
    auto inside = [&](int x, int y) { return x >= 0 && y >= 0 && x < w && y < h; };
    auto add = [&](int color, DistanceConstraint::Kind kind, int x0, int y0, int x1, int y1, float compliance, float strength) {
        if (!inside(x0, y0) || !inside(x1, y1)) return;
        DistanceConstraint c;
        c.a = first + x0 + w * y0;
        c.b = first + x1 + w * y1;
        c.restLength = length(rest[x0 + w * y0] - rest[x1 + w * y1]);
        c.compliance = compliance;
        c.strength = strength;
        c.kind = kind;
        c.x = int16_t(std::min(x0, x1));
        c.y = int16_t(std::min(y0, y1));
        byColor[color].push_back(c);
    };
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x) {
            float warp = m.strengthWarp * s * (contains(m.seamColumns, x) ? m.seamStrength : 1.0f);
            float weft = m.strengthWeft * s * (contains(m.seamRows, y) ? m.seamStrength : 1.0f);
            add(0 + x % 2, DistanceConstraint::Warp, x, y, x + 1, y, threadCompliance, warp);
            add(2 + y % 2, DistanceConstraint::Weft, x, y, x, y + 1, threadCompliance, weft);
            add(4 + x % 2, DistanceConstraint::Shear, x, y, x + 1, y + 1, shearCompliance, 0);
            add(6 + x % 2, DistanceConstraint::Shear, x + 1, y, x, y + 1, shearCompliance, 0);
            add(8 + x % 4, DistanceConstraint::Bend, x, y, x + 2, y, m.bendCompliance, 0);
            add(12 + y % 4, DistanceConstraint::Bend, x, y, x, y + 2, m.bendCompliance, 0);
        }
    cloth.constraints.clear();
    cloth.batchStart.assign(1, 0);
    for (const auto& batch : byColor) {
        cloth.constraints.insert(cloth.constraints.end(), batch.begin(), batch.end());
        cloth.batchStart.push_back(int(cloth.constraints.size()));
    }
    cloth.rest = rest;
    cloth.constraintIndex.clear();
    for (int k = 0; k < int(cloth.constraints.size()); ++k)
        cloth.constraintIndex[pairKey(cloth.constraints[k].a, cloth.constraints[k].b)] = k;
    cloth.cellIntact.assign(size_t(w - 1) * (h - 1), 1);
    cloth.warpBroken.assign(size_t(w - 1) * h, 0);
    cloth.weftBroken.assign(size_t(w) * (h - 1), 0);
    cloth.tornThreads = 0;
    cloth.temperature.assign(size_t(w) * h, 0.0f);
    cloth.unburnt.assign(size_t(w) * h, 1.0f);
    cloth.burntThreads = 0;
}

// Tether lengths: the distance along the cloth from every anchor - the straight rest distance on an
// intact sheet, shortest paths along the intact threads and shear links (Dijkstra) once torn; a
// particle cut off from an anchor is not held by it. A cloth that can tear gets 15 % slack - more
// than its threads stretch before they break - so the tethers are only a safety net and the
// threads carry (and concentrate) the load; one that cannot tear keeps them tight.
static void measureTethers(Cloth& cloth) {
    const int count = cloth.width * cloth.height, first = cloth.firstParticle, per = cloth.tethersPerParticle;
    if (per == 0) return;
    const bool tearable = cloth.material.strengthWarp > 0 || cloth.material.strengthWeft > 0;
    const float slack = tearable ? 1.15f : 1.0f;
    if (cloth.tornThreads == 0 && cloth.burntThreads == 0) { // intact (threads also break by burning)
        for (int i = 0; i < count; ++i)
            for (int k = 0; k < per; ++k) {
                Tether& t = cloth.tethers[size_t(i) * per + k];
                t.maxLength = length(cloth.rest[i] - cloth.rest[t.anchor - first]) * slack;
            }
        return;
    }
    // Torn: one multi-source shortest-path search (Dijkstra) along the intact threads and shear
    // links finds, for every particle, the nearest anchor it is still connected to and how far along
    // the cloth it is; the particle is tied to that anchor only (unreachable: not tied at all).
    std::vector<std::vector<std::pair<int, float>>> adj(count);
    for (const DistanceConstraint& c : cloth.constraints) {
        if (c.broken || c.kind == DistanceConstraint::Bend) continue;
        adj[c.a - first].push_back({c.b - first, c.restLength});
        adj[c.b - first].push_back({c.a - first, c.restLength});
    }
    std::vector<float> dist(count, kInf);
    std::vector<int> source(count, -1);
    using Item = std::pair<float, int>;
    std::priority_queue<Item, std::vector<Item>, std::greater<Item>> open;
    for (int a : cloth.anchors) {
        dist[a] = 0;
        source[a] = a;
        open.push({0.0f, a});
    }
    while (!open.empty()) {
        auto [d, u] = open.top();
        open.pop();
        if (d > dist[u]) continue;
        for (auto [v, len] : adj[u])
            if (d + len < dist[v]) {
                dist[v] = d + len;
                source[v] = source[u];
                open.push({dist[v], v});
            }
    }
    for (int i = 0; i < count; ++i)
        for (int k = 0; k < per; ++k) {
            Tether& t = cloth.tethers[size_t(i) * per + k];
            if (k == 0 && source[i] >= 0) {
                t.anchor = first + source[i];
                t.maxLength = dist[i] * slack;
            } else {
                t.maxLength = kInf;
            }
        }
}

void updateTethers(Cloth& cloth) {
    if (!cloth.tethersDirty) return;
    measureTethers(cloth);
    cloth.tethersDirty = false;
}

void buildTethers(Cloth& cloth, const std::vector<float>& invMass) {
    const int count = cloth.width * cloth.height, first = cloth.firstParticle;
    std::vector<int>& anchors = cloth.anchors;
    anchors.clear();
    for (int a = 0; a < count; ++a)
        if (invMass[first + a] == 0) anchors.push_back(a);
    // Every particle is tied to its (at most) four nearest anchors: for a cloth pinned along a whole
    // edge the nearest ones hold it, and the cost stays proportional to the particle count.
    const int per = std::min<int>(4, int(anchors.size()));
    cloth.tethers.clear();
    cloth.tethersPerParticle = per;
    for (int i = 0; i < count; ++i) {
        std::vector<int> nearest = anchors;
        std::partial_sort(nearest.begin(), nearest.begin() + per, nearest.end(), [&](int a, int b) {
            return length2(cloth.rest[a] - cloth.rest[i]) < length2(cloth.rest[b] - cloth.rest[i]);
        });
        for (int k = 0; k < per; ++k) cloth.tethers.push_back({first + i, first + nearest[k], kInf});
    }
    measureTethers(cloth);
}

static void solveConstraint(DistanceConstraint& c, std::vector<Vector3>& p, const std::vector<float>& invMass, float invDt2) {
    const float wa = invMass[c.a], wb = invMass[c.b];
    if (c.broken || wa + wb == 0) return;
    Vector3 d = p[c.a] - p[c.b];
    float len = length(d);
    if (len < 1e-9f) return;
    Vector3 n = d / len;
    float C = len - c.restLength;
    float alpha = c.compliance * invDt2;
    float dl = (-C - alpha * c.lambda) / (wa + wb + alpha);
    c.lambda += dl;
    p[c.a] += n * (wa * dl);
    p[c.b] -= n * (wb * dl);
}

void solveCloth(Cloth& cloth, std::vector<Vector3>& p, const std::vector<float>& invMass, float dt) {
    const float invDt2 = 1.0f / (dt * dt);
    if (cloth.constraints.size() < 60000) { // up to ~10k particles one thread beats the batches' sync overhead
        for (DistanceConstraint& c : cloth.constraints) solveConstraint(c, p, invMass, invDt2);
    } else {
        for (size_t b = 0; b + 1 < cloth.batchStart.size(); ++b) {
            const int begin = cloth.batchStart[b], count = cloth.batchStart[b + 1] - begin;
            parallelFor(count, [&](int k) { solveConstraint(cloth.constraints[begin + k], p, invMass, invDt2); }, 256);
        }
    }
    // Tethers: one-sided, back onto the sphere of the allowed distance around each anchor.
    const int per = cloth.tethersPerParticle;
    if (per == 0) return;
    parallelFor(int(cloth.tethers.size() / per), [&](int i) {
        for (int k = 0; k < per; ++k) {
            const Tether& t = cloth.tethers[size_t(i) * per + k];
            if (invMass[t.particle] == 0) continue;
            Vector3 d = p[t.particle] - p[t.anchor];
            float len = length(d);
            if (len > t.maxLength) p[t.particle] = p[t.anchor] + d * (t.maxLength / len);
        }
    }, 256);
}

float threadTension(const DistanceConstraint& c, const std::vector<Vector3>& p, float dt) {
    // Elastic thread: T = k * elongation; inextensible: the constraint force -lambda / dt^2.
    if (c.compliance > 0) return (length(p[c.a] - p[c.b]) - c.restLength) / c.compliance;
    return -c.lambda / (dt * dt);
}

// A thread (warp or weft) breaks: its bending links across the break go with it, and the cells on
// both sides fall apart once both of their threads in this direction are gone.
static void breakThread(Cloth& cloth, DistanceConstraint& c) {
    const int w = cloth.width, h = cloth.height, first = cloth.firstParticle;
    auto breakPair = [&](int xa, int ya, int xb, int yb) {
        if (xa < 0 || ya < 0 || xb < 0 || yb < 0 || xa >= w || xb >= w || ya >= h || yb >= h) return;
        auto it = cloth.constraintIndex.find(pairKey(first + xa + w * ya, first + xb + w * yb));
        if (it != cloth.constraintIndex.end()) cloth.constraints[it->second].broken = true;
    };
    auto warpBroken = [&](int x, int y) { return x >= 0 && y >= 0 && x < w - 1 && y < h && cloth.warpBroken[x + (w - 1) * y]; };
    auto weftBroken = [&](int x, int y) { return x >= 0 && y >= 0 && x < w && y < h - 1 && cloth.weftBroken[x + w * y]; };
    auto cutCell = [&](int x, int y) { // the crack went through: the cell falls apart
        if (x < 0 || y < 0 || x >= w - 1 || y >= h - 1 || !cloth.cellIntact[x + (w - 1) * y]) return;
        cloth.cellIntact[x + (w - 1) * y] = 0;
        breakPair(x, y, x + 1, y + 1);
        breakPair(x + 1, y, x, y + 1);
    };

    c.broken = true;
    cloth.tethersDirty = true;
    const int x = c.x, y = c.y;
    if (c.kind == DistanceConstraint::Warp) {
        cloth.warpBroken[x + (w - 1) * y] = 1;
        breakPair(x - 1, y, x + 1, y); // bending links across the break
        breakPair(x, y, x + 2, y);
        // The cells above and below: cut once both of their warp threads are gone.
        if (warpBroken(x, y - 1)) cutCell(x, y - 1);
        if (warpBroken(x, y + 1)) cutCell(x, y);
    } else {
        cloth.weftBroken[x + w * y] = 1;
        breakPair(x, y - 1, x, y + 1);
        breakPair(x, y, x, y + 2);
        if (weftBroken(x - 1, y)) cutCell(x - 1, y);
        if (weftBroken(x + 1, y)) cutCell(x, y);
    }
}

// Part of a thread's strength left after burning (1 for fabric that does not burn).
static float strengthLeft(const Cloth& cloth, const DistanceConstraint& c) {
    if (!cloth.material.flammable) return 1.0f;
    const float u = std::min(cloth.unburnt[c.a - cloth.firstParticle], cloth.unburnt[c.b - cloth.firstParticle]);
    return cloth.material.charStrength + (1.0f - cloth.material.charStrength) * u;
}

int tearCloth(Cloth& cloth, const std::vector<Vector3>& p, float dt) {
    int broken = 0;
    for (DistanceConstraint& c : cloth.constraints) {
        if (c.broken || c.strength <= 0 || threadTension(c, p, dt) <= c.strength * strengthLeft(cloth, c)) continue;
        breakThread(cloth, c);
        ++broken;
    }
    cloth.tornThreads += broken;
    return broken;
}

int burnCloth(Cloth& cloth, const std::vector<float>& gasTemperature, const std::vector<float>& irradiance,
              float ambientTemperature, float dt, std::vector<float>& heatToGas, std::vector<float>& fuelToGas) {
    const ClothMaterial& m = cloth.material;
    if (!m.flammable) return 0;
    const float sigma = 5.670e-8f; // Stefan-Boltzmann [W/(m^2 K^4)]
    // The patch of fabric one particle stands for (the same share of the sheet as its mass).
    const float area = cloth.particleArea > 0 ? cloth.particleArea : cloth.spacing * cloth.spacing;
    const float freshMass = m.areaDensity * area;
    const float hA = 2.0f * m.heatTransfer * area;    // convection [W/K], both faces
    const float T0 = ambientTemperature;
    auto heatCapacityOf = [&](size_t k) { // [J/K] of a patch, lighter as it burns
        return freshMass * (m.charMassFraction + (1.0f - m.charMassFraction) * cloth.unburnt[k]) * m.specificHeat;
    };

    // Conduction (Fourier): between two patches joined by an intact thread flows k t (Tb - Ta) -
    // a strip of cross-section s x t and length s. Explicit: k t dt / C is far below 1.
    {
        const float conductance = m.conductivity * m.thickness; // [W/K]
        std::vector<float> dT(cloth.temperature.size(), 0.0f);
        const int first = cloth.firstParticle;
        for (const DistanceConstraint& c : cloth.constraints) {
            if (c.broken || (c.kind != DistanceConstraint::Warp && c.kind != DistanceConstraint::Weft)) continue;
            const int a = c.a - first, b = c.b - first;
            const float Q = conductance * (cloth.temperature[b] - cloth.temperature[a]) * dt; // [J] a <- b
            dT[a] += Q / heatCapacityOf(a);
            dT[b] -= Q / heatCapacityOf(b);
        }
        for (size_t k = 0; k < dT.size(); ++k) cloth.temperature[k] += dT[k];
    }

    const float R = 8.314f; // gas constant [J/(mol K)]
    auto pyrolysisRate = [&](float Tabove) { // Arrhenius [1/s]
        return m.pyrolysisPreExponential * std::exp(-m.pyrolysisActivationEnergy / (R * (Tabove + T0)));
    };
    for (size_t k = 0; k < cloth.unburnt.size(); ++k) {
        float& T = cloth.temperature[k];
        float& unburnt = cloth.unburnt[k];
        const float Tgas = gasTemperature[k];
        const float heatCapacity = heatCapacityOf(k);

        // 1) Heat in: radiation absorbed from the flame (the patch taken as facing it) minus its own
        //    emission from both faces; convection relaxes the patch towards the gas temperature
        //    shifted by that balance - exact for the step with the radiation held constant.
        const float Tabs = T + T0;
        const float radiation = m.emissivity * area * (irradiance[k] - 2.0f * sigma * (Tabs * Tabs * Tabs * Tabs - T0 * T0 * T0 * T0));
        float Theated = hA > 0 ? Tgas + radiation / hA + (T - Tgas - radiation / hA) * std::exp(-hA / heatCapacity * dt)
                               : T + radiation * dt / heatCapacity; // no convection: radiation alone
        heatToGas[k] += hA * (0.5f * (T + Theated) - Tgas) * dt;    // a hotter patch warms the gas

        // 2) Pyrolysis: the part u (1 - exp(-k(T) dt)) of what is left decomposes and takes its heat
        //    from the patch. k(T) is stiff (x10 per ~50 K), so the end temperature is solved from
        //    the energy balance C (T - Theated) + H m u (1 - exp(-k(T) dt)) = 0 by bisection: the
        //    left side grows with T, it is >= 0 at Theated and <= 0 if everything decomposed.
        const float volatileMass = freshMass * (1.0f - m.charMassFraction) * unburnt; // [kg] left to decompose
        float decomposed = 0;                                                         // fraction of `unburnt`
        if (unburnt > 0 && pyrolysisRate(Theated) * dt > 1e-7f) {
            auto balance = [&](float Tend, float& fraction) {
                fraction = 1.0f - std::exp(-pyrolysisRate(Tend) * dt);
                return heatCapacity * (Tend - Theated) + m.pyrolysisHeat * volatileMass * fraction;
            };
            float lo = std::max(-T0 + 1.0f, Theated - m.pyrolysisHeat * volatileMass / heatCapacity), hi = Theated;
            for (int it = 0; it < 30; ++it) {
                const float mid = 0.5f * (lo + hi);
                float f;
                (balance(mid, f) > 0 ? hi : lo) = mid;
            }
            balance(hi, decomposed);
            Theated = hi;
        }
        T = std::max(0.0f, Theated);
        if (decomposed > 0) {
            fuelToGas[k] += m.fuelYield * volatileMass * decomposed;
            unburnt *= 1.0f - decomposed;
            if (unburnt < 0.02f) unburnt = 0; // charred through
        }
    }
    int burntThrough = 0;
    const int first = cloth.firstParticle;
    for (DistanceConstraint& c : cloth.constraints) {
        if (c.broken || (c.kind != DistanceConstraint::Warp && c.kind != DistanceConstraint::Weft)) continue;
        if (cloth.unburnt[c.a - first] > 0 || cloth.unburnt[c.b - first] > 0) continue;
        breakThread(cloth, c);
        ++burntThrough;
    }
    cloth.burntThreads += burntThrough;
    return burntThrough;
}

} // namespace rf
