// Standalone profiler for the contact-only torus sample. Frame mode follows Simulation;
// substep mode replays its rigid steps, observing raw convex contacts and the solver contacts.
// Build: c++ -std=c++17 -O3 -DNDEBUG -march=native -Isrc -I. tools/TorusProfile.cpp
//        build-core/librfsamples.a build-core/librfcore.a -pthread -o /tmp/rf_torus_profile
#include "samples/Samples.h"
#include "samples/TorusChainsScene.h"
#include "core/Format.h"
#include "core/Probe.h"
#include "rigid/NarrowPhase.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

using namespace rf;
namespace {

struct Options {
    int frames = 120, steps = 10, iterations = 6, pairA = 67, pairB = 68, links = 40;
    float ccdThreshold = 0.5f, load = 10, kick = 1, cfm = 0.001f, slop = 0.004f;
    bool substep = false, shock = true, lock = true, warm = true, ccd = true, block = true;
    std::filesystem::path output;
};

Options options(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        const std::string key = argv[i];
        if (++i == argc) throw std::runtime_error("missing option value");
        const std::string value = argv[i];
        if (key == "--output") o.output = value;
        else if (key == "--mode") {
            if (value != "frame" && value != "substep") throw std::runtime_error("mode must be frame or substep");
            o.substep = value == "substep";
        }
        else {
            float v = 0;
            if (!parseNumber(value, v) || !std::isfinite(v)) throw std::runtime_error("invalid number");
            if (key == "--frames" || key == "--substeps" || key == "--iterations" ||
                key == "--pair-a" || key == "--pair-b" || key == "--links")
                if (v < 0 || v > 100000 || std::floor(v) != v) throw std::runtime_error("invalid integer option");
            if (key == "--frames") o.frames = int(v);
            else if (key == "--substeps") o.steps = int(v);
            else if (key == "--iterations") o.iterations = int(v);
            else if (key == "--pair-a") o.pairA = int(v);
            else if (key == "--pair-b") o.pairB = int(v);
            else if (key == "--links") o.links = int(v);
            else if (key == "--ccd-threshold") o.ccdThreshold = float(v);
            else if (key == "--cfm") o.cfm = v;
            else if (key == "--slop") o.slop = v;
            else if (key == "--load") o.load = float(v);
            else if (key == "--kick") o.kick = float(v);
            else if (key == "--shock") o.shock = v != 0;
            else if (key == "--lock") o.lock = v != 0;
            else if (key == "--warm") o.warm = v != 0;
            else if (key == "--ccd") o.ccd = v != 0;
            else if (key == "--block") o.block = v != 0;
            else throw std::runtime_error("unknown option: " + key);
        }
    }
    if (o.output.empty() || o.frames < 1 || o.steps < 1 || o.iterations < 1 ||
        o.links < 8 || o.links > 100 || o.pairA < 0 || o.pairB < 0 ||
        o.pairA >= 3 * o.links || o.pairB >= 3 * o.links || o.pairA == o.pairB ||
        o.ccdThreshold < 0 || o.cfm < 0 || o.slop < 0 || o.load < 1 || o.load > 50 || o.kick < 0 || o.kick > 5)
        throw std::runtime_error("invalid options; provide --output <new-directory>");
    if (std::filesystem::exists(o.output)) throw std::runtime_error("output already exists");
    return o;
}

struct Csv {
    std::ofstream stream;
    bool first = true;
    explicit Csv(const std::filesystem::path& path) : stream(path) {
        if (!stream) throw std::runtime_error("cannot write CSV");
    }
    void add(double value) { if (!first) stream << ','; stream << numberText(value); first = false; }
    void end() { stream << '\n'; first = true; }
};

std::vector<PosedShape> partsOf(const RigidBody& body) {
    std::vector<PosedShape> parts;
    const auto& compound = static_cast<const CompoundShape&>(*body.shape);
    const Matrix3x3 rotation = body.rotation();
    for (const auto& c : compound.children())
        parts.push_back({c.shape.get(), rotation * c.R, body.pos + rotation * c.t});
    return parts;
}

bool otherContains(const std::vector<PosedShape>& parts, size_t skip, const Vector3& point) {
    for (size_t i = 0; i < parts.size(); ++i)
        if (i != skip && parts[i].shape->contains(parts[i].R.transposed() * (point - parts[i].p))) return true;
    return false;
}

struct Raw {
    int pairs = 0, rejected = 0, retained = 0, points = 0, safePoints = 0, lostSafePoints = 0;
    float depth = -kInf, retainedDepth = -kInf;
};

// Replay the current seam predicate against raw narrow-phase contacts. This is an observation,
// not a substitute solver: the actual contact count/impulses are read independently below.
Raw rawPair(const RigidBody& a, const RigidBody& b, Csv* details, int step) {
    Raw raw;
    const auto pa = partsOf(a), pb = partsOf(b);
    NarrowPhase narrow;
    for (size_t u = 0; u < pa.size(); ++u)
        for (size_t v = 0; v < pb.size(); ++v) {
            ContactManifold manifold;
            if (!narrow.collide(pa[u], pb[v], manifold) || manifold.points.empty()) continue;
            const auto& p = manifold.points.front();
            const Vector3 onA = p.position - p.normal * (0.5f * p.depth), onB = p.position + p.normal * (0.5f * p.depth);
            const bool seamA = otherContains(pa, u, onA - p.normal * 0.0005f);
            const bool seamB = otherContains(pb, v, onB + p.normal * 0.0005f);
            ++raw.pairs;
            raw.points += int(manifold.points.size());
            raw.rejected += seamA || seamB;
            raw.retained += !(seamA || seamB);
            int safe = 0;
            for (const auto& c : manifold.points) {
                const Vector3 surfaceA = c.position - c.normal * (0.5f * c.depth);
                const Vector3 surfaceB = c.position + c.normal * (0.5f * c.depth);
                safe += !otherContains(pa, u, surfaceA - c.normal * 0.0005f) &&
                        !otherContains(pb, v, surfaceB + c.normal * 0.0005f);
            }
            raw.safePoints += safe;
            if (seamA || seamB) raw.lostSafePoints += safe;
            for (const auto& c : manifold.points) raw.depth = std::max(raw.depth, c.depth);
            if (!seamA && !seamB)
                for (const auto& c : manifold.points) raw.retainedDepth = std::max(raw.retainedDepth, c.depth);
            if (details) {
                const double values[] = {double(step), double(u), double(v), double(manifold.points.size()),
                    p.depth, double(seamA), double(seamB), p.normal.x, p.normal.y, p.normal.z, double(safe)};
                for (double value : values) details->add(value);
                details->end();
            }
        }
    return raw;
}

struct Topology {
    int broken[3] = {}, first[3] = {-1, -1, -1};
    std::vector<int> initial;
    void observe(const RigidWorld& world, int count, int step) {
        size_t pair = 0;
        for (int chain = 0; chain < 3; ++chain) {
            broken[chain] = 0;
            for (int k = 1; k < count; ++k, ++pair) {
                const int a = chain * count + k - 1, b = a + 1;
                const int linking = TorusChainsScene::linkingNumber(world.bodies()[size_t(a)], world.bodies()[size_t(b)]);
                if (step == 0) initial.push_back(linking);
                else if (linking != initial[pair] || std::abs(linking) != 1) {
                    ++broken[chain];
                    if (first[chain] < 0) std::printf("first broken chain %d, pair %d-%d, observed substep %d\n", chain, a, b, step);
                }
            }
            if (broken[chain] && first[chain] < 0) first[chain] = step;
        }
    }
};

uint64_t stateHash(const Simulation& sim) {
    uint64_t hash = 14695981039346656037ull;
    auto add = [&](float value) {
        const auto* bytes = reinterpret_cast<const unsigned char*>(&value);
        for (size_t i = 0; i < sizeof(value); ++i) hash = (hash ^ bytes[i]) * 1099511628211ull;
    };
    for (const auto& b : sim.rigid.bodies()) {
        for (float v : {b.pos.x, b.pos.y, b.pos.z, b.rot.x, b.rot.y, b.rot.z, b.rot.w,
                        b.vel.x, b.vel.y, b.vel.z, b.angVel.x, b.angVel.y, b.angVel.z}) add(v);
    }
    return hash;
}

void configure(Simulation& sim, const Options& o) {
    loadSample(sim, Preset::TorusChains);
    sim.scene()->setParam(0, float(o.links)); sim.scene()->setParam(1, o.load); sim.scene()->setParam(2, o.kick);
    sim.reset();
    auto& p = sim.rigid.params;
    p.substeps = o.steps; p.iterations = o.iterations; p.ccdThreshold = o.ccdThreshold;
    p.shockPropagation = o.shock; p.rotationalLock = o.lock; p.warmStarting = o.warm; p.ccd = o.ccd;
    p.blockCfm = o.cfm; p.slop = o.slop; p.blockSolver = o.block;
}

void pairRow(Csv& csv, Csv& details, Csv& poses, Simulation& sim, const Options& o, int step) {
    const auto &a = sim.rigid.bodies()[size_t(o.pairA)], &b = sim.rigid.bodies()[size_t(o.pairB)];
    const Raw raw = rawPair(a, b, step >= 400 && step <= 510 ? &details : nullptr, step);
    int contacts = 0;
    float impulse = 0, deepest = -kInf, residual = 0;
    for (const auto& c : sim.rigid.debugContacts())
        if ((c.a == o.pairA && c.b == o.pairB) || (c.a == o.pairB && c.b == o.pairA)) {
            ++contacts;
            impulse += c.impulse;
            deepest = std::max(deepest, c.depth);
            const auto &A = sim.rigid.bodies()[size_t(c.a)], &B = sim.rigid.bodies()[size_t(c.b)];
            residual = std::min(residual, dot(A.velocityAt(c.position) - B.velocityAt(c.position), c.normal));
        }
    const double values[] = {double(step), double(step) / (60 * o.steps), double(TorusChainsScene::linkingNumber(a, b)),
        length(a.pos - b.pos), double(raw.pairs), double(raw.rejected), double(raw.retained), raw.depth, raw.retainedDepth,
        double(contacts), deepest, impulse, residual, length(a.vel), length(b.vel), length(a.angVel), length(b.angVel),
        double(sim.rigid.ccdHits()), double(raw.safePoints), double(raw.lostSafePoints)};
    for (double v : values) csv.add(std::isfinite(v) ? v : 0);
    csv.end();
    poses.add(step);
    for (const auto* body : {&a, &b})
        for (float value : {body->pos.x, body->pos.y, body->pos.z, body->rot.w, body->rot.x, body->rot.y, body->rot.z,
                           body->vel.x, body->vel.y, body->vel.z, body->angVel.x, body->angVel.y, body->angVel.z})
            poses.add(value);
    poses.end();
}

double mechanicalEnergy(const RigidWorld& world) {
    double energy = world.kineticEnergy();
    for (const auto& body : world.bodies())
        if (body.invMass > 0) energy -= double(body.mass) * dot(world.params.gravity, body.pos);
    return energy;
}

void writeBodyPoses(Csv& csv, const RigidWorld& world, int frame) {
    for (size_t i = 0; i < world.bodies().size(); ++i) {
        const auto& b = world.bodies()[i];
        csv.add(frame); csv.add(double(i));
        for (float value : {b.pos.x, b.pos.y, b.pos.z, b.rot.w, b.rot.x, b.rot.y, b.rot.z,
                            b.vel.x, b.vel.y, b.vel.z, b.angVel.x, b.angVel.y, b.angVel.z}) csv.add(value);
        csv.end();
    }
}

int profile(const Options& o) {
    std::filesystem::create_directories(o.output);
    Simulation sim;
    configure(sim, o);
    Topology topology;
    topology.observe(sim.rigid, o.links, 0);
    Csv frames(o.output / "frames.csv"), pairs(o.output / "pair.csv"), details(o.output / "raw-parts.csv");
    Csv poses(o.output / "pair-poses.csv");
    Csv bodies(o.output / "body-poses.csv");
    const double initialEnergy = mechanicalEnergy(sim.rigid);
    double maxEnergyGain = 0;
    bodies.stream << "frame,body,x,y,z,qw,qx,qy,qz,vx,vy,vz,wx,wy,wz\n";
    frames.stream << "frame,time_s,wall_ms,collide_ms,solve_ms,ccd_ms,integrate_ms,islands_ms,controllers_ms,contacts,manifolds,pairs,ccd_hits,broken_hanging,broken_released,broken_bridge,toi_body_pairs,toi_part_pairs,toi_gjk_calls,toi_initial_overlap,toi_initial_near,toi_iteration_limit,energy_j,energy_gain_j\n";
    pairs.stream << "substep,time_s,linking,distance_m,raw_pairs,rejected_pairs,retained_pairs,raw_max_depth_m,retained_max_depth_m,solver_contacts,solver_max_depth_m,normal_impulse_ns,min_residual_vn_m_s,speed_a,speed_b,spin_a,spin_b,ccd_hits_last_substep,safe_points,lost_safe_points\n";
    details.stream << "substep,part_a,part_b,points,depth_m,seam_a,seam_b,nx,ny,nz,safe_points\n";
    poses.stream << "substep,ax,ay,az,aqw,aqx,aqy,aqz,avx,avy,avz,awx,awy,awz,bx,by,bz,bqw,bqx,bqy,bqz,bvx,bvy,bvz,bwx,bwy,bwz\n";
    const float dt = sim.frameDt / float(o.steps);
    for (int frame = 1; frame <= o.frames; ++frame) {
        double wall = 0;
        if (!o.substep) {
            auto start = std::chrono::steady_clock::now();
            sim.stepFrame();
            wall = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
            topology.observe(sim.rigid, o.links, frame * o.steps);
        } else {
            Probe::beginFrame();
            for (int sub = 1; sub <= o.steps; ++sub) {
                auto start = std::chrono::steady_clock::now();
                sim.rigid.step(dt);
                wall += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
                const int step = (frame - 1) * o.steps + sub;
                topology.observe(sim.rigid, o.links, step);
                pairRow(pairs, details, poses, sim, o, step);
            }
        }
        const double energy = mechanicalEnergy(sim.rigid);
        maxEnergyGain = std::max(maxEnergyGain, energy - initialEnergy);
        writeBodyPoses(bodies, sim.rigid, frame);
        const auto p = Probe::snapshot();
        const double values[] = {double(frame), double(frame) / 60, wall, p.value("rigid/collide ms"), p.value("rigid/solve ms"),
            p.value("rigid/ccd ms"), p.value("rigid/integrate ms"), p.value("rigid/islands ms"), p.value("scene/controllers ms"),
            double(sim.rigid.contactCount()), double(sim.rigid.manifoldCount()), double(sim.rigid.pairCount()), p.value("rigid/ccd hits"),
            double(topology.broken[0]), double(topology.broken[1]), double(topology.broken[2]),
            p.value("toi/body pairs", -1), p.value("toi/part pairs", -1), p.value("toi/gjk calls", -1),
            p.value("toi/initial overlap", -1), p.value("toi/initial near", -1), p.value("toi/iteration limit", -1), energy, maxEnergyGain};
        for (double value : values) frames.add(value);
        frames.end();
        frames.stream.flush(); bodies.stream.flush();
        if (frame % 30 == 0) std::printf("frame %d: broken %d/%d/%d\n", frame, topology.broken[0], topology.broken[1], topology.broken[2]);
    }
    std::ofstream summary(o.output / "summary.txt");
    summary << "mode " << (o.substep ? "substep" : "Simulation") << "\nframes " << o.frames << "\nsubsteps " << o.steps
        << "\niterations " << o.iterations << "\nccdThreshold " << numberText(o.ccdThreshold) << "\nlinks " << o.links
        << "\nload " << numberText(o.load) << "\nkick " << numberText(o.kick) << "\nstateHash " << std::hex << stateHash(sim) << std::dec << '\n';
    summary << "pairA " << o.pairA << "\npairB " << o.pairB << "\nshock " << o.shock << "\nlock " << o.lock << "\nwarm " << o.warm << "\nccd " << o.ccd
        << "\nblock " << o.block << "\ncfm " << numberText(o.cfm) << "\nslop " << numberText(o.slop) << '\n';
    summary << "initialEnergyJ " << numberText(initialEnergy) << "\nmaxEnergyGainJ " << numberText(maxEnergyGain) << '\n';
    for (int c = 0; c < 3; ++c) summary << "chain " << c << " firstBrokenSubstep " << topology.first[c]
        << " brokenFinal " << topology.broken[c] << '\n';
    return 0; // this is a measurement tool, not a passing physics acceptance test
}

} // namespace
int main(int argc, char** argv) {
    try { return profile(options(argc, argv)); }
    catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 2; }
}
