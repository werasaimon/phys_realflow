// Measure settling of the editor's eight-torus chain, including an identical mouse excitation.
// U = -sum(m g.x); K is the world's translational plus rotational kinetic energy. Static
// Coulomb friction permits equilibria above the frictionless minimum (Drake friction model:
// https://drake.mit.edu/doxygen_cxx/group__friction__model.html). This is an experiment, not
// a global minimizer. Branches restore poses/velocities with cold contact caches, including control.
// Build: c++ -std=c++17 -O3 -DNDEBUG -march=native -Isrc -I. tools/ChainRelaxation.cpp
//        build-core/librfsamples.a build-core/librfcore.a -pthread -o /tmp/rf_chain_relaxation
// Run: RF_THREADS=1 /tmp/rf_chain_relaxation default 3600 NEW_DIR shake
// Args: mode, settling frames, NEW_DIR, shake|straight|SNAPSHOT_CSV. Zero frames inspects a snapshot.
// No solver settings are saved globally.
#include "samples/Samples.h"
#include "samples/HangingTorusScene.h"
#include "samples/TorusChainsScene.h"
#include "core/Format.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace rf;
namespace {

void row(std::ostream& out, std::initializer_list<double> values) {
    bool first = true;
    for (double value : values) {
        if (!std::isfinite(value)) throw std::runtime_error("non-finite observation");
        if (!first) out << ',';
        out << numberText(value);
        first = false;
    }
    out << '\n';
}

void settings(RigidWorld& w, const std::string& mode) {
    if (mode == "default") return;
    if (mode == "no-roll") w.params.rollingResistance = 0;
    else if (mode == "no-lock") w.params.rotationalLock = false;
    else if (mode == "no-shock") w.params.shockPropagation = false;
    else if (mode == "fine") { w.params.substeps = 80; w.params.iterations = 24; }
    else if (mode == "no-coulomb") {
        for (auto& b : w.bodies()) b.friction = b.staticFriction = 0;
    } else if (mode == "free" || mode == "free-damped" || mode == "free-fine") {
        for (auto& b : w.bodies()) b.friction = b.staticFriction = 0;
        w.params.rollingResistance = 0;
        w.params.rotationalLock = false;
        if (mode == "free-damped" || mode == "free-fine")
            w.params.linearDamping = w.params.angularDamping = 1;
        if (mode == "free-fine") { w.params.substeps = 80; w.params.iterations = 24; }
    } else throw std::runtime_error("unknown mode");
}

Vector3 rim(const RigidBody& body) {
    const auto& shape = static_cast<const CompoundShape&>(*body.shape);
    const Matrix3x3 r = body.rotation() * shape.principalRotation().transposed();
    return body.pos + r * (Vector3(0, 0, TorusChainsScene::majorRadius) - shape.centerOfMass());
}

void mouse(Simulation& sim, int f, Vector3& start) {
    auto& w = sim.rigid;
    if (f == 60) {
        const Vector3 direction = normalize(Vector3(1, 0, 1)), origin = rim(w.bodies().back()) + 2 * direction;
        int hit = -1; float distance = 0; Vector3 normal;
        if (!w.raycast(origin, -direction, 3, hit, distance, normal) || hit != 7)
            throw std::runtime_error("mouse missed lower rim");
        start = origin - direction * distance;
        w.grab(hit, start);
    }
    if (f >= 60 && f < 240) {
        const float phase = 2 * kPi * float(f - 60) / 72;
        const Vector3 offset = f >= 150 && f < 180 ? Vector3(-0.22f, 0.18f, -0.08f)
            : Vector3(0.18f * rf::sin(phase), 0.12f * (1 - rf::cos(phase)), 0.08f * rf::sin(0.5f * phase));
        w.setGrabTarget(start + offset);
    }
}

void poses(std::ostream& out, const RigidWorld& w, int frame) {
    for (size_t i = 0; i < w.bodies().size(); ++i) {
        const auto& b = w.bodies()[i];
        row(out, {double(frame), double(i), b.pos.x, b.pos.y, b.pos.z, b.rot.w, b.rot.x, b.rot.y, b.rot.z,
            b.vel.x, b.vel.y, b.vel.z, b.angVel.x, b.angVel.y, b.angVel.z});
    }
}

void restore(RigidWorld& w, const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open snapshot");
    std::string line;
    std::getline(input, line);
    size_t index = 0;
    while (std::getline(input, line)) {
        std::istringstream words(line); std::string word; std::vector<float> v;
        while (std::getline(words, word, ',')) {
            float x = 0;
            if (!parseNumber(word, x) || !std::isfinite(x)) throw std::runtime_error("invalid snapshot number");
            v.push_back(x);
        }
        if (v.size() != 15 || index >= w.bodies().size() || v[1] != float(index))
            throw std::runtime_error("invalid snapshot row");
        auto& b = w.bodies()[index++];
        b.pos = b.prevPos = {v[2], v[3], v[4]}; b.rot = b.prevRot = {v[5], v[6], v[7], v[8]};
        b.vel = {v[9], v[10], v[11]}; b.angVel = {v[12], v[13], v[14]};
        b.updateInertia();
    }
    if (index != w.bodies().size()) throw std::runtime_error("incomplete snapshot");
}

void observe(std::ostream& out, const Simulation& sim, int frame, double time) {
    const auto& w = sim.rigid;
    double u = 0, mass = 0, cx = 0, cy = 0, cz = 0, maxSpeed = 0, maxSpin = 0, radial2 = 0;
    int broken = 0;
    for (size_t i = 0; i < w.bodies().size(); ++i) {
        const auto& b = w.bodies()[i];
        if (i) broken += std::abs(TorusChainsScene::linkingNumber(w.bodies()[i - 1], b)) != 1;
        if (b.invMass == 0) continue;
        u -= double(b.mass) * dot(w.params.gravity, b.pos); mass += b.mass;
        cx += double(b.mass) * b.pos.x; cy += double(b.mass) * b.pos.y; cz += double(b.mass) * b.pos.z;
        radial2 += double(b.mass) * (b.pos.x * b.pos.x + b.pos.z * b.pos.z);
        maxSpeed = std::max(maxSpeed, double(length(b.vel))); maxSpin = std::max(maxSpin, double(length(b.angVel)));
    }
    if (broken) throw std::runtime_error("chain link broken");
    const auto& tip = w.bodies().back();
    row(out, {double(frame), time, w.kineticEnergy(), u, mass, cx / mass, cy / mass, cz / mass,
        std::sqrt(radial2 / mass), tip.pos.x, tip.pos.y, tip.pos.z, maxSpeed, maxSpin,
        w.deepestPenetration(), double(broken), double(w.ccdHits()), double(w.sleepingCount())});
}

void experiment(Simulation& sim, int frames, const std::filesystem::path& dir, const std::string& initial) {
    std::ofstream csv(dir / "energy.csv"), body(dir / "poses.csv"), state(dir / "final.csv"), release(dir / "release.csv");
    if (!csv || !body || !state || !release) throw std::runtime_error("cannot write output");
    csv << "frame,time_s,kinetic_j,potential_j,mass_kg,com_x,com_y,com_z,radial_rms_m,tip_x,tip_y,tip_z,max_speed_mps,max_spin_radps,contact_depth_m,broken,ccd_hits_last_substep,sleepers\n";
    const char* header = "frame,body,x,y,z,qw,qx,qy,qz,vx,vy,vz,wx,wy,wz\n";
    body << header; state << header; release << header;
    const int prep = initial == "shake" ? 240 : 0;
    Vector3 start;
    observe(csv, sim, 0, -double(prep) / 60);
    poses(body, sim.rigid, 0);
    for (int f = 0; f < prep + frames; ++f) {
        if (f < prep) mouse(sim, f, start);
        if (f == prep) { sim.rigid.releaseGrab(); poses(release, sim.rigid, f); }
        sim.stepFrame();
        observe(csv, sim, f + 1, double(f + 1 - prep) / 60);
        if ((f + 1) % 60 == 0 || f + 1 == prep + frames) poses(body, sim.rigid, f + 1);
        if ((f + 1) % 600 == 0) { csv.flush(); body.flush(); }
    }
    poses(state, sim.rigid, prep + frames);
    if (!csv || !body || !state || !release) throw std::runtime_error("failed writing output");
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 5) throw std::runtime_error("MODE FRAMES NEW_DIR shake|straight|SNAPSHOT_CSV");
        float n = 0;
        if (!parseNumber(argv[2], n) || !std::isfinite(n) || n < 0 || n > 100000 || n != std::floor(n))
            throw std::runtime_error("invalid frame count");
        Simulation sim;
        loadSample(sim, Preset::HangingTorus);
        const std::string initial = argv[4];
        if (initial != "shake" && initial != "straight") restore(sim.rigid, initial);
        settings(sim.rigid, argv[1]);
        if (!std::filesystem::create_directory(argv[3])) throw std::runtime_error("output directory exists");
        experiment(sim, int(n), argv[3], initial);
        std::cout << "Completed " << argv[1] << ": " << argv[3] << '\n';
        return 0;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
