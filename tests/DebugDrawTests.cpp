// The research layers of the debug drawing (core/Probe.h): with every layer off a frame draws
// nothing and allocates nothing more; with a layer on, what it draws matches what the solver
// computed - one point per contact point, one box per node of the world's AABB tree (2 n - 1 for
// n bodies), one line per cloth thread, the EPA depth of a watched pair equal to the contact's
// depth; the gas and plasma layers appear on the slice; labels and contacts reach the snapshot;
// a layer stops at Probe::kLayerCap primitives with a label saying so.
#include "TestRunner.h"
#include "Tests.h"

#include "core/Probe.h"

#include <cmath>
#include <string>

namespace {

size_t countLines(const Probe::Snapshot& s, DrawLayer l) {
    size_t n = 0;
    for (const Probe::Line& line : s.lines) n += line.layer == l;
    return n;
}

size_t countPoints(const Probe::Snapshot& s, DrawLayer l) {
    size_t n = 0;
    for (const Probe::Point& p : s.points) n += p.layer == l;
    return n;
}

size_t countLabels(const Probe::Snapshot& s, DrawLayer l) {
    size_t n = 0;
    for (const Probe::Label& label : s.labels) n += label.layer == l;
    return n;
}

// Heap allocations of `frames` frames of the simulation, from the census in main.cpp.
double frameAllocations(Simulation& sim, int frames) {
    double sum = 0;
    for (int f = 0; f < frames; ++f) {
        sim.stepFrame();
        sum += Probe::snapshot().value("memory/allocations per frame");
    }
    return sum;
}

void resetProbe() {
    Probe::clear();
    Probe::setLayers(0);
    Probe::watchPair(-1, -1);
}

} // namespace

void testDebugLayersOff() {
    // The tower of 100 boxes, 40 frames, all layers off: the reference.
    resetProbe();
    Simulation plain;
    loadSample(plain, Preset::RigidTower);
    frameAllocations(plain, 30);
    const double reference = frameAllocations(plain, 10);
    CHECK(Probe::primitiveCount() == 0, "layers off: %zu primitives", Probe::primitiveCount());
    // The same tower with every layer on and a watched pair for 15 frames, then off again: from
    // then on the same frames draw nothing and allocate no more than the reference.
    Simulation drawn;
    loadSample(drawn, Preset::RigidTower);
    Probe::setLayers(~0u);
    Probe::watchPair(0, 1);
    frameAllocations(drawn, 15);
    const size_t whileOn = Probe::primitiveCount();
    Probe::setLayers(0);
    frameAllocations(drawn, 15);
    const double afterwards = frameAllocations(drawn, 10);
    std::printf("  layers on: %zu primitives a frame; off: %.0f allocations in 10 frames (reference %.0f)\n", whileOn,
                afterwards, reference);
    // Drawing only looks: the tower that was drawn stands exactly where the plain one does.
    float drift = 0;
    for (size_t i = 0; i < plain.rigid.bodies().size(); ++i)
        drift = std::max(drift, length(plain.rigid.bodies()[i].pos - drawn.rigid.bodies()[i].pos));
    std::printf("  the drawn tower vs the plain one: largest difference of a body's position %.1e m\n", double(drift));
    CHECK(whileOn > 0, "every layer on draws something");
    CHECK(Probe::primitiveCount() == 0, "layers off again: %zu primitives", Probe::primitiveCount());
    CHECK(drift == 0, "drawing changed the physics: a body moved %e m", double(drift));
    // The count of allocations is not the same bit for bit from run to run: a worker thread makes
    // its scratch buffer when it first gets work, and which worker gets work first is up to the
    // operating system (measured: 600 .. 615 against a reference of 603 .. 606). So "no more than
    // the reference" is read with that noise: 3 % of it. A layer that leaked would cost thousands.
    CHECK(afterwards <= 1.03 * reference, "layers off cost %.0f allocations more than never on", afterwards - reference);
    resetProbe();
}

void testDebugContactLayers() {
    // A stack of 6 boxes at rest: one point, one normal arrow (3 segments) per contact point.
    resetProbe();
    RigidWorld w;
    w.setDomain(AABB({-5, 0, -5}, {5, 10, 5}));
    for (int i = 0; i < 6; ++i) w.addBox({0.01f * float(i % 2), 0.25f + 0.5f * float(i), 0}, Vector3(0.25f), Quaternion(), 500, Vector3(1));
    for (int f = 0; f < 60; ++f) w.step(1.0f / 240);
    Probe::setLayers(Probe::bit(DrawLayer::ContactPoints) | Probe::bit(DrawLayer::ContactNormals));
    Probe::beginFrame();
    w.step(1.0f / 240);
    w.step(1.0f / 240); // a second substep replaces the first one's drawing
    const Probe::Snapshot s = Probe::snapshot();
    const size_t points = countPoints(s, DrawLayer::ContactPoints), normals = countLines(s, DrawLayer::ContactNormals);
    std::printf("  %zu contact points, %zu points drawn, %zu normal segments\n", w.contactCount(), points, normals);
    CHECK(w.contactCount() > 0, "the stack has contacts");
    CHECK(points == w.contactCount(), "one point per contact point: %zu vs %zu", points, w.contactCount());
    CHECK(normals == 3 * w.contactCount(), "one arrow (3 segments) per contact point: %zu", normals);
    CHECK(countLines(s, DrawLayer::BodyAabbs) == 0, "a layer that is off draws nothing");

    // The same through the Simulation: the snapshot's contact list and the drawing of the frame.
    Simulation sim;
    loadSample(sim, Preset::RigidTower);
    for (int f = 0; f < 20; ++f) sim.stepFrame();
    RenderSnapshot snap;
    sim.fillSnapshot(snap);
    CHECK(snap.contacts.size() == sim.rigid.contactCount(), "snapshot contacts %zu vs %zu", snap.contacts.size(), sim.rigid.contactCount());
    CHECK(countPoints(snap.probe, DrawLayer::ContactPoints) == sim.rigid.contactCount(), "one frame of substeps draws each point once: %zu",
          countPoints(snap.probe, DrawLayer::ContactPoints));
    for (const RenderSnapshot::ContactInfo& c : snap.contacts)
        CHECK(std::fabs(length(c.normal) - 1.0f) < 1e-3f && c.normalImpulse >= 0 && c.bodyA >= -1, "a contact as the solver has it");
    Probe::setLayers(0);
    sim.fillSnapshot(snap);
    CHECK(snap.contacts.empty(), "no contact list while ContactPoints is off: %zu", snap.contacts.size());
    resetProbe();
}

void testDebugWatchedPair() {
    // Two boxes of 1 m overlapping by 0.1 m in y, no gravity: GJK finds the overlap, EPA its depth,
    // which must equal the depth of the contact the narrow phase made (SAT for boxes).
    resetProbe();
    RigidWorld w;
    w.params.gravity = Vector3(0.0f);
    w.setDomain(AABB({-10, -10, -10}, {10, 10, 10}));
    const int a = w.addBox({0, 1, 0}, Vector3(0.5f), Quaternion(), 500, Vector3(1));
    const int b = w.addBox({0.1f, 1.9f, 0}, Vector3(0.5f), Quaternion(), 500, Vector3(1));
    Probe::setLayers(Probe::bits(DrawLayer::GjkSimplex, DrawLayer::WitnessPoints));
    Probe::watchPair(a, b);
    w.step(1.0f / 240);
    const RigidWorld::WatchReport& r = w.watchReport();
    float contactDepth = 0;
    for (const RigidWorld::DebugContact& c : w.debugContacts()) contactDepth = std::max(contactDepth, c.depth);
    const Probe::Snapshot s = Probe::snapshot();
    float arrowLength = 0; // the first WitnessPoints line is the shaft of the penetration arrow
    for (const Probe::Line& l : s.lines)
        if (l.layer == DrawLayer::WitnessPoints) { arrowLength = length(l.b - l.a); break; }
    std::printf("  GJK %d simplices, EPA %d faces, depth %.5f (contact %.5f, arrow %.5f)\n", r.gjkIterations, r.epaFaces,
                r.depth, contactDepth, arrowLength);
    CHECK(r.valid && r.intersect, "the pair overlaps");
    CHECK(r.gjkIterations >= 1, "GJK iterations %d", r.gjkIterations);
    CHECK(r.epaFaces > 0, "EPA faces %d", r.epaFaces);
    CHECK(std::fabs(r.depth - contactDepth) < 1e-4f, "EPA depth %.6f vs contact depth %.6f", r.depth, contactDepth);
    CHECK(std::fabs(arrowLength - contactDepth) < 1e-4f, "penetration arrow %.6f vs contact depth %.6f", arrowLength, contactDepth);
    CHECK(countLabels(s, DrawLayer::GjkSimplex) == 1 && countLabels(s, DrawLayer::EpaPolytope) == 1, "GJK and EPA labels");
    CHECK(countLines(s, DrawLayer::EpaPolytope) == size_t(3 * r.epaFaces), "3 edges per EPA face: %zu", countLines(s, DrawLayer::EpaPolytope));
    resetProbe();
}

void testDebugWorldTree() {
    // 7 bodies far apart: the dynamic AABB tree is a binary tree with 7 leaves, 13 nodes, 12 edges each.
    resetProbe();
    RigidWorld w;
    w.params.gravity = Vector3(0.0f);
    w.setDomain(AABB({-20, -20, -20}, {20, 20, 20}));
    const int n = 7;
    for (int i = 0; i < n; ++i) w.addBox({2.0f * float(i) - 6.0f, 1, 0}, Vector3(0.3f), Quaternion(), 500, Vector3(1));
    Probe::setLayers(Probe::bit(DrawLayer::WorldTree));
    w.step(1.0f / 240);
    const size_t boxes = countLines(Probe::snapshot(), DrawLayer::WorldTree) / 12;
    CHECK(boxes == size_t(2 * n - 1), "world tree: %zu nodes for %d bodies", boxes, n);
    resetProbe();
}

void testDebugClothTension() {
    // A curtain on a rod: one line per warp and weft thread (shear and bend are not threads).
    resetProbe();
    ParticleSystem s;
    s.reset(AABB({-1, 0, -1}, {1, 2, 1}));
    s.addCloth({-0.3f, 1.5f, 0}, {0.6f, 0, 0}, {0, -0.6f, 0}, ClothMaterial(), 16, Vector3(1));
    for (int k = 0; k < 30; ++k) s.step(1.0f / 180);
    size_t threads = 0;
    for (const DistanceConstraint& d : s.cloths()[0].constraints)
        threads += !d.broken && (d.kind == DistanceConstraint::Warp || d.kind == DistanceConstraint::Weft);
    Probe::setLayers(Probe::bit(DrawLayer::ClothTension));
    Probe::beginFrame();
    s.step(1.0f / 180);
    const size_t lines = countLines(Probe::snapshot(), DrawLayer::ClothTension);
    std::printf("  %zu threads, %zu lines\n", threads, lines);
    CHECK(threads > 0 && lines == threads, "one line per thread: %zu vs %zu", lines, threads);
    resetProbe();
}

void testDebugGasLayers() {
    // The sphere in the wind tunnel: grid, velocity, -grad p / rho and div u on the slice, drawn
    // once a frame however often the snapshot is taken; the labels come through the snapshot.
    resetProbe();
    Simulation sim;
    loadSample(sim, Preset::TunnelSphere);
    for (int f = 0; f < 5; ++f) sim.stepFrame();
    Probe::setLayers(Probe::bits(DrawLayer::GasGrid, DrawLayer::Vorticity));
    sim.stepFrame();
    RenderSnapshot s;
    sim.fillSnapshot(s);
    const size_t gradient = countLines(s.probe, DrawLayer::PressureGradient);
    std::printf("  grid %zu, velocity %zu, gradient %zu, divergence %zu, vorticity %zu, labels %zu\n",
                countLines(s.probe, DrawLayer::GasGrid), countLines(s.probe, DrawLayer::GasVelocity), gradient,
                countPoints(s.probe, DrawLayer::Divergence), countLines(s.probe, DrawLayer::Vorticity), s.probe.labels.size());
    CHECK(countLines(s.probe, DrawLayer::GasGrid) > 0, "the grid of the slice");
    CHECK(countLines(s.probe, DrawLayer::GasVelocity) > 0, "velocity arrows");
    CHECK(gradient > 0 && gradient % 3 == 0, "pressure gradient arrows: %zu segments", gradient);
    CHECK(countPoints(s.probe, DrawLayer::Divergence) > 0, "divergence points");
    CHECK(countLines(s.probe, DrawLayer::Vorticity) > 0, "vorticity arrows around the sphere");
    CHECK(countLabels(s.probe, DrawLayer::PressureGradient) == 1, "the gradient's label reaches the snapshot");
    sim.fillSnapshot(s);
    CHECK(countLines(s.probe, DrawLayer::PressureGradient) == gradient, "a second snapshot of the frame draws nothing twice");
    CHECK(s.contacts.empty(), "no rigid contacts in the tunnel");
    // The magnetosphere: field lines of B and the current density.
    Simulation plasma;
    loadSample(plasma, Preset::Magnetosphere);
    Probe::setLayers(Probe::bit(DrawLayer::FieldLinesB) | Probe::bit(DrawLayer::CurrentDensity));
    for (int f = 0; f < 3; ++f) plasma.stepFrame();
    plasma.fillSnapshot(s);
    std::printf("  field lines B %zu segments, current density %zu\n", countLines(s.probe, DrawLayer::FieldLinesB),
                countLines(s.probe, DrawLayer::CurrentDensity));
    CHECK(countLines(s.probe, DrawLayer::FieldLinesB) > 0, "field lines of B");
    CHECK(countLabels(s.probe, DrawLayer::FieldLinesB) == 1, "|B| label");
    resetProbe();
}

void testDebugLayerCap() {
    // A layer stops at kLayerCap primitives a frame and says so once; another layer is not affected.
    resetProbe();
    Probe::setLayers(Probe::bit(DrawLayer::Misc) | Probe::bit(DrawLayer::Velocities));
    for (int i = 0; i < Probe::kLayerCap + 100; ++i) Probe::point(Vector3(float(i), 0, 0), Vector3(1.0f));
    Probe::point(DrawLayer::Velocities, Vector3(0.0f), Vector3(1.0f));
    Probe::label(DrawLayer::Velocities, Vector3(0.0f), "v");
    Probe::point(DrawLayer::Islands, Vector3(0.0f), Vector3(1.0f)); // off: not stored
    const Probe::Snapshot s = Probe::snapshot();
    CHECK(countPoints(s, DrawLayer::Misc) == size_t(Probe::kLayerCap), "capped at %d: %zu", Probe::kLayerCap, countPoints(s, DrawLayer::Misc));
    CHECK(countLabels(s, DrawLayer::Misc) == 1 && s.labels[0].text.find("обрезано") != std::string::npos, "one cut-off label");
    CHECK(countPoints(s, DrawLayer::Velocities) == 1 && countLabels(s, DrawLayer::Velocities) == 1, "another layer goes on");
    CHECK(countPoints(s, DrawLayer::Islands) == 0, "a layer that is off stores nothing");
    Probe::clearLayers(Probe::bit(DrawLayer::Misc));
    CHECK(Probe::primitiveCount() == 2, "clearLayers removes only its layers: %zu left", Probe::primitiveCount());
    resetProbe();
}
