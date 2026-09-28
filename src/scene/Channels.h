#pragma once
// Channels: the physical quantities a researcher plots - of the whole scene, of one object, of a
// group of objects and of the fields at a point or along a line. «Графики — это система
// исследований всего и вся, как сцены, так и частных элементов».
//
// A channel is a number with a name card: a stable ASCII id ("rigid/kinetic energy", the Probe key
// and the CSV column), a unit ("Дж") and a name for people ("Твёрдые тела: энергия движения").
// The same functions give the numbers to the window, to the tests and to rf_verify, so what a
// chart shows is exactly what a test checks.
//
// Nothing is measured unless asked. Simulation::channels says what to measure after every frame
// (the scene, some objects, some groups, some probes); while it is empty the frame does not touch
// any of this - no time, no memory - and measuring never changes the physics (it only reads).
//
//   sim.channels.scene = true;                                      the scene's energy, momentum ...
//   sim.channels.objects.push_back({ObjectRef::Kind::Body, 3, "obj3"}); body 3: speed, contact force ...
//   sim.channels.probes.push_back({"probe1", a, b, 11, FieldQuantity::VelocityX}); a profile along a..b
//   sim.stepFrame();
//   for (const Measurement& m : sim.measurements()) ...            every value with its unit
//
// Every value is also reported to the Probe under its id, so the charts and the CSV export pick it
// up by name. Tests call measureScene / measureObject / measureGroup / measureProbe directly.
#include "math/Vector3.h"

#include <string>
#include <vector>

namespace rf {

class Simulation;

// Where a channel belongs in the viewer's list of channels.
enum class ChannelGroup { Scene, Object, Group, Field };

// The name card of a channel.
struct ChannelInfo {
    std::string id;   // ASCII "part/quantity": the Probe key and the CSV column
    std::string unit; // "Дж", "м/с", "Н", ... (UTF-8)
    std::string name; // for people, in Russian: «Твёрдые тела: энергия движения»
    ChannelGroup group = ChannelGroup::Scene;
};

// One measured value - or, for a probe line, the profile of one quantity along it.
struct Measurement {
    ChannelInfo info;
    double value = 0;              // the number (for a profile: its mean)
    std::vector<double> profile;   // probe line: the values along it ...
    std::vector<double> distance;  // ... at these distances from its start [m]
};

// A thing in the scene the viewer points at: a rigid body (index in RigidWorld::bodies()) or a
// particle group (a soft body, a cloth or a block of liquid: ParticleSystem::groupOf).
struct ObjectRef {
    enum class Kind { Body, ParticleGroup };
    Kind kind = Kind::Body;
    int index = -1;
    std::string label; // ASCII prefix of its channel ids: "obj3" -> "obj3/speed"
};

// Several objects measured as one: their total mass, centre of mass, momentum and energy.
struct GroupRef {
    std::string label; // ASCII prefix: "group1" -> "group1/momentum"
    std::vector<ObjectRef> members;
};

// What a probe in the fields can read.
enum class FieldQuantity {
    VelocityX, VelocityY, VelocityZ, Speed, // the gas velocity [m/s]
    Pressure,                               // the gas pressure [Pa] (defined up to a constant in a closed box)
    Temperature,                            // above the ambient temperature [K]
    Smoke,                                  // smoke density (0 .. 1)
    MagneticX, MagneticY, MagneticZ, MagneticMagnitude, // B [T]
    CurrentDensity                          // |J| [A/m^2]
};

// A probe in the gas and the magnetic field. One point (points = 1, `to` unused): every quantity
// at it. A line of `points` points from `from` to `to`: the profile of `quantity` along it - a
// velocity profile across a channel to lay beside Poiseuille's parabola.
struct FieldProbe {
    std::string label; // ASCII prefix: "probe1"
    Vector3 from{0.0f}, to{0.0f};
    int points = 1;
    FieldQuantity quantity = FieldQuantity::Speed;
};

// What to measure after every frame. Empty: nothing is measured.
struct ChannelRequests {
    bool scene = false;
    std::vector<ObjectRef> objects;
    std::vector<GroupRef> groups;
    std::vector<FieldProbe> probes;
    bool empty() const { return !scene && objects.empty() && groups.empty() && probes.empty(); }
};

// --- Measuring (read only: the physics never changes) --------------------------------------------
// The scene: energy, momentum, angular momentum of the bodies; energy of the particles, density
// error of the liquid; energy, divergence and speed of the gas; the magnetic energy. Only what the
// scene has (no gas: no gas channels).
std::vector<Measurement> measureScene(const Simulation& sim);
// One object: where it is, how fast, its energy; a body also the contact force on it.
std::vector<Measurement> measureObject(const Simulation& sim, const ObjectRef& object);
// A group of objects: total mass, centre of mass, momentum and energy.
std::vector<Measurement> measureGroup(const Simulation& sim, const GroupRef& group);
// A probe in the fields (only while the scene has a gas grid).
std::vector<Measurement> measureProbe(const Simulation& sim, const FieldProbe& probe);
// Everything `requests` asks for, in order: scene, objects, groups, probes.
std::vector<Measurement> measureAll(const Simulation& sim, const ChannelRequests& requests);

// The name cards of the scene channels (for a viewer's list before anything was measured), and the
// card of any scene channel or engine channel by its id (nullptr: not a known physical channel).
const std::vector<ChannelInfo>& sceneChannelCatalog();
const ChannelInfo* describeChannel(const std::string& id);
// The unit and the Russian name of a probe quantity.
ChannelInfo fieldQuantityInfo(FieldQuantity q);

} // namespace rf
