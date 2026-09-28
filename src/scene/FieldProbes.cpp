// Probes in the fields (see Channels.h): what the gas and the magnetic field are at a point, or
// along a line - the profile a researcher lays beside a formula (the velocity across a channel
// beside Poiseuille's parabola, the pressure down a column beside ρ g h).
//
// The gas state is read from the gas cells only (GasSolver::fluidVelocityAt / fluidPressureAt):
// right at a wall the plain interpolation would mix in the solid. A point exactly at a cell's
// centre reads that cell's own value - no interpolation error at all.
#include "scene/ChannelMath.h"
#include "scene/Simulation.h"

namespace rf {

// The name card of a probe quantity: the id suffix ("velocity x"), the unit and the name.
ChannelInfo fieldQuantityInfo(FieldQuantity q) {
    switch (q) {
    case FieldQuantity::VelocityX: return {"velocity x", "м/с", "Скорость газа, x", ChannelGroup::Field};
    case FieldQuantity::VelocityY: return {"velocity y", "м/с", "Скорость газа, y", ChannelGroup::Field};
    case FieldQuantity::VelocityZ: return {"velocity z", "м/с", "Скорость газа, z", ChannelGroup::Field};
    case FieldQuantity::Speed: return {"speed", "м/с", "Скорость газа |u|", ChannelGroup::Field};
    case FieldQuantity::Pressure: return {"pressure", "Па", "Давление газа", ChannelGroup::Field};
    case FieldQuantity::Temperature: return {"temperature", "К", "Температура над окружающей", ChannelGroup::Field};
    case FieldQuantity::Smoke: return {"smoke", "доля", "Плотность дыма", ChannelGroup::Field};
    case FieldQuantity::MagneticX: return {"B x", "Тл", "Магнитное поле, x", ChannelGroup::Field};
    case FieldQuantity::MagneticY: return {"B y", "Тл", "Магнитное поле, y", ChannelGroup::Field};
    case FieldQuantity::MagneticZ: return {"B z", "Тл", "Магнитное поле, z", ChannelGroup::Field};
    case FieldQuantity::MagneticMagnitude: return {"B", "Тл", "Магнитное поле |B|", ChannelGroup::Field};
    case FieldQuantity::CurrentDensity: return {"J", "А/м²", "Плотность тока |J|", ChannelGroup::Field};
    }
    return {};
}

namespace {

bool isMagnetic(FieldQuantity q) {
    return q == FieldQuantity::MagneticX || q == FieldQuantity::MagneticY || q == FieldQuantity::MagneticZ ||
           q == FieldQuantity::MagneticMagnitude || q == FieldQuantity::CurrentDensity;
}

// One quantity at one point of the world.
double sample(const Simulation& sim, FieldQuantity q, const Vector3& x) {
    const GasSolver& g = sim.grid;
    switch (q) {
    case FieldQuantity::VelocityX: return g.fluidVelocityAt(x).x;
    case FieldQuantity::VelocityY: return g.fluidVelocityAt(x).y;
    case FieldQuantity::VelocityZ: return g.fluidVelocityAt(x).z;
    case FieldQuantity::Speed: return length(g.fluidVelocityAt(x));
    case FieldQuantity::Pressure: return g.fluidPressureAt(x);
    case FieldQuantity::Temperature: return g.temperatureAt(x);
    case FieldQuantity::Smoke: return g.smoke().sample((x - g.origin()) / g.dx());
    case FieldQuantity::MagneticX: return g.magnetic.fieldAt(x).x;
    case FieldQuantity::MagneticY: return g.magnetic.fieldAt(x).y;
    case FieldQuantity::MagneticZ: return g.magnetic.fieldAt(x).z;
    case FieldQuantity::MagneticMagnitude: return length(g.magnetic.fieldAt(x));
    case FieldQuantity::CurrentDensity: return length(g.magnetic.currentAt(x));
    }
    return 0;
}

// A point probe: every quantity the scene has at the point.
void measurePoint(const Simulation& sim, const FieldProbe& p, std::vector<Measurement>& out) {
    for (int k = 0; k <= int(FieldQuantity::CurrentDensity); ++k) {
        const FieldQuantity q = FieldQuantity(k);
        if (isMagnetic(q) && !sim.grid.magnetic.enabled) continue;
        const ChannelInfo card = fieldQuantityInfo(q);
        out.push_back(measured(p.label + "/" + card.id, card.unit.c_str(), card.name, ChannelGroup::Field, sample(sim, q, p.from)));
    }
}

// A line probe: `points` samples evenly from `from` to `to` (both ends included); the value is
// their mean, the profile the samples with their distance from `from`.
void measureLine(const Simulation& sim, const FieldProbe& p, std::vector<Measurement>& out) {
    const ChannelInfo card = fieldQuantityInfo(p.quantity);
    Measurement m = measured(p.label + "/" + card.id + " profile", card.unit.c_str(), card.name + " вдоль линии",
                             ChannelGroup::Field, 0);
    const Vector3 step = (p.to - p.from) / float(p.points - 1);
    double sum = 0;
    for (int k = 0; k < p.points; ++k) {
        const double v = sample(sim, p.quantity, p.from + step * float(k));
        m.profile.push_back(v);
        m.distance.push_back(double(length(step)) * k);
        sum += v;
    }
    m.value = sum / p.points;
    out.push_back(std::move(m));
}

} // namespace

std::vector<Measurement> measureProbe(const Simulation& sim, const FieldProbe& probe) {
    std::vector<Measurement> out;
    if (sim.mode() != SimMode::WindTunnel) return out; // no gas grid, no fields to probe
    if (probe.points <= 1) measurePoint(sim, probe, out);
    else if (!isMagnetic(probe.quantity) || sim.grid.magnetic.enabled) measureLine(sim, probe, out);
    return out;
}

} // namespace rf
