// Bounded actuator-energy study of preset 42. Usage: FRAMES LOAD SUBSTEPS MEASURE CSV.
// CSV is observational: E-E0-Wmotor includes dissipation and all numerical energy defects.
#include "samples/ChainWheelScene.h"
#include "samples/Samples.h"
#include "core/Format.h"

#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace {

double energy(const rf::Simulation& sim) {
    double e = sim.rigid.kineticEnergy();
    for (const auto& b : sim.rigid.bodies()) if (b.invMass > 0)
        e -= double(b.mass) * dot(sim.rigid.params.gravity, b.pos);
    return e;
}

void hashFrame(uint64_t& hash, const rf::RigidWorld& w) {
    for (const auto& b : w.bodies())
        for (float v : {b.pos.x,b.pos.y,b.pos.z,b.rot.w,b.rot.x,b.rot.y,b.rot.z,
                        b.vel.x,b.vel.y,b.vel.z,b.angVel.x,b.angVel.y,b.angVel.z}) {
            uint32_t bits; std::memcpy(&bits,&v,sizeof(bits)); hash=(hash^bits)*1099511628211ull;
        }
}

void study(int frames, float load, int substeps, bool measure, const char* path) {
    rf::Simulation sim; rf::loadSample(sim, rf::Preset::ChainDrive);
    sim.scene()->setParam(0,load); sim.rigid.params.substeps=substeps;
    if (std::getenv("RF_CHAIN_UNLOCK")) sim.rigid.params.rotationalLock=false;
    if (std::getenv("RF_CHAIN_NO_SHOCK")) sim.rigid.params.shockPropagation=false;
    sim.rigid.params.measureMotorWork=measure; sim.reset();
    const double initial=energy(sim);
    double peakResidual=0, peakFrameResidual=0, previous=initial, settledResidual=0, peakAfterSettle=0;
    uint64_t hash=14695981039346656037ull;
    std::ofstream csv(path); if (!csv) throw std::runtime_error("cannot write CSV");
    csv << "frame,time_s,energy_J,motor_work_J,residual_J,angle_rad,speed_rad_s,tip_y_m,broken_pairs\n";
    for (int f=0;f<frames;++f) {
        const double beforeWork=sim.rigid.motorWork(); sim.stepFrame(); hashFrame(hash,sim.rigid);
        const double e=energy(sim), residual=e-initial-sim.rigid.motorWork();
        peakResidual=std::max(peakResidual,residual);
        if (f==119) settledResidual=residual;
        if (f>=120) peakAfterSettle=std::max(peakAfterSettle,residual-settledResidual);
        peakFrameResidual=std::max(peakFrameResidual,e-previous-(sim.rigid.motorWork()-beforeWork));
        previous=e;
        const auto* scene=static_cast<const rf::ChainWheelScene*>(sim.scene());
        if (!std::isfinite(e) || scene->brokenPairs()) throw std::runtime_error("invalid state or broken chain");
        const auto& motor=static_cast<const rf::HingeJoint&>(*sim.rigid.joints().front());
        csv << f+1 << ',' << rf::numberText(double(f+1)/60) << ',' << rf::numberText(e) << ','
            << rf::numberText(sim.rigid.motorWork()) << ',' << rf::numberText(residual) << ','
            << rf::numberText(motor.angle(sim.rigid.bodies())) << ',' << rf::numberText(sim.rigid.bodies()[0].angVel.z)
            << ',' << rf::numberText(sim.rigid.bodies().back().pos.y) << ',' << scene->brokenPairs() << '\n';
    }
    std::cout << "{\"frames\":" << frames << ",\"load\":" << rf::numberText(load) << ",\"substeps\":" << substeps
        << ",\"measured\":" << (measure ? "true" : "false") << ",\"motor_work_J\":" << rf::numberText(sim.rigid.motorWork())
        << ",\"peak_residual_J\":" << rf::numberText(peakResidual) << ",\"peak_frame_residual_J\":" << rf::numberText(peakFrameResidual)
        << ",\"peak_after_settle_J\":" << rf::numberText(peakAfterSettle)
        << ",\"final_residual_J\":" << rf::numberText(energy(sim)-initial-sim.rigid.motorWork())
        << ",\"trajectory_hash\":\"" << std::hex << hash << std::dec << "\",\"broken_pairs\":0}\n";
    if (measure && (peakResidual>0.02 || peakAfterSettle>0.02))
        throw std::runtime_error("unexplained energy gain above 0.02 J (initial or settled baseline)");
}

} // namespace

int main(int argc,char** argv) {
    try {
        float f=0,l=0,s=0,m=0;
        if (argc!=6 || !rf::parseNumber(argv[1],f) || !rf::parseNumber(argv[2],l) || !rf::parseNumber(argv[3],s)
            || !rf::parseNumber(argv[4],m) || !std::isfinite(f+l+s+m) || f<120 || f>3600 || f!=std::floor(f)
            || l<1 || l>3 || s<1 || s>160 || s!=std::floor(s) || (m!=0 && m!=1))
            throw std::runtime_error("FRAMES 120..3600 LOAD 1..3 SUBSTEPS 1..160 MEASURE 0/1 CSV");
        study(int(f),l,int(s),m!=0,argv[5]); return 0;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 2; }
}
