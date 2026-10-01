// Free nut on a faceted square thread: initial clearance and contact-produced screw travel.
// The ideal pitch relation is an independent kinematic oracle, with an explicit backlash band.
#include "TestRunner.h"
#include "Tests.h"
#include "samples/NutBoltScene.h"
#include "core/Format.h"

#include <filesystem>
#include <fstream>

namespace {

std::ofstream nutCsv(const char* name, const char* header) {
    const char* dir = std::getenv("RF_PLOT_DIR");
    if (!dir || !*dir) return {};
    std::filesystem::create_directories(dir);
    std::ofstream out(std::filesystem::path(dir) / name);
    CHECK(bool(out), "cannot write %s", name);
    out << header << '\n'; return out;
}

void nutPoses(std::ofstream& csv, int frame, const Simulation& sim) {
    if (!csv) return;
    for (size_t i=0; i<sim.rigid.bodies().size(); ++i) {
        const auto& b=sim.rigid.bodies()[i];
        csv << frame << ',' << i;
        for (float x : {b.pos.x,b.pos.y,b.pos.z,b.rot.w,b.rot.x,b.rot.y,b.rot.z,
                        b.vel.x,b.vel.y,b.vel.z,b.angVel.x,b.angVel.y,b.angVel.z}) csv << ',' << numberText(x);
        csv << '\n';
    }
}

double nutEnergy(const Simulation& sim) {
    const auto& b=sim.rigid.bodies()[NutBoltScene::nutBody];
    return sim.rigid.kineticEnergy()-double(b.mass)*dot(sim.rigid.params.gravity,b.pos);
}

} // namespace

void testNutBoltAssembly() {
    Simulation sim; loadSample(sim,Preset::NutBolt);
    CHECK(sim.rigid.bodies().size()==2 && sim.rigid.joints().empty(), "contact-only nut, no screw joint");
    CHECK(sim.rigid.bodies()[0].invMass==0 && sim.rigid.bodies()[1].invMass>0, "fixed bolt, dynamic nut");
    auto& b=sim.rigid.bodies()[1];
    const auto& shape=static_cast<const CompoundShape&>(*b.shape);
    const Vector3 origin=NutBoltScene::modelOrigin(b);
    for (float angle : {-0.4f,0.0f,0.4f}) {
        const Quaternion rotation=Quaternion::fromAxisAngle({0,1,0},angle);
        b.rot=(rotation*Quaternion::fromMatrix3x3(shape.principalRotation())).normalized();
        b.pos=origin+Vector3(0,NutBoltScene::pitch*angle/(2*kPi),0)+rotation.rotate(shape.centerOfMass());
        b.updateInertia();
        CHECK(maxPartOverlap(sim.rigid)<2e-6f, "ideal screw clearance lost at %.2f rad",angle);
    }
    CHECK(shape.children().size()>24, "nut needs an internal thread, not just a hole");
    const auto& bolt=static_cast<const CompoundShape&>(*sim.rigid.bodies()[0].shape);
    std::printf("  bolt/nut convex pieces %zu/%zu, pitch %.1f mm\n",bolt.children().size(),shape.children().size(),1000*NutBoltScene::pitch);
}

void testNutBoltCycle() {
    Simulation sim; loadSample(sim,Preset::NutBolt);
    const auto* scene=static_cast<const NutBoltScene*>(sim.scene());
    auto csv=nutCsv("nut-bolt.csv","frame,time_s,height_m,angle_rad,speed_rad_s,energy_J,drive_work_J,residual_J,screw_error_m,sampled_depth_m");
    auto poses=nutCsv("nut-bolt-poses.csv","frame,body,x,y,z,qw,qx,qy,qz,vx,vy,vz,wx,wy,wz");
    const double initial=nutEnergy(sim);
    const Vector3 fixed=sim.rigid.bodies()[0].pos;
    float baseY=0,baseAngle=0,lift=0,returnError=0,worstPitch=0,depth=0;
    double peakResidual=0;
    nutPoses(poses,0,sim);
    for (int f=1; f<=720; ++f) {
        sim.stepFrame();
        const auto& b=sim.rigid.bodies()[1];
        const float y=NutBoltScene::modelOrigin(b).y, angle=NutBoltScene::modelAngle(b);
        CHECK(std::isfinite(length2(b.pos)+length2(b.vel)+length2(b.angVel)), "non-finite nut at %d",f);
        CHECK(length2(sim.rigid.bodies()[0].pos-fixed)==0 && b.worldBounds().lo.y>0.08f, "support/floor changed at %d",f);
        if (f==180) { baseY=y; baseAngle=angle; }
        const float error=f>=180 ? y-baseY-NutBoltScene::pitch*(angle-baseAngle)/(2*kPi) : 0;
        if (f>=180) worstPitch=std::max(worstPitch,std::fabs(error));
        if (f==360) lift=y-baseY;
        if (f==540) returnError=std::fabs(y-baseY);
        if (f%60==0) { depth=std::max(depth,maxPartOverlap(sim.rigid)); nutPoses(poses,f,sim); }
        const double e=nutEnergy(sim), residual=e-initial-scene->driveWork();
        peakResidual=std::max(peakResidual,residual);
        if (csv) csv << f << ',' << numberText(double(f)/60) << ',' << numberText(y) << ',' << numberText(angle) << ','
            << numberText(b.angVel.y) << ',' << numberText(e) << ',' << numberText(scene->driveWork()) << ','
            << numberText(residual) << ',' << numberText(error) << ',' << numberText(depth) << '\n';
    }
    std::printf("  nut: lift %.3f mm, return %.3f mm, pitch error %.3f mm, sampled depth %.3f mm, unexplained gain %.6g J\n",
                1000*lift,1000*returnError,1000*worstPitch,1000*depth,peakResidual);
    CHECK(lift>0.003f && returnError<0.0015f && worstPitch<0.0015f, "contact did not follow pitch within backlash");
    CHECK(depth<0.0005f && peakResidual<0.005, "deep penetration or unexplained energy gain");
}
