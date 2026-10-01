// Observe the real short-chain preset without changing the solver or enabling sleep.
// GCC -fno-access-control exposes contact metadata only for this external diagnostic.
#include "tests/TestRunner.h"
#include "samples/HangingTorusScene.h"
#include "samples/TorusChainsScene.h"
#include "scene/ChannelMath.h"
#include "core/Format.h"
#include <fstream>
#include <iostream>
#include <stdexcept>

struct Energy { double kinetic=0, potential=0; double total() const { return kinetic+potential; } };

Energy energy(const RigidWorld& w) {
    Energy e;
    for (const auto& b : w.bodies()) if (b.alive && b.invMass>0) {
        const Double3 v(b.vel), p(b.pos), g(w.params.gravity);
        const auto R=b.rotation(); const auto a=b.angVel;
        const double x=double(R.m[0][0])*a.x+double(R.m[1][0])*a.y+double(R.m[2][0])*a.z;
        const double y=double(R.m[0][1])*a.x+double(R.m[1][1])*a.y+double(R.m[2][1])*a.z;
        const double z=double(R.m[0][2])*a.x+double(R.m[1][2])*a.y+double(R.m[2][2])*a.z;
        e.kinetic+=0.5*double(b.mass)*dot(v,v)+0.5*(x*x/b.invInertiaLocal.x+y*y/b.invInertiaLocal.y+z*z/b.invInertiaLocal.z);
        e.potential-=double(b.mass)*dot(g,p);
    }
    return e;
}

Vector3 startGrab(Simulation& sim) {
    const auto& b=sim.rigid.bodies().back();
    const auto& shape=static_cast<const CompoundShape&>(*b.shape);
    const Vector3 rim=b.pos+b.rotation()*shape.principalRotation().transposed()*(Vector3(0,0,TorusChainsScene::majorRadius)-shape.centerOfMass());
    const auto d=normalize(Vector3(1,0,1)); const auto origin=rim+2*d;
    int hit=-1;float distance=0;Vector3 normal;
    if (!sim.rigid.raycast(origin,-d,3,hit,distance,normal) || hit!=7) throw std::runtime_error("rim grab failed");
    const auto p=origin-d*distance;sim.rigid.grab(hit,p);return p;
}

void mouse(Simulation& sim, int frame, Vector3& start) {
    if (frame==60) start=startGrab(sim);
    if (frame>=60 && frame<240) {
        const float phase=2*kPi*float(frame-60)/72;
        const auto d=frame>=150 && frame<180 ? Vector3(-0.22f,0.18f,-0.08f)
            : Vector3(0.18f*rf::sin(phase),0.12f*(1-rf::cos(phase)),0.08f*rf::sin(0.5f*phase));
        sim.rigid.setGrabTarget(start+d);
    }
    if (frame==240) sim.rigid.releaseGrab();
}

int main(int argc,char** argv) {
    try {
        if (argc!=5) throw std::runtime_error("FRAMES PULL_0_OR_1 VARIANT_0_TO_4 CSV");
        const int frames=std::stoi(argv[1]),pull=std::stoi(argv[2]),variant=std::stoi(argv[3]);
        if (frames<600 || frames>7200 || pull<0 || pull>1 || variant<0 || variant>4) throw std::runtime_error("bounds");
        Simulation sim;loadSample(sim,Preset::HangingTorus);auto& w=sim.rigid;
        if (variant==1) w.params.shockPropagation=false;
        if (variant==2) w.params.rotationalLock=false;
        if (variant==3) w.params.warmStarting=false;
        if (variant==4) w.params.splitImpulse=false;
        if (w.params.sleeping || !w.joints().empty()) throw std::runtime_error("unexpected active scene");
        std::ofstream csv(argv[4]);if(!csv) throw std::runtime_error("CSV");
        csv << "frame,time_s,kinetic_J,potential_J,total_J,delta_E_J,ui_total_J,depth_m,contacts,matched,locked,max_speed_m_s,max_spin_rad_s,max_move_m,tip_y_m,broken\n";
        auto previous=energy(w);const double initial=previous.total();Vector3 start;
        double peakE=initial,lateGain=0,lateKinetic=0,peakLateKinetic=0,peakLateMove=0,depth=0;
        int worstFrame=0,peakFrame=0,broken=0;StateHash trajectory;
        for (int f=0;f<frames;++f) {
            if(pull)mouse(sim,f,start);
            std::vector<Vector3> before;for(const auto& b:w.bodies())before.push_back(b.pos);
            sim.stepFrame();const auto e=energy(w);if(!std::isfinite(e.total()))throw std::runtime_error("nonfinite energy");
            double v=0,spin=0,move=0;int locked=0;
            for(size_t i=0;i<w.bodies().size();++i) { const auto& b=w.bodies()[i];
                v=std::max(v,double(length(b.vel)));spin=std::max(spin,double(length(b.angVel)));
                move=std::max(move,double(length(b.pos-before[i])));trajectory.add(b.pos);trajectory.add(b.rot);trajectory.add(b.vel);trajectory.add(b.angVel);
            }
            for(const auto& m:w.manifolds_)locked+=m.locked;
            double ui=0;for(const auto& m:measureScene(sim))if(m.info.id=="rigid/mechanical energy")ui=m.value;
            const double gain=e.total()-previous.total();peakE=std::max(peakE,e.total());depth=std::max(depth,double(w.deepestPenetration()));
            if(f>=600) { if(gain>lateGain){lateGain=gain;worstFrame=f+1;}if(e.kinetic>peakLateKinetic){peakLateKinetic=e.kinetic;peakFrame=f+1;}
                peakLateMove=std::max(peakLateMove,move);if(f>=frames-60)lateKinetic+=e.kinetic/60; }
            broken=std::max(broken,static_cast<const HangingTorusScene*>(sim.scene())->brokenPairs());
            csv << f+1 << ',' << numberText(double(f+1)/60) << ',' << numberText(e.kinetic) << ',' << numberText(e.potential) << ','
                << numberText(e.total()) << ',' << numberText(gain) << ',' << numberText(ui) << ',' << numberText(w.deepestPenetration()) << ','
                << w.contactCount() << ',' << w.warmStartStats().matched << ',' << locked << ',' << numberText(v) << ',' << numberText(spin) << ','
                << numberText(move) << ',' << numberText(w.bodies().back().pos.y) << ',' << broken << '\n';
            previous=e;
        }
        std::cout << "{\"frames\":" << frames << ",\"pull\":" << pull << ",\"variant\":" << variant << ",\"initial_E_J\":" << numberText(initial)
            << ",\"peak_E_J\":" << numberText(peakE) << ",\"late_frame_gain_J\":" << numberText(lateGain) << ",\"gain_frame\":" << worstFrame
            << ",\"late_peak_T_J\":" << numberText(peakLateKinetic) << ",\"T_frame\":" << peakFrame << ",\"last_second_T_J\":" << numberText(lateKinetic)
            << ",\"late_max_frame_move_m\":" << numberText(peakLateMove) << ",\"max_contact_depth_m\":" << numberText(depth) << ",\"broken\":" << broken
            << ",\"trajectory_hash\":\"" << std::hex << trajectory.h << std::dec << "\"}\n";
        return broken ? 1 : 0;
    }catch(const std::exception& e){std::cerr << e.what() << '\n';return 2;}
}
