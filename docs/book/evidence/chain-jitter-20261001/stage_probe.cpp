// Shadow-replay one frame through the compiled SDK stages and check every substep's mechanical
// state against the actual atomic step. Restore the reference checkpoint after each shadow.
#define main frameProbeMain
#include "jitter_probe.cpp"
#undef main
#include "rigid/RigidStepCheckpoint.h"

uint64_t motionHash(const RigidWorld& w) {
    StateHash h;for(const auto& b:w.bodies()){h.add(b.pos);h.add(b.rot);h.add(b.vel);h.add(b.angVel);}return h.h;
}

void row(std::ofstream& csv,const RigidWorld& w,int step,const char* stage,int iteration=0) {
    static Energy previous;
    const auto e=energy(w);int locks=0,corrections=0;double splitU=0;
    if(std::string(stage)=="start")previous=e;
    for(const auto& b:w.bodies())if(b.invMass>0)splitU-=double(b.mass)*dot(Double3(w.params.gravity),Double3(b.biasVel))*double(w.lastDt_);
    for(const auto& m:w.manifolds_){locks+=m.locked;for(const auto& p:m.points)corrections+=p.positionBias>0;}
    csv << step << ',' << stage << ',' << iteration << ',' << numberText(e.kinetic) << ',' << numberText(e.potential) << ','
        << numberText(e.total()) << ',' << numberText(splitU) << ',' << numberText(w.deepestPenetration()) << ','
        << w.contactCount() << ',' << w.warmStartStats().matched << ',' << locks << ',' << corrections << ','
        << numberText(e.kinetic-previous.kinetic) << ',' << numberText(e.potential-previous.potential) << ','
        << numberText(e.total()-previous.total()) << '\n';
    previous=e;
}

bool shadow(RigidWorld& w,float h,int step,std::ofstream& csv) {
    RigidStepCheckpoint before(w);w.step(h);RigidStepCheckpoint after(w);const auto expected=motionHash(w);
    before.restore(w);w.lastDt_=w.correctionDt_=h;w.ccdHits_=0;w.ccdDiagnostics_=CcdDiagnostics();
    row(csv,w,step,"start");w.beginStep();w.integrateVelocities(h);row(csv,w,step,"gravity-kick");
    w.collide();row(csv,w,step,"collide");w.prepare(h);w.prepareGrab(h);row(csv,w,step,"warm-start");
    for(int i=0;i<w.params.iterations;++i){w.solve();w.solveGrab(h);row(csv,w,step,"contact-iteration",i+1);}
    w.applyRestitution();row(csv,w,step,"restitution");
    if(w.params.shockPropagation && !w.manifolds_.empty())w.propagateShock();
    row(csv,w,step,"shock");
    w.rememberContactImpulses();w.dampRestingBodies(h);row(csv,w,step,"rest-damping");
    RigidStepCheckpoint driftStart(w);
    for(auto& b:w.bodies()){b.biasVel=Vector3(0);b.biasAngVel=Vector3(0);}
    w.integratePoses(h);const auto noBias=energy(w);driftStart.restore(w);
    w.integratePoses(h);const auto withBias=energy(w);row(csv,w,step,"pose-drift");
    csv << step << ",split-counterfactual,0," << numberText(withBias.kinetic) << ',' << numberText(withBias.potential)
        << ',' << numberText(withBias.total()) << ",0,0,0,0,0,0," << numberText(withBias.kinetic-noBias.kinetic)
        << ',' << numberText(withBias.potential-noBias.potential) << ',' << numberText(withBias.total()-noBias.total()) << '\n';
    w.finishStep(h,std::chrono::steady_clock::now());row(csv,w,step,"finish");
    const bool same=motionHash(w)==expected && w.ccdHits_==0;after.restore(w);return same;
}

int main(int argc,char** argv) {
    try {
        if(argc!=4)throw std::runtime_error("FRAME PULL_0_OR_1 CSV");
        const int frame=std::stoi(argv[1]),pull=std::stoi(argv[2]);
        if(frame<601 || frame>7200 || pull<0 || pull>1)throw std::runtime_error("bounds");
        Simulation sim;loadSample(sim,Preset::HangingTorus);Vector3 start;
        for(int f=0;f<frame-1;++f){if(pull)mouse(sim,f,start);sim.stepFrame();}
        auto& w=sim.rigid;const auto initial=energy(w);const float h=sim.frameDt/float(w.params.substeps);
        if(w.grabJoint().active || !w.joints().empty() || w.params.sleeping)throw std::runtime_error("unsupported shadow");
        std::ofstream csv(argv[3]);if(!csv)throw std::runtime_error("CSV");
        csv << "substep,stage,iteration,kinetic_J,potential_J,total_J,predicted_split_U_J,depth_m,contacts,matched,locked,corrected_points,delta_T_J,delta_U_J,delta_E_J\n";
        for(int step=0;step<w.params.substeps;++step)if(!shadow(w,h,step,csv))throw std::runtime_error("shadow does not match atomic reference");
        const auto final=energy(w);
        std::cout << "{\"frame\":" << frame << ",\"pull\":" << pull << ",\"matched_substeps\":" << w.params.substeps
            << ",\"delta_E_J\":" << numberText(final.total()-initial.total()) << ",\"delta_T_J\":" << numberText(final.kinetic-initial.kinetic)
            << ",\"delta_U_J\":" << numberText(final.potential-initial.potential) << "}\n";
        return 0;
    }catch(const std::exception& e){std::cerr << e.what() << '\n';return 2;}
}
