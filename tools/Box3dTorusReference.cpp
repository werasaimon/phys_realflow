// Optional Box3D backend runner for SDK sample 39; the SDK core has no Box3D dependency.
// Usage: rf_box3d_reference NEW_DIRECTORY [substeps=10] [frames=600] [full_steps=0] [length_scale=1]
// full_steps=1 refreshes native collisions at every 1/(60*substeps) s. Scale 10 means native slop 0.5 mm in SI.
#include "samples/TorusChainsScene.h"
#include "samples/Samples.h"
#include "core/Format.h"
#include "rigid_reference/Box3dRigidBackend.h"
#include "rigid_reference/ReferenceTiming.h"
#include <filesystem>
#include <fstream>
#include <iostream>

namespace {
void pose(std::ostream& out,const rf::RigidBody& b) {
    const float fields[]={b.pos.x,b.pos.y,b.pos.z,b.rot.w,b.rot.x,b.rot.y,b.rot.z,
        b.vel.x,b.vel.y,b.vel.z,b.angVel.x,b.angVel.y,b.angVel.z};
    for (float v:fields) out << ',' << rf::numberText(v);
}
void checkImport(const std::vector<rf::RigidBody>& initial,const rf::RigidWorld& world) {
    for (size_t i=0;i<initial.size();++i) {
        const auto& a=initial[i]; const auto& b=world.bodies()[i];
        if (rf::length(a.pos-b.pos)>1e-6f || rf::length(a.vel-b.vel)>1e-6f ||
            rf::length(a.rot.rotate({1,2,3})-b.rot.rotate({1,2,3}))>1e-5f)
            throw std::runtime_error("Box3D pose/velocity import roundtrip failed");
    }
}
}

int main(int argc,char** argv) {
    try {
        if (argc<2 || argc>6) throw std::runtime_error("NEW_DIRECTORY [substeps] [frames] [full_steps] [length_scale]");
        const int steps=argc>2 ? std::stoi(argv[2]) : 10;
        const int frames=argc>3 ? std::stoi(argv[3]) : 600;
        const bool fullSteps=argc>4 && std::stoi(argv[4])!=0;
        float scale=1;
        if (argc>5 && !rf::parseNumber(argv[5],scale)) throw std::runtime_error("invalid length scale");
        if (steps<1 || frames<1) throw std::runtime_error("positive counts required");
        const std::filesystem::path output(argv[1]);
        if (!std::filesystem::create_directory(output)) throw std::runtime_error("output must be a new directory");
        rf::Simulation sim; rf::loadSample(sim,rf::Preset::TorusChains);
        auto* scene=static_cast<rf::TorusChainsScene*>(sim.scene());
        const auto initial=sim.rigid.bodies();
        rf::reference::Box3dRigidBackend reference(scale); reference.load(sim.rigid);
        checkImport(initial,sim.rigid); scene->afterStep(sim);
        if (scene->brokenPairs()) throw std::runtime_error("Box3D import lost linking");
        std::ofstream csv(output/"frames.csv"), summary(output/"summary.txt"), poses(output/"pair-poses.csv"), all(output/"body-poses.csv");
        if (!csv || !summary || !poses || !all) throw std::runtime_error("cannot write report");
        csv << "frame,time_s,broken_pairs,deepest_box3d_contact_m,engine_ms,observer_ms,pairs_ms,collide_ms,solve_ms,bullets_ms\n";
        poses << "substep,ax,ay,az,aqw,aqx,aqy,aqz,avx,avy,avz,awx,awy,awz,bx,by,bz,bqw,bqx,bqy,bqz,bvx,bvy,bvz,bwx,bwy,bwz\n";
        all << "frame,body,x,y,z,qw,qx,qy,qz,vx,vy,vz,wx,wy,wz\n";
        int first=-1,worst=0;
        for (int frame=1;frame<=frames;++frame) {
            ReferenceTiming timing;
            double pairsMs=0,collideMs=0,solveMs=0,bulletsMs=0;
            const int calls=fullSteps ? steps : 1;
            for (int sub=0;sub<calls;++sub) {
                reference.step(sim.rigid,(1.0f/60)/float(calls),fullSteps ? 1 : steps);
                const auto elapsed=reference.timings(); const auto p=reference.profile();
                timing.engineMs+=elapsed.engine; timing.observerMs+=elapsed.input+elapsed.output;
                pairsMs+=p.pairs; collideMs+=p.collide; solveMs+=p.solve; bulletsMs+=p.bullets;
                timing.start(); scene->afterStep(sim);
                worst=std::max(worst,scene->brokenPairs());
                if (scene->brokenPairs() && first<0) first=frame;
                poses << ((frame-1)*steps+(fullSteps ? sub+1 : steps));
                pose(poses,sim.rigid.bodies()[67]); pose(poses,sim.rigid.bodies()[68]); poses << '\n';
                timing.finishObserver();
            }
            csv << frame << ',' << rf::numberText(double(frame)/60) << ',' << scene->brokenPairs()
                << ',' << rf::numberText(reference.deepestContact()) << timing.csv() << ',' << rf::numberText(pairsMs)
                << ',' << rf::numberText(collideMs) << ',' << rf::numberText(solveMs) << ',' << rf::numberText(bulletsMs) << '\n';
            for (size_t i=0;i<sim.rigid.bodies().size();++i) { all << frame << ',' << i; pose(all,sim.rigid.bodies()[i]); all << '\n'; }
            csv.flush(); poses.flush(); all.flush();
            if (frame%30==0) std::cout << "frame " << frame << ", worst broken " << worst << std::endl;
        }
        summary << "Box3D commit 9f998c862d54c03a633ecea3831937385c78b532\nframes " << frames << "\nsubsteps " << steps
            << "\nfull_collision_steps " << fullSteps << "\nnumerical_length_scale " << rf::numberText(scale)
            << "\nnominal_native_slop_in_SI_m " << rf::numberText(.005f/scale)
            << "\nfirst_break_frame " << first << "\nworst_broken " << worst
            << "\nmax_reported_contact_depth_m " << rf::numberText(reference.deepestContact())
            << "\nhull_face_support_error_m " << rf::numberText(reference.supportImportError())
            << "\nimported_convex_parts " << reference.importedPartCount()
            << "\nexact union partition for native hull limits; SDK mass/inertia and initial poses preserved\n"
            << "single thread; native solver defaults; numerical scaling preserves free dynamics and scales native geometric tolerances\n"
            << "native manifold depth is not a swept or global geometry certificate\n";
        std::cout << "first break " << first << ", worst broken " << worst << ", contact depth " << rf::numberText(reference.deepestContact()) << " m\n";
        return worst ? 1 : 0;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 2; }
}
