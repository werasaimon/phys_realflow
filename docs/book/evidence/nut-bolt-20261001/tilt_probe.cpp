// Static assembly admissibility under small transverse tilts, before any solver step.
// Rotate around the nut model origin; this is not a dynamic threading experiment.
#include "tests/TestRunner.h"
#include "samples/NutBoltScene.h"
#include "core/Format.h"
#include <fstream>

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    std::ofstream summary(argv[1]), poses(argv[2]);
    if (!summary || !poses) return 2;
    summary << "tilt_deg,initial_gjk_depth_m\n";
    poses << "frame,body,x,y,z,qw,qx,qy,qz,vx,vy,vz,wx,wy,wz\n";
    int frame = 0;
    for (float degrees : {0.0f, 0.25f, 0.5f, 1.0f, 2.0f, 5.0f}) {
        Simulation sim; loadSample(sim, Preset::NutBolt);
        auto& nut = sim.rigid.bodies()[NutBoltScene::nutBody];
        const auto& shape = static_cast<const CompoundShape&>(*nut.shape);
        const Vector3 origin = NutBoltScene::modelOrigin(nut);
        const Quaternion tilt = Quaternion::fromAxisAngle({1, 0, 0}, degrees * kPi / 180);
        nut.rot = (tilt * Quaternion::fromMatrix3x3(shape.principalRotation())).normalized();
        nut.pos = origin + tilt.rotate(shape.centerOfMass());
        nut.updateInertia();
        summary << numberText(degrees) << ',' << numberText(maxPartOverlap(sim.rigid)) << '\n';
        for (size_t i = 0; i < sim.rigid.bodies().size(); ++i) {
            const auto& b = sim.rigid.bodies()[i];
            poses << frame << ',' << i;
            for (float x : {b.pos.x,b.pos.y,b.pos.z,b.rot.w,b.rot.x,b.rot.y,b.rot.z,
                            b.vel.x,b.vel.y,b.vel.z,b.angVel.x,b.angVel.y,b.angVel.z})
                poses << ',' << numberText(x);
            poses << '\n';
        }
        ++frame;
    }
    return summary && poses ? 0 : 2;
}
