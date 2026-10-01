// Trace the half-step mixed CCD failure without changing SDK state or solver settings.
#include "scene/Simulation.h"
#include "samples/Samples.h"
#include <cstdio>
#include <cstdlib>
using namespace rf;
static double kinetic(const RigidBody& b) {
    const Vector3 u = b.rotation().transposed() * b.angVel;
    return 0.5 * b.mass * length2(b.vel) + 0.5 * (u.x*u.x/b.invInertiaLocal.x + u.y*u.y/b.invInertiaLocal.y + u.z*u.z/b.invInertiaLocal.z);
}
int main(int argc, char**) {
    Simulation sim; loadSample(sim, Preset::RigidCcd);
    auto& w = sim.rigid;
    if (argc > 1) w.params.ccd = false;
    for (int k=0; k<60; ++k) {
        auto before = w.bodies();
        w.step(1.0f/1200);
        int hottest=-1; double heat=0;
        for (int i=0;i<int(before.size());++i) {
            const auto& b=w.bodies()[i];
            if (b.invMass <= 0) continue;
            double e=kinetic(b);
            if (e>heat || !std::isfinite(e)) { hottest=i; heat=e; }
        }
        const auto& b=w.bodies()[hottest]; const auto& old=before[hottest];
        std::printf("step %d id %d K %.9g -> %.9g p %.9g %.9g %.9g v %.9g %.9g %.9g w %.9g %.9g %.9g hits %zu unresolved %zu contacts %zu depth %.9g\n", k+1,hottest,kinetic(old),heat,b.pos.x,b.pos.y,b.pos.z,b.vel.x,b.vel.y,b.vel.z,b.angVel.x,b.angVel.y,b.angVel.z,w.ccdHits(),w.ccdDiagnostics().unresolved,w.contactCount(),w.deepestPenetration());
        if (!std::isfinite(heat)) break;
    }
}
