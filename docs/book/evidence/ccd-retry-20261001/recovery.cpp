// Three spheres isolate a CCD retry whose split correction did not shrink with dt.
#include "rigid/RigidWorld.h"
#include <cstdio>
using namespace rf;
int main() {
    RigidWorld w;
    w.params.gravity = Vector3(0);
    w.params.sleeping = w.params.collideWithDomain = false;
    w.params.linearDamping = w.params.angularDamping = 0;
    w.params.ccdAllMoving = true;
    w.addSphere({0,0,0},0.1f,0,Vector3(1));
    const int moving=w.addSphere({0.1f,0,0},0.1f,1000,Vector3(1));
    w.addSphere({0.32f,0,0},0.1f,0,Vector3(1));
    const bool ok=w.tryStep(0.01f);
    const auto& r=w.lastStepResult();
    std::printf("completed=%d status=%d accepted=%zu rejected=%zu x=%.9g energy=%.9g\n",ok,int(r.status),r.acceptedSubsteps,r.rejectedTrials,w.bodies()[moving].pos.x,w.kineticEnergy());
    return ok ? 0 : 1;
}
