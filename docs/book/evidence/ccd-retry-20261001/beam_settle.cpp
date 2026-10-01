// Bounded numerical profile of the beam and the mixed CCD sample at three time steps.
// Compile with C++17, -O3 -march=native, -Isrc -I., librfsamples.a, librfcore.a and -lpthread.
#include "scene/Simulation.h"
#include "samples/Samples.h"
#include "core/Format.h"
#include <chrono>
#include <iostream>
#include <cstdlib>
using namespace rf;

static double energy(const RigidWorld& world) {
    double total = world.kineticEnergy();
    for (const RigidBody& body : world.bodies())
        if (body.mass > 0) total -= body.mass * dot(world.params.gravity, body.pos);
    return total;
}

static void beam(RigidWorld& world) {
    world.clear();
    world.setStaticMesh(nullptr);
    world.setDomain(AABB({-5, 0, -5}, {5, 10, 5}));
    world.params.sleeping = false;
    for (int i = 0; i < 6; ++i)
        world.addBox({-0.75f + 0.3f * i, 0.1f, 0}, Vector3(0.1f), Quaternion(), 500, Vector3(1));
    const int id = world.addBox({0, 1.5f, 0}, {1, 0.03f, 0.03f}, Quaternion::fromAxisAngle({0, 1, 0}, 0.3f), 800, Vector3(1));
    world.bodies()[id].vel = {0, -25, 0};
    world.bodies()[id].angVel = {0, 40, 0};
}

int main() {
 for(bool all : {false,true}) { RigidWorld w;beam(w);w.params.ccdAllMoving=all;
 for(int i=0;i<3000;++i) {w.step(1.f/600);if((i+1)%300==0) {auto&b=w.bodies().back();std::printf("all=%d t=%.1f x=%.6g y=%.6g z=%.6g v=%.8g w=%.8g E=%.8g depth=%.7g\n",int(all),(i+1)/600.,b.pos.x,b.pos.y,b.pos.z,length(b.vel),length(b.angVel),energy(w),w.deepestPenetration());}}
 }
}
