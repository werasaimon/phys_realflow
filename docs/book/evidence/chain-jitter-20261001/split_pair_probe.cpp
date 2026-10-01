// Analytical dynamic-pair check of zero-target split contact rows. Real velocities
// must remain unchanged; the paired pseudo impulse must conserve linear momentum.
#include "rigid/RigidWorld.h"
#include "core/Format.h"
#include <iostream>
using namespace rf;
int main() {
    RigidWorld w;
    for(int i=0;i<2;++i){RigidBody b;b.mass=1;b.invMass=1;b.pos={float(i),0,0};
        b.invInertiaLocal={1,1,1};b.updateInertia();w.bodies().push_back(b);}
    auto& a=w.bodies()[0];auto& b=w.bodies()[1];a.biasVel={0,-1,0};
    RigidWorld::Manifold m;m.a=0;m.b=1;
    RigidWorld::SolverPoint p;p.position={0.5f,0,0};p.normal={0,1,0};
    p.depth=w.params.slop/2;p.positionBias=0;p.massN=1/2.5f;m.points.push_back(p);
    const auto beforeMomentum=a.mass*a.biasVel+b.mass*b.biasVel;
    const auto beforeL=cross(a.pos,a.mass*a.biasVel)+cross(b.pos,b.mass*b.biasVel)
        +a.biasAngVel+b.biasAngVel;
    w.solveSplitImpulse(m);
    const float afterSpeed=dot(a.biasVel+cross(a.biasAngVel,p.position-a.pos)
        -b.biasVel-cross(b.biasAngVel,p.position-b.pos),p.normal);
    const auto momentum=a.mass*a.biasVel+b.mass*b.biasVel;
    const auto afterL=cross(a.pos,a.mass*a.biasVel)+cross(b.pos,b.mass*b.biasVel)
        +a.biasAngVel+b.biasAngVel;
    const float realChange=length(a.vel)+length(b.vel)+length(a.angVel)+length(b.angVel);
    std::cout << "{\"relative_split_speed_before\":-1,\"relative_split_speed_after\":"
        << numberText(afterSpeed) << ",\"split_impulse\":" << numberText(m.points[0].jp)
        << ",\"expected_impulse\":0.4,\"momentum_error\":" << numberText(length(momentum-beforeMomentum))
        << ",\"angular_momentum_error\":" << numberText(length(afterL-beforeL))
        << ",\"real_velocity_change\":" << numberText(realChange) << "}\n";
    return !std::isfinite(afterSpeed);
}
