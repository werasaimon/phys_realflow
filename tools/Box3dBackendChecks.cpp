// Analytical SI tests for the optional native Box3D backend, including numerical length rescaling.
// Cold torus fixture is the actual SDK contact-starvation capture; it is already overlapped at release.
#include "rigid_reference/Box3dRigidBackend.h"
#include "rigid_reference/SatCcd.h"
#include "samples/Models.h"
#include "core/Format.h"
#include <iostream>
#include <stdexcept>

namespace {
using namespace rf;
using rf::reference::Box3dRigidBackend;
int checks=0;
void require(bool ok,const char* message) { ++checks; if (!ok) throw std::runtime_error(message); }
void configure(RigidWorld& w) {
    w.params.collideWithDomain=false; w.params.sleeping=false;
    w.params.linearDamping=0; w.params.angularDamping=0;
}
void freeFall(float scale) {
    RigidWorld w; configure(w);
    w.addBox({0,2,0},{.1f,.2f,.3f},Quaternion(),1000,Vector3(1));
    Box3dRigidBackend backend(scale); backend.load(w);
    const float dt=1.0f/600; const int n=120;
    for (int i=0;i<n;++i) backend.step(w,dt);
    const auto& b=w.bodies()[0];
    const Vector3 expectedPos=Vector3(0,2,0)+w.params.gravity*(dt*dt*float(n*(n+1))/2);
    require(length(b.vel-w.params.gravity*(dt*n))<2e-5f,"free-fall velocity differs from g*t");
    require(length(b.pos-expectedPos)<2e-5f,"free-fall position differs from semi-implicit integration");
    require(length(b.angVel)<1e-7f,"free fall created angular motion");
    require(length(b.prevPos-(b.pos-b.vel*dt))<2e-6f,"previous pose was not preserved");
}
void forceAndTorque(float scale) {
    RigidWorld w; configure(w); w.params.gravity=Vector3(0);
    w.addBox({0,0,0},{.1f,.2f,.3f},Quaternion(),1000,Vector3(1));
    auto& b=w.bodies()[0]; const float dt=.001f;
    const Vector3 force(3,-2,1),torque(.1f,.2f,-.3f);
    const Vector3 expected=force*(b.invMass*dt);
    const Vector3 angular(torque.x*b.invInertiaLocal.x*dt,torque.y*b.invInertiaLocal.y*dt,torque.z*b.invInertiaLocal.z*dt);
    Box3dRigidBackend backend(scale); backend.load(w); b.force=force; b.torque=torque;
    backend.step(w,dt);
    require(length(b.vel-expected)<1e-7f,"scaled force did not preserve F/m acceleration");
    require(length(b.angVel-angular)<1e-6f,"scaled torque did not preserve angular acceleration");
    require(length(b.force)+length(b.torque)==0,"loads were not consumed");
    const auto velocity=b.vel; backend.step(w,dt);
    require(length(b.vel-velocity)<1e-7f,"force was applied twice");
}
void momentum(float scale) {
    RigidWorld w; configure(w); w.params.gravity=Vector3(0);
    w.addBox({-.2f,0,0},Vector3(.1f),Quaternion(),1000,Vector3(1));
    w.addBox({.2f,0,0},Vector3(.1f),Quaternion(),2000,Vector3(1));
    for (auto& b:w.bodies()) b.friction=b.restitution=0;
    w.bodies()[0].vel={1,0,0}; w.bodies()[1].vel={-1,0,0};
    const auto total=[&w] { Vector3 p(0); for(const auto& b:w.bodies()) p+=b.vel*b.mass; return p; };
    const auto initial=total(); Box3dRigidBackend backend(scale); backend.load(w); size_t contacts=0;
    for(int i=0;i<240;++i) { backend.step(w,1.0f/600); contacts=std::max(contacts,backend.contactPointCount()); }
    require(contacts>0,"head-on boxes never produced contacts");
    require(length(total()-initial)<1e-4f,"contact response changed total linear momentum");
    std::cout << "collision scale " << numberText(scale) << ": velocities " << numberText(w.bodies()[0].vel.x)
        << '/' << numberText(w.bodies()[1].vel.x) << ", contact depth " << numberText(backend.deepestContact()) << '\n';
    // Soft contact bias may leave a small separating speed even with restitution zero.
    // Check response against the common-velocity solution with an explicit 0.1 m/s residual bound.
    require(std::fabs(w.bodies()[0].vel.x+1.0f/3)<.1f && std::fabs(w.bodies()[1].vel.x+1.0f/3)<.1f,
        "inelastic unequal-mass collision exceeds allowed soft-contact residual");
}
float compoundOverlap(const RigidWorld& w) {
    float worst=0;
    const auto& a=w.bodies()[0]; const auto& b=w.bodies()[1];
    for(const auto& ca:static_cast<const CompoundShape&>(*a.shape).children())
        for(const auto& cb:static_cast<const CompoundShape&>(*b.shape).children()) {
            const PosedShape pa{ca.shape.get(),a.rotation()*ca.R,a.pos+a.rotation()*ca.t};
            const PosedShape pb{cb.shape.get(),b.rotation()*cb.R,b.pos+b.rotation()*cb.t};
            const auto sat=rf::experimental::polytopeSat(pa,pb,0,true);
            require(sat.valid,"cold torus SAT audit failed");
            worst=std::max(worst,float(-sat.gap));
        }
    return worst;
}
void coldTori(float scale) {
    RigidWorld w; configure(w); w.params.gravity=Vector3(0);
    const auto ring=ringShape();
    w.addBody(ring,{.49678934f,1.689123f,.86407435f},{.87695223f,.1506506f,.07180688f,-.45066956f},7800,Vector3(1));
    w.addBody(ring,{.5551778f,1.5995147f,.8649265f},{.5809721f,.6733655f,-.35156918f,-.2923175f},7800,Vector3(1));
    const float initial=compoundOverlap(w); Box3dRigidBackend backend(scale); backend.load(w);
    backend.step(w,1.0f/600); const size_t first=backend.contactPointCount();
    for(int i=1;i<120;++i) backend.step(w,1.0f/600);
    const float final=compoundOverlap(w);
    std::cout << "cold tori scale " << numberText(scale) << ": initial/final m " << numberText(initial)
        << '/' << numberText(final) << ", first contacts " << first << '\n';
    require(initial>.013f && initial<.015f,"cold fixture has wrong initial overlap");
    require(first>0,"native backend lost contacts on stationary intersecting tori");
    require(final<.5f*initial,"native backend did not substantially resolve cold overlap");
    require(backend.supportImportError()<5e-6f,"torus import changed hull support");
}
void invalidInputs() {
    RigidWorld w; configure(w); w.addSphere({0,0,0},.1f,1000,Vector3(1));
    Box3dRigidBackend backend; backend.load(w);
    w.bodies()[0].mass*=2;
    bool rejected=false; try { backend.step(w,.001f); } catch(const std::runtime_error&) { rejected=true; }
    require(rejected,"mass mutation was silently ignored");
    w.bodies()[0].mass/=2; w.addSphere({1,0,0},.1f,1000,Vector3(1));
    rejected=false; try { backend.step(w,.001f); } catch(const std::runtime_error&) { rejected=true; }
    require(rejected,"body-list mutation was silently ignored");
    w.addBallJoint(0,1,{.5f,0,0});
    rejected=false; try { backend.load(w); } catch(const std::runtime_error&) { rejected=true; }
    require(rejected,"unsupported joints were silently ignored");
}
}
int main(int argc,char** argv) {
    try {
        if(argc==2 && std::string(argv[1])=="--cold") {
            int failures=0;
            for(float scale:{1.0f,10.0f}) {
                try { coldTori(scale); } catch(const std::exception& e) { ++failures; std::cerr << "cold FAIL: " << e.what() << '\n'; }
            }
            std::cout << "cold fixture failures " << failures << '\n'; return failures ? 1 : 0;
        }
        if(argc!=1) throw std::runtime_error("optional argument: --cold");
        for(float scale:{1.0f,10.0f}) { freeFall(scale); forceAndTorque(scale); momentum(scale); }
        invalidInputs(); std::cout << checks << " backend checks PASS\n"; return 0;
    } catch(const std::exception& e) { std::cerr << "FAIL: " << e.what() << '\n'; return 1; }
}
