// Analytical checks for the experimental polytope contact/CCD path, independent of torus success.
// Unit tetrahedron x,y,z >= 0, x+y+z <= 1; restore its model frame after hull mass centering.
#include "SatCcd.h"
#include "samples/Samples.h"
#include <iostream>
#include <stdexcept>

using namespace rf;
using namespace rf::experimental;
namespace {
int checks=0;
void require(bool yes,const char* label) { ++checks; if (!yes) throw std::runtime_error(label); }
TriMesh tetrahedron() {
    TriMesh mesh;
    mesh.positions={{0,0,0},{1,0,0},{0,1,0},{0,0,1}};
    mesh.triangles={{0,2,1},{0,1,3},{0,3,2},{1,2,3}};
    return mesh;
}
PosedShape posed(const ConvexHullShape& shape,Vector3 translation={}) {
    return {&shape,shape.principalRotation(),shape.centerOfMass()+translation};
}
SweptPose sweep(const ConvexHullShape& shape,Vector3 from,Vector3 to) {
    SweptPose s;
    s.shape=&shape; s.p0=shape.centerOfMass()+from; s.p1=shape.centerOfMass()+to;
    s.partR=shape.principalRotation(); return s;
}
void tetraChecks() {
    ConvexHullShape t(tetrahedron());
    const auto a=posed(t);
    auto far=polytopeSat(a,posed(t,{.5f,.5f,.5f}));
    require(far.valid && std::abs(far.gap-.5/std::sqrt(3.0))<1e-6,"tetra diagonal separation");
    auto overlap=polytopeSat(a,posed(t,{.75f,0,0}));
    require(std::abs(overlap.gap+.25/std::sqrt(3.0))<1e-6,"tetra analytical penetration depth");
    auto reverse=polytopeSat(posed(t,{.75f,0,0}),a);
    require(std::abs(reverse.gap-overlap.gap)<1e-7 && dot(reverse.normal,overlap.normal)<-.999f,"SAT swap symmetry");
    Matrix3x3 r=Quaternion::fromAxisAngle({.3f,.5f,.7f},.7f).toMatrix3x3();
    auto b=posed(t,{.75f,0,0}), c=a;
    for (auto* p : {&b,&c}) { p->R=r*p->R; p->p=r*p->p+Vector3(3,2,-1); }
    require(std::abs(polytopeSat(c,b).gap-overlap.gap)<2e-6,"common rigid transform");
    auto hit=sweptSat(sweep(t,{-2,0,0},{2,0,0}),sweep(t,{},{ }),1e-5f);
    require(hit.status==SweepStatus::Hit && std::abs(hit.s-.25f)<1e-4f,"tetra swept translation exact t=1/4");
    ContactManifold m;
    require(collideSat(a,posed(t,{.75f,0,0}),m) && !m.points.empty(),"tetra SAT contact manifold");
    float deepest=-kInf;
    for (const auto& p : m.points) deepest=std::max(deepest,p.depth);
    require(std::abs(deepest-.25f/std::sqrt(3.0f))<1e-5f,"tetra deepest contact matches analytical depth");
}
void boxChecks() {
    BoxShape large(Vector3(1)), small(Vector3(.25f));
    PosedShape a{&large,{},{}}, b{&small,{},{}};
    require(std::abs(polytopeSat(a,b).gap+1.25)<1e-12,"nested interval escape depth");
    b={&large,{},Vector3(1.9995f,0,0)};
    require(std::abs(polytopeSat(a,b).gap+.0005)<1e-7,"half millimetre penetration");
    Matrix3x3 rotation;
    const Vector3 columns[]={{.475433528f,.400452136f,-.783326910f},
        {-.644217687f,.764842187f,0},{.599121467f,.504633050f,.621609968f}};
    for (int i=0; i<3; ++i) for (int j=0; j<3; ++j) rotation.m[j][i]=columns[i][j];
    b={&large,rotation,{-2.18290062f,.14231013f,-2.09478714f}};
    require(polytopeSat(a,b).gap>.17,"edge-cross-only OBB separating axis");
    BoxShape half(Vector3(.5f));
    Matrix3x3 rz=Quaternion::fromAxisAngle({0,0,1},kPi/4).toMatrix3x3();
    Matrix3x3 rx=Quaternion::fromAxisAngle({1,0,0},kPi/4).toMatrix3x3();
    a={&half,rz,{}}; b={&half,rx,{0,std::sqrt(2.0f)-.05f,0}};
    ContactManifold m;
    require(collideSat(a,b,m) && m.points.size()==1,"edge-edge contact count");
    require(std::abs(m.points[0].depth-.05f)<1e-5f && length(m.points[0].position-Vector3(0,.6821f,0))<.01f,"edge interior witness");
}
void sweepChecks() {
    BoxShape rod({1,.01f,.01f}), target(Vector3(.05f));
    SweptPose a; a.shape=&rod; a.dTheta={0,0,kPi};
    SweptPose b; b.shape=&target; b.p0=b.p1={0,.9f,0};
    require(polytopeSat(a.at(0),b.at(0)).gap>0 && polytopeSat(a.at(1),b.at(1)).gap>0,"rotation separated endpoints");
    require(polytopeSat(a.at(.5f),b.at(.5f)).gap<0,"rotation midpoint overlaps");
    auto r=sweptSat(a,b,1e-5f,128);
    require(r.status==SweepStatus::Hit && r.s>.1f && r.s<.5f,"rotational tunneling caught");
    require(polytopeSat(a.at(r.s),b.at(r.s)).gap>=-1e-6,"CCD stops before overlap");
    require(sweptSat(a,b,1e-5f,1).status==SweepStatus::Unresolved,"iteration exhaustion is unknown");
    require(conservativeToi(a,b,1e-5f,1).hit,"world must not treat unknown as no collision");
    a.dTheta={}; b.p0=b.p1={};
    require(sweptSat(a,b,1e-5f).status==SweepStatus::InitialContact,"initial overlap belongs to discrete contact");
    b.p0=b.p1={0,2,0};
    require(sweptSat(a,b,1e-5f).status==SweepStatus::Clear,"zero motion separated");
}
void candidateChecks() {
    TriMesh left=primitives::box(Vector3(.1f)), right=left;
    left.translate({-5,0,0});
    CompoundShape two({left,right},primitives::merge({left,right}));
    BoxShape moving(Vector3(.1f));
    SweptPose a; a.shape=&moving; a.p0={-3,0,0}; a.p1={3,0,0};
    SweptPose b; b.shape=&two; b.p0=b.p1=two.centerOfMass(); b.q0=Quaternion::fromMatrix3x3(two.principalRotation());
    auto result=compoundToi(a,b,1e-5f);
    require(result.hit && std::abs(result.s-2.8f/6)<1e-4f,"farther current child has the future impact");
    b.p0.y+=2; b.p1.y+=2;
    require(!compoundToi(a,b,1e-5f).hit,"support-plane candidate culling is safe");
}
} // namespace
void edgePruningChecks() {
    Simulation sim; loadSample(sim,Preset::TorusChains);
    const auto& body=sim.rigid.bodies()[67];
    const auto& children=static_cast<const CompoundShape&>(*body.shape).children();
    for (int i=0; i<64; ++i) {
        const auto& ca=children[size_t(i)%children.size()];
        const auto& cb=children[size_t(i*3)%children.size()];
        PosedShape a{ca.shape.get(),body.rotation()*ca.R,body.pos+body.rotation()*ca.t}, b=a;
        b.shape=cb.shape.get();
        b.R=Quaternion::fromAxisAngle({.3f,.5f,.7f},i*.131f).toMatrix3x3()*cb.R;
        b.p=a.p+Vector3(.008f*rf::sin(i*.71f),.015f*rf::cos(i*.37f),.035f*rf::sin(i*.17f));
        const auto fast=polytopeSat(a,b), full=polytopeSat(a,b,1e100,true);
        require((fast.gap>1e-7)==(full.gap>1e-7),"edge cone pruning preserves separation");
        if (full.gap<0) require(std::abs(fast.gap-full.gap)<1e-7,"edge cone pruning preserves convex overlap depth");
    }
}
int main() {
    try {
        tetraChecks(); boxChecks(); sweepChecks(); candidateChecks(); edgePruningChecks();
        std::cout << checks << " analytical checks PASS\n"; writeStats(std::cout); return 0;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
