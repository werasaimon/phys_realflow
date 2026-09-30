// Finite-difference check of contact gap vs solver normal for an actual point/plane manifold.
// Also checks whether EPA depth bracketing moves the reconstructed witnesses off the surface.
#include "rigid/NarrowPhase.h"
#include "core/Format.h"
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
using namespace rf;
float gap(const PosedShape& a,const PosedShape& b,Vector3 normal) {
    std::vector<ContactPoint> points;
    if(!NarrowPhase::faceManifold(a,b,normal,-kInf,kInf,points) || points.empty()) throw std::runtime_error("no plane manifold");
    return -points[0].depth;
}
void check(bool requireConsistent) {
    BoxShape upper({.2f,.1f,.2f}),lower({1,.1f,1});
    const PosedShape a{&upper,Matrix3x3::identity(),{0,.198f,0}},b{&lower,Matrix3x3::identity(),Vector3(0)};
    const Vector3 input=normalize(Vector3(.01f,1,0));
    std::vector<ContactPoint> points; NarrowPhase::faceManifold(a,b,input,-kInf,kInf,points);
    if(points.empty()) throw std::runtime_error("empty manifold");
    Vector3 gradient(0); constexpr float eps=1e-4f;
    for(int axis=0;axis<3;++axis) {
        auto plus=a,minus=a; plus.p[axis]+=eps; minus.p[axis]-=eps;
        gradient[axis]=(gap(plus,b,input)-gap(minus,b,input))/(2*eps);
    }
    const float error=length(gradient-points[0].normal);
    std::cout << "normal vs gap-gradient error " << numberText(error) << '\n';
    if(requireConsistent && error>2e-4f) throw std::runtime_error("contact normal differs from gap Jacobian");
    points.clear(); NarrowPhase::faceManifold(a,b,input,.0028f,.0029f,points);
    float residual=0;
    for(const auto& c:points) {
        const auto pa=c.position-c.normal*(.5f*c.depth),pb=c.position+c.normal*(.5f*c.depth);
        Vector3 ignored;
        residual=std::max(residual,std::fabs(upper.signedDistance(pa-a.p,ignored)));
        residual=std::max(residual,std::fabs(lower.signedDistance(pb-b.p,ignored)));
    }
    std::cout << "bracketed manifold max surface residual m " << numberText(residual) << '\n';
    if(requireConsistent && residual>1e-5f) throw std::runtime_error("depth bracket created synthetic surface points");
}
}
int main(int argc,char** argv) {
    try {
        const bool expect=argc==2 && std::string(argv[1])=="--expect-consistent";
        if(argc>2 || (argc==2 && !expect)) throw std::runtime_error("optional argument: --expect-consistent");
        check(expect); std::cout << (expect ? "Face geometry/gradient checks PASS\n" : "Baseline geometry/gradient discrepancies measured\n"); return 0;
    } catch(const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
