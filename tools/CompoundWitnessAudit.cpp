// Audit compound contacts against actual convex surface planes and exact line/solid intervals.
// This observer applies no impulses. A line interval supplies real surface witnesses, not a new collision solver.
#include "samples/Models.h"
#include "rigid/RigidWorld.h"
#include "core/Format.h"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
using namespace rf;
struct Interval { double lo=-1e100,hi=1e100; bool valid=true; };
double scalar(Vector3 a,Vector3 b) { return double(a.x)*b.x+double(a.y)*b.y+double(a.z)*b.z; }
// Intersect q+t*n with all outward half-spaces of the solid; exact algebra in ordinary floating point.
Interval lineInterval(const PosedShape& p,Vector3 q,Vector3 n) {
    const auto x=p.R.transposed()*(q-p.p),d=p.R.transposed()*n;
    Interval result;
    auto plane=[&](Vector3 normal,Vector3 vertex) {
        const double speed=scalar(normal,d),gap=scalar(normal,vertex-x);
        if(std::fabs(speed)<1e-12) { if(gap<-1e-7) result.valid=false; return; }
        if(speed>0) result.hi=std::min(result.hi,gap/speed); else result.lo=std::max(result.lo,gap/speed);
    };
    if(p.shape->type()==ShapeType::Box) {
        const auto h=static_cast<const BoxShape&>(*p.shape).halfExtents();
        for(int axis=0;axis<3;++axis) for(int sign:{-1,1}) { Vector3 normal(0),v(0); normal[axis]=float(sign); v[axis]=sign*h[axis]; plane(normal,v); }
    } else {
        const auto& hull=static_cast<const ConvexHullShape&>(*p.shape); const auto& mesh=*hull.mesh();
        for(size_t f=0;f<mesh.triangles.size();++f) plane(mesh.faceNormal(f),hull.vertices()[mesh.triangles[f][0]]);
    }
    result.valid=result.valid && result.lo<=result.hi; return result;
}
Interval unionComponent(const std::vector<PosedShape>& parts,size_t child,Vector3 q,Vector3 n) {
    auto seed=lineInterval(parts[child],q,n); if(!seed.valid) return seed;
    std::vector<Interval> intervals; for(const auto& p:parts) { const auto r=lineInterval(p,q,n); if(r.valid) intervals.push_back(r); }
    std::sort(intervals.begin(),intervals.end(),[](Interval a,Interval b) { return a.lo<b.lo; });
    std::vector<Interval> merged;
    for(const auto r:intervals) {
        if(merged.empty() || r.lo>merged.back().hi+1e-7) merged.push_back(r);
        else merged.back().hi=std::max(merged.back().hi,r.hi);
    }
    for(const auto r:merged) if(r.lo<=seed.lo+1e-7 && r.hi>=seed.hi-1e-7) return r;
    throw std::runtime_error("source interval missing from union");
}
float residual(const PosedShape& p,Vector3 x) { Vector3 normal; return p.shape->signedDistance(p.R.transposed()*(x-p.p),normal); }
bool inOther(const std::vector<PosedShape>& parts,size_t skip,Vector3 x,float interior=0) {
    for(size_t k=0;k<parts.size();++k) if(k!=skip && residual(parts[k],x)<=-interior) return true;
    return false;
}
float coneError(const PosedShape& p,Vector3 x,Vector3 outward) {
    // m in N_K(x) iff m.(v-x)<=0 for every vertex v of convex K.
    return dot(outward,p.support(outward)-x);
}
std::vector<PosedShape> expand(const RigidBody& b) {
    std::vector<PosedShape> result;
    for(const auto& c:static_cast<const CompoundShape&>(*b.shape).children()) result.push_back({c.shape.get(),b.rotation()*c.R,b.pos+b.rotation()*c.t});
    return result;
}
void analyticChecks() {
    BoxShape shape({1,1,1}); const PosedShape p{&shape,Matrix3x3::identity(),Vector3(0)};
    const auto inside=lineInterval(p,{0,0,0},{1,0,0});
    const auto outside=lineInterval(p,{0,2,0},{1,0,0});
    if(!inside.valid || std::fabs(inside.lo+1)>1e-12 || std::fabs(inside.hi-1)>1e-12 || outside.valid)
        throw std::runtime_error("analytic line/box test failed");
    const PosedShape adjacent{&shape,Matrix3x3::identity(),{2,0,0}};
    const auto merged=unionComponent({p,adjacent},0,{0,0,0},{1,0,0});
    if(!merged.valid || std::fabs(merged.lo+1)>1e-12 || std::fabs(merged.hi-3)>1e-12)
        throw std::runtime_error("analytic touching union test failed");
    const PosedShape gap{&shape,Matrix3x3::identity(),{3,0,0}};
    if(std::fabs(unionComponent({p,gap},0,{0,0,0},{1,0,0}).hi-1)>1e-12)
        throw std::runtime_error("line union filled a non-convex gap");
}
void audit(const RigidWorld& world,std::ostream& csv,std::ostream& summary) {
    const auto a=expand(world.bodies()[0]),b=expand(world.bodies()[1]); NarrowPhase narrow;
    int points=0,manifolds=0,seamManifolds=0,individualExternal=0,badAnchors=0,rays=0,externalRays=0;
    float maxResidual=0,maxRayShift=0;
    for(size_t u=0;u<a.size();++u) for(size_t v=0;v<b.size();++v) {
        ContactManifold m; if(!narrow.collide(a[u],b[v],m) || m.points.empty()) continue;
        ++manifolds;
        const auto& first=m.points[0]; const auto n=first.normal;
        const bool discarded=inOther(a,u,first.position-n*(.5f*first.depth+.0005f)) || inOther(b,v,first.position+n*(.5f*first.depth+.0005f));
        seamManifolds+=discarded;
        for(size_t i=0;i<m.points.size();++i) {
            const auto& c=m.points[i]; const auto na=c.normal;
            const auto onA=c.position-na*(.5f*c.depth),onB=c.position+na*(.5f*c.depth);
            const float ra=residual(a[u],onA),rb=residual(b[v],onB);
            const bool external=!inOther(a,u,onA-na*.0005f) && !inOther(b,v,onB+na*.0005f);
            const auto ia=lineInterval(a[u],c.position,na),ib=lineInterval(b[v],c.position,na);
            ++points; individualExternal+=external;
            maxResidual=std::max(maxResidual,std::max(std::fabs(ra),std::fabs(rb))); badAnchors+=std::fabs(ra)>1e-5f || std::fabs(rb)>1e-5f;
            csv << u << ',' << v << ',' << i << ',' << discarded << ',' << external << ',' << numberText(c.depth)
                << ',' << numberText(ra) << ',' << numberText(rb) << ',' << (ia.valid && ib.valid);
            if(ia.valid && ib.valid) {
                ++rays; const auto pa=c.position+na*float(ia.lo),pb=c.position+na*float(ib.hi);
                const float depth=float(ib.hi-ia.lo); maxRayShift=std::max(maxRayShift,std::fabs(depth-c.depth));
                externalRays+=!inOther(a,u,pa,1e-6f) && !inOther(b,v,pb,1e-6f);
                const auto ua=unionComponent(a,u,c.position,na),ub=unionComponent(b,v,c.position,na);
                csv << ',' << numberText(depth) << ',' << numberText(residual(a[u],pa)) << ',' << numberText(residual(b[v],pb))
                    << ',' << numberText(coneError(a[u],pa,-na)) << ',' << numberText(coneError(b[v],pb,na)) << ',' << numberText(ub.hi-ua.lo);
            } else csv << ",,,,,,";
            csv << '\n';
        }
    }
    summary << "manifolds " << manifolds << "\npoints " << points << "\nseam_manifolds " << seamManifolds
        << "\nindividually_external_original_probe " << individualExternal << "\nsynthetic_anchors_off_surface_over_10um " << badAnchors
        << "\nmax_surface_plane_residual_m " << numberText(maxResidual) << "\nvalid_line_witness_pairs " << rays
        << "\nline_witness_pairs_not_inside_other_children " << externalRays << "\nmax_assigned_vs_line_depth_m " << numberText(maxRayShift) << '\n';
}
}
int main(int argc,char** argv) {
    try {
        if(argc!=2 || !std::filesystem::create_directory(argv[1])) throw std::runtime_error("NEW_DIRECTORY required");
        analyticChecks(); RigidWorld world; const auto ring=ringShape();
        world.addBody(ring,{.49678934f,1.689123f,.86407435f},{.87695223f,.1506506f,.07180688f,-.45066956f},7800,Vector3(1));
        world.addBody(ring,{.5551778f,1.5995147f,.8649265f},{.5809721f,.6733655f,-.35156918f,-.2923175f},7800,Vector3(1));
        std::ofstream csv(std::filesystem::path(argv[1])/"witnesses.csv"),summary(std::filesystem::path(argv[1])/"summary.txt");
        if(!csv || !summary) throw std::runtime_error("cannot write audit");
        csv << "child_a,child_b,point,discarded_manifold,external_original_probe,assigned_depth_m,residual_a_m,residual_b_m,valid_line_pair,line_depth_m,line_residual_a_m,line_residual_b_m,normal_cone_error_a_m,normal_cone_error_b_m,union_line_depth_m\n";
        audit(world,csv,summary); std::cout << "Analytical line/union checks PASS; witness audit written; no solver change.\n"; return 0;
    } catch(const std::exception& e) { std::cerr << e.what() << '\n'; return 2; }
}
