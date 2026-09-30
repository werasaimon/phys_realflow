// Independent sampled-pose SAT audit of the actual convex children of sample 39.
// Double precision, triangle face normals and all edge cross products; no GJK/EPA/seam filter.
// SAT theorem: https://www.geometrictools.com/Documentation/MethodOfSeparatingAxes.pdf
// Usage: rf_torus_geometry_audit NEW_DIRECTORY [POSES_CSV [CHECKPOINT [PRESET]]]
// PRESET defaults to 39; use 40 for the short hanging chain.
// Pair CSV: body 67/68 at defect substeps. Body CSV: all pairs at selected frame checkpoints.
// CHECKPOINT=-2 audits every body frame, writing positive overlaps plus per-frame summaries.
#include "samples/Samples.h"
#include "core/Format.h"
#include <array>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace {
struct V {
    double x = 0, y = 0, z = 0;
    V() = default;
    V(double a, double b, double c) : x(a), y(b), z(c) {}
    explicit V(rf::Vector3 p) : x(p.x), y(p.y), z(p.z) {}
    V operator+(V b) const { return {x+b.x,y+b.y,z+b.z}; }
    V operator-(V b) const { return {x-b.x,y-b.y,z-b.z}; }
    V operator*(double s) const { return {x*s,y*s,z*s}; }
};
double dot(V a, V b) { return a.x*b.x+a.y*b.y+a.z*b.z; }
V cross(V a, V b) { return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x}; }
void direction(std::vector<V>& dirs, V v) {
    const double len = std::sqrt(dot(v,v));
    if (len < 1e-14) return;
    v = v*(1/len);
    // Only remove numerically identical directions; do not merge curved facets.
    for (V d : dirs) if (dot(d-v,d-v) < 1e-24 || dot(d+v,d+v) < 1e-24) return;
    dirs.push_back(v);
}
struct Poly {
    std::vector<V> vertices, normals, edges;
    V lo{1e100,1e100,1e100}, hi{-1e100,-1e100,-1e100};
    void bound() {
        for (V p : vertices) {
            lo = {std::min(lo.x,p.x),std::min(lo.y,p.y),std::min(lo.z,p.z)};
            hi = {std::max(hi.x,p.x),std::max(hi.y,p.y),std::max(hi.z,p.z)};
        }
    }
};
V transform(const rf::Matrix3x3& r, V v, V t) {
    // Promote before arithmetic; input rotations/meshes remain SDK float data.
    const V x(r*rf::Vector3(1,0,0)), y(r*rf::Vector3(0,1,0)), z(r*rf::Vector3(0,0,1));
    return t+x*v.x+y*v.y+z*v.z;
}
std::vector<Poly> geometry(const rf::RigidBody& body) {
    const auto& compound = static_cast<const rf::CompoundShape&>(*body.shape);
    std::vector<Poly> result;
    for (const auto& part : compound.children()) {
        const auto* hull = dynamic_cast<const rf::ConvexHullShape*>(part.shape.get());
        if (!hull) throw std::runtime_error("expected hull child");
        Poly p;
        for (auto v : hull->vertices())
            p.vertices.push_back(transform(body.rotation(),transform(part.R,V(v),V(part.t)),V(body.pos)));
        std::set<std::pair<unsigned,unsigned>> edges;
        for (const auto& f : hull->mesh()->triangles) {
            direction(p.normals,cross(p.vertices[f[1]]-p.vertices[f[0]],p.vertices[f[2]]-p.vertices[f[0]]));
            for (int k = 0; k < 3; ++k) edges.emplace(std::min(f[k],f[(k+1)%3]),std::max(f[k],f[(k+1)%3]));
        }
        for (auto edge : edges) direction(p.edges,p.vertices[edge.second]-p.vertices[edge.first]);
        p.bound();
        result.push_back(std::move(p));
    }
    return result;
}
double axisGap(const Poly& a, const Poly& b, V n) {
    const double length = std::sqrt(dot(n,n));
    if (length < 1e-14) return -1e100;
    n = n*(1/length);
    double amin=1e100,amax=-1e100,bmin=1e100,bmax=-1e100;
    for (V p : a.vertices) { double d=dot(p,n); amin=std::min(amin,d); amax=std::max(amax,d); }
    for (V p : b.vertices) { double d=dot(p,n); bmin=std::min(bmin,d); bmax=std::max(bmax,d); }
    return std::max(bmin-amax,amin-bmax);
}
double depth(const Poly& a, const Poly& b) {
    if (a.lo.x>b.hi.x || b.lo.x>a.hi.x || a.lo.y>b.hi.y || b.lo.y>a.hi.y || a.lo.z>b.hi.z || b.lo.z>a.hi.z) return 0;
    double gap=-1e100;
    for (const auto* directions : {&a.normals,&b.normals})
        for (V n : *directions) { gap=std::max(gap,axisGap(a,b,n)); if (gap>=0) return 0; }
    for (V u : a.edges) for (V v : b.edges) {
        gap=std::max(gap,axisGap(a,b,cross(u,v))); if (gap>=0) return 0;
    }
    return -gap;
}
Poly cube(V center, double half) {
    Poly p;
    for (int x : {-1,1}) for (int y : {-1,1}) for (int z : {-1,1})
        p.vertices.push_back(center+V(x*half,y*half,z*half));
    p.normals=p.edges={{1,0,0},{0,1,0},{0,0,1}};
    p.bound(); return p;
}
void edgeAxisTest() {
    const Poly a=cube({},1);
    Poly b=cube({},1);
    const V center(-2.18290062173539,.14231012575111368,-2.0947871410170906);
    const std::vector<V> axes={{.47543352776997644,.4004521361232219,-.7833269096274834},
        {-.644217687237691,.7648421872844885,0}, {.5991214669182833,.5046330500712651,.6216099682706644}};
    for (auto& v : b.vertices) v=center+axes[0]*v.x+axes[1]*v.y+axes[2]*v.z;
    b.normals=b.edges=axes;
    b.lo={1e100,1e100,1e100}; b.hi={-1e100,-1e100,-1e100}; b.bound();
    for (const auto* ns : std::array<const std::vector<V>*,2>{&a.normals,&b.normals}) for (V n : *ns)
        if (axisGap(a,b,n)>=0) throw std::runtime_error("edge-only fixture has face separation");
    // Projection of OBB centers and half extents gives >0.17 m on this edge-cross axis.
    if (axisGap(a,b,cross(a.edges[1],b.edges[0]))<.17 || depth(a,b)!=0)
        throw std::runtime_error("edge-cross separating axis missing");
}
void selfTest() {
    edgeAxisTest();
    const auto a=cube({0,0,0},1);
    for (const auto& test : std::vector<std::pair<Poly,double>>{
            {cube({3,0,0},1),0}, {cube({2,0,0},1),0}, {cube({1.75,0,0},1),.25},
            {cube({0,0,0},.25),1.25}, {cube({1.9995,0,0},1),.0005}})
        if (std::abs(depth(a,test.first)-test.second)>1e-12) throw std::runtime_error("analytic SAT check failed");
    // A common rigid transform must preserve the measured overlap.
    Poly b=cube({0,0,0},1), c=b;
    const rf::Matrix3x3 r=rf::Quaternion::fromAxisAngle({.3f,.5f,.7f},.7f).toMatrix3x3();
    for (auto& v : b.vertices) v=transform(r,v,V(3,2,-1));
    for (auto& v : c.vertices) v=transform(r,v+V(1.75,0,0),V(3,2,-1));
    for (auto& n : b.normals) n=transform(r,n,{});
    for (auto& n : c.normals) n=transform(r,n,{});
    b.edges=b.normals; c.edges=c.normals;
    b.lo=c.lo={1e100,1e100,1e100}; b.hi=c.hi={-1e100,-1e100,-1e100}; b.bound(); c.bound();
    if (std::abs(depth(b,c)-.25)>1e-6) throw std::runtime_error("transformed SAT check failed");
}
double row(std::ofstream& csv, int step, int a, int b, const std::vector<Poly>& pa, const std::vector<Poly>& pb, bool positiveOnly=false) {
    double worst=0; int overlaps=0;
    for (const auto& x : pa) for (const auto& y : pb) {
        const double d=depth(x,y); worst=std::max(worst,d); overlaps+=d>1e-7;
    }
    if (positiveOnly && worst==0) return 0;
    csv << step << ',' << a << ',' << b << ',' << overlaps << ',' << rf::numberText(worst)
        << ',' << (worst>.0001+1e-7) << ',' << (worst>.0005+1e-7) << '\n';
    csv.flush();
    return worst;
}
std::vector<float> numbers(std::string line) {
    std::istringstream in(line); std::string cell; std::vector<float> values;
    while (std::getline(in,cell,',')) {
        float v=0; if (!rf::parseNumber(cell,v) || !std::isfinite(v)) throw std::runtime_error("invalid CSV number");
        values.push_back(v);
    }
    return values;
}
void bodyCheckpoint(int frame,rf::Simulation& sim,std::ofstream& csv,std::ofstream& summary,bool positiveOnly) {
    std::vector<std::vector<Poly>> shapes;
    for(const auto& b:sim.rigid.bodies()) shapes.push_back(geometry(b));
    size_t checked=0,overlaps=0; double maximum=0;
    for(size_t a=0;a<shapes.size();++a) for(size_t b=a+1;b<shapes.size();++b) {
        const double d=row(csv,frame,int(a),int(b),shapes[a],shapes[b],positiveOnly);
        ++checked; overlaps+=d>0; maximum=std::max(maximum,d);
    }
    summary << frame << ',' << checked << ',' << overlaps << ',' << rf::numberText(maximum) << '\n';
    summary.flush();
}
void replayBodies(std::istream& input,rf::Simulation& sim,const std::filesystem::path& output,int extra) {
    std::ofstream csv(output/"body-audit.csv");
    std::ofstream summary(output/"body-checkpoints.csv");
    if(!csv || !summary) throw std::runtime_error("cannot write body audit");
    summary << "frame,pairs_checked,pairs_overlapping,max_child_sat_depth_m\n";
    csv << "frame,body_a,body_b,overlapping_child_pairs,max_child_sat_depth_m,over_0_1mm,over_0_5mm\n";
    int frame=-1,used=0; std::set<int> seen;
    auto finish=[&] {
        if(frame<0) return;
        if(seen.size()!=sim.rigid.bodies().size()) throw std::runtime_error("incomplete body checkpoint");
        if(extra==-2 || (extra>=0 ? frame==extra : std::set<int>{1,30,60,120,300,600}.count(frame)!=0)) {
            bodyCheckpoint(frame,sim,csv,summary,extra==-2); ++used;
        }
    };
    std::string line;
    while(std::getline(input,line)) {
        const auto v=numbers(line);
        if(v.size()!=15 || v[0]!=float(int(v[0])) || v[1]!=float(int(v[1]))) throw std::runtime_error("invalid body pose row");
        const int next=int(v[0]),index=int(v[1]);
        if(index<0 || size_t(index)>=sim.rigid.bodies().size()) throw std::runtime_error("invalid body index");
        if(next!=frame) { finish(); if(next<=frame) throw std::runtime_error("non-monotonic frames"); frame=next; seen.clear(); }
        if(!seen.insert(index).second) throw std::runtime_error("duplicate body pose");
        auto& b=sim.rigid.bodies()[size_t(index)]; b.pos={v[2],v[3],v[4]}; b.rot={v[5],v[6],v[7],v[8]};
    }
    finish(); if(!used) throw std::runtime_error("no body checkpoints found");
}
void replay(const char* path, rf::Simulation& sim, std::ofstream& csv, int extra,const std::filesystem::path& output) {
    std::ifstream input(path); if (!input) throw std::runtime_error("cannot read pose CSV");
    std::string line; std::getline(input,line);
    if(line.rfind("frame,body,x,y,z,qw,qx,qy,qz,",0)==0) { replayBodies(input,sim,output,extra); return; }
    if (line.rfind("substep,ax,ay,az,aqw,aqx,aqy,aqz,",0)!=0) throw std::runtime_error("unexpected pose CSV schema");
    int used=0;
    while (std::getline(input,line)) {
        auto v=numbers(line);
        if (v.size()!=27) throw std::runtime_error("expected 27 pose fields");
        int step=int(v[0]);
        if (step!=extra && !std::set<int>{1,371,372,373,460,461,462,489,490,491,492}.count(step)) continue;
        for (int i=0; i<2; ++i) {
            auto& b=sim.rigid.bodies()[size_t(67+i)]; size_t o=size_t(1+13*i);
            b.pos={v[o],v[o+1],v[o+2]}; b.rot={v[o+3],v[o+4],v[o+5],v[o+6]};
        }
        row(csv,step,67,68,geometry(sim.rigid.bodies()[67]),geometry(sim.rigid.bodies()[68])); ++used;
    }
    if (!used) throw std::runtime_error("no checkpoints found");
}
}
int main(int argc, char** argv) {
    try {
        if (argc<2 || argc>5) throw std::runtime_error("NEW_DIRECTORY [POSES_CSV [CHECKPOINT [PRESET]]]");
        selfTest();
        if (!std::filesystem::create_directory(argv[1])) throw std::runtime_error("output exists");
        std::ofstream csv(std::filesystem::path(argv[1])/"audit.csv");
        if (!csv) throw std::runtime_error("cannot write audit");
        csv << "substep,body_a,body_b,overlapping_child_pairs,max_child_sat_depth_m,over_0_1mm,over_0_5mm\n";
        const int preset=argc==5 ? std::stoi(argv[4]) : int(rf::Preset::TorusChains);
        if(preset!=int(rf::Preset::TorusChains) && preset!=int(rf::Preset::HangingTorus))
            throw std::runtime_error("expected torus-chain preset 39 or 40");
        rf::Simulation sim; rf::loadSample(sim,rf::Preset(preset));
        std::vector<std::vector<Poly>> shapes;
        for (const auto& b : sim.rigid.bodies()) shapes.push_back(geometry(b));
        for (size_t a=0; a<shapes.size(); ++a) for (size_t b=a+1; b<shapes.size(); ++b)
            row(csv,0,int(a),int(b),shapes[a],shapes[b]);
        if (argc>=3) replay(argv[2],sim,csv,argc>=4 ? std::stoi(argv[3]) : -1,argv[1]);
        std::cout << "Analytic checks PASS; sampled geometry written. This is not swept CCD or an exact union penetration depth.\n";
        return 0;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 2; }
}
