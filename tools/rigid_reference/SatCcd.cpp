// Experimental SAT uses actual convex mesh facets and edge cross products, including tetrahedra.
// GJK orders compound pairs, but only a support-plane lower bound may reject swept candidates.
// Conservative advancement recomputes geometry during rotation; a finite budget is not "clear".
#include "SatCcd.h"
#include "core/Format.h"
#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <map>
#include <ostream>
#include <set>

namespace rf::experimental {
namespace {
struct Vec {
    double x = 0, y = 0, z = 0;
    Vec() = default;
    Vec(double a, double b, double c) : x(a), y(b), z(c) {}
    explicit Vec(Vector3 v) : x(v.x), y(v.y), z(v.z) {}
    Vec operator+(Vec b) const { return {x+b.x,y+b.y,z+b.z}; }
    Vec operator-(Vec b) const { return {x-b.x,y-b.y,z-b.z}; }
    Vec operator*(double s) const { return {x*s,y*s,z*s}; }
    Vector3 single() const { return {float(x),float(y),float(z)}; }
};
double dot(Vec a, Vec b) { return a.x*b.x+a.y*b.y+a.z*b.z; }
Vec cross(Vec a, Vec b) { return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x}; }
Vec turn(const Matrix3x3& r, Vec v) {
    return {r.m[0][0]*v.x+r.m[0][1]*v.y+r.m[0][2]*v.z,
            r.m[1][0]*v.x+r.m[1][1]*v.y+r.m[1][2]*v.z,
            r.m[2][0]*v.x+r.m[2][1]*v.y+r.m[2][2]*v.z};
}
void direction(std::vector<Vec>& out, Vec v) {
    double len = std::sqrt(dot(v,v));
    if (len < 1e-14) return;
    v = v*(1/len);
    for (Vec n : out) if (dot(n-v,n-v)<1e-24 || dot(n+v,n+v)<1e-24) return;
    out.push_back(v); // retain distinct curved facets, not a coarse angular merge
}
struct Edge { unsigned a=0,b=0,faceA=0,faceB=0; bool closed=false; };
struct Cone { Vec direction,sideA,sideB; bool closed=false,flat=false; };
struct Geometry {
    std::vector<Vec> vertices, normals, edges, faceNormals;
    std::vector<std::array<uint32_t,3>> triangles;
    std::vector<Edge> adjacency;
    std::vector<Cone> cones;
};
struct Counters {
    std::atomic<uint64_t> sat{0}, axes{0}, narrow{0}, rejected{0}, single{0};
    std::atomic<uint64_t> sweeps{0}, iterations{0}, exhausted{0}, initial{0};
    std::atomic<uint64_t> candidates{0}, culled{0};
} counts;

Geometry meshGeometry(const TriMesh& mesh) {
    Geometry g;
    for (Vector3 p : mesh.positions) g.vertices.emplace_back(p);
    g.triangles=mesh.triangles;
    std::map<std::pair<unsigned,unsigned>,std::vector<unsigned>> edges;
    unsigned face=0;
    for (const auto& f : mesh.triangles) {
        direction(g.normals,cross(g.vertices[f[1]]-g.vertices[f[0]],g.vertices[f[2]]-g.vertices[f[0]]));
        for (int i=0; i<3; ++i) edges[{std::min(f[i],f[(i+1)%3]),std::max(f[i],f[(i+1)%3])}].push_back(face);
        ++face;
    }
    // Triangulation diagonals add valid projection axes; retaining them avoids approximate welding.
    for (const auto& e : edges) {
        direction(g.edges,g.vertices[e.first.second]-g.vertices[e.first.first]);
        Edge edge; edge.a=e.first.first; edge.b=e.first.second; edge.faceA=e.second.front();
        edge.closed=e.second.size()==2; edge.faceB=edge.closed ? e.second.back() : edge.faceA;
        g.adjacency.push_back(edge);
    }
    return g;
}
const Geometry& hullGeometry(const ConvexHullShape& hull) {
    using Owner = std::shared_ptr<const TriMesh>;
    // Keep ownership: a reused raw shape pointer must never retrieve another shape's geometry.
    thread_local std::map<Owner,Geometry,std::owner_less<Owner>> cache;
    auto it = cache.find(hull.mesh());
    if (it == cache.end()) it = cache.emplace(hull.mesh(),meshGeometry(*hull.mesh())).first;
    return it->second;
}
Geometry boxGeometry(const BoxShape& box) {
    Geometry g;
    Vector3 h = box.halfExtents();
    for (int x : {-1,1}) for (int y : {-1,1}) for (int z : {-1,1}) g.vertices.emplace_back(x*h.x,y*h.y,z*h.z);
    g.normals = g.edges = {{1,0,0},{0,1,0},{0,0,1}};
    return g;
}
void worldGeometry(const PosedShape& pose, Geometry& out) {
    Geometry box;
    const Geometry* local = nullptr;
    if (pose.shape->type()==ShapeType::Box) { box=boxGeometry(static_cast<const BoxShape&>(*pose.shape)); local=&box; }
    else local=&hullGeometry(static_cast<const ConvexHullShape&>(*pose.shape));
    out.vertices.clear(); out.normals.clear(); out.edges.clear(); out.faceNormals.clear(); out.cones.clear();
    out.triangles=local->triangles; out.adjacency=local->adjacency;
    for (Vec p : local->vertices) out.vertices.push_back(turn(pose.R,p)+Vec(pose.p));
    if (local->triangles.empty()) for (Vec n : local->normals) out.normals.push_back(turn(pose.R,n));
    for (const auto& f : local->triangles) {
        Vec n=cross(out.vertices[f[1]]-out.vertices[f[0]],out.vertices[f[2]]-out.vertices[f[0]]);
        double len=std::sqrt(dot(n,n)); n=n*(1/std::max(len,1e-100));
        out.faceNormals.push_back(n); out.normals.push_back(n);
    }
    for (Vec e : local->edges) out.edges.push_back(turn(pose.R,e));
    for (const Edge& edge : out.adjacency) {
        Cone cone; cone.direction=out.vertices[edge.b]-out.vertices[edge.a]; cone.closed=edge.closed;
        Vec a=out.faceNormals[edge.faceA], b=out.faceNormals[edge.faceB], arc=cross(a,b);
        double len=std::sqrt(dot(arc,arc)); cone.flat=len<1e-14;
        if (!cone.flat) { cone.sideA=cross(arc,a)*(1/len); cone.sideB=cross(b,arc)*(1/len); }
        out.cones.push_back(cone);
    }
}
bool testAxis(const Geometry& a, const Geometry& b, Vec axis, SatResult& result, double reject);
bool edgeSupports(const Cone& cone, Vec n) {
    if (!cone.closed) return true; // unknown topology: retain this candidate
    if (cone.flat) return false; // coplanar triangulation diagonal: its face axes were tested
    // The supporting edge's outward direction lies in the cone of its two adjacent normals.
    // These oriented triple products test that cone without dividing by a nearly-zero Gram determinant.
    // Eberly SAT section 4: only supporting edge pairs contribute Minkowski difference facets.
    return dot(n,cone.sideA)>=-1e-10 && dot(n,cone.sideB)>=-1e-10;
}
bool edgeAxes(const Geometry& a, const Geometry& b, SatResult& result, double reject, bool full) {
    if (full || a.adjacency.empty() || b.adjacency.empty()) {
        for (Vec x : a.edges) for (Vec y : b.edges)
            if (!testAxis(a,b,cross(x,y),result,reject)) return false;
        return true;
    }
    for (const Cone& x : a.cones) for (const Cone& y : b.cones) {
        if ((x.closed && x.flat) || (y.closed && y.flat)) continue;
        Vec n=cross(x.direction,y.direction);
        double len=std::sqrt(dot(n,n));
        if (len<1e-14) continue;
        n=n*(1/len);
        bool support=(edgeSupports(x,n*(-1)) && edgeSupports(y,n)) ||
                     (edgeSupports(x,n) && edgeSupports(y,n*(-1)));
        if (support && !testAxis(a,b,n,result,reject)) return false;
    }
    return true;
}
bool testAxis(const Geometry& a, const Geometry& b, Vec axis, SatResult& result, double reject) {
    double len = std::sqrt(dot(axis,axis));
    if (len<1e-14) return true;
    axis=axis*(1/len);
    double amin=1e100,amax=-1e100,bmin=1e100,bmax=-1e100;
    for (Vec p : a.vertices) { double d=dot(p,axis); amin=std::min(amin,d); amax=std::max(amax,d); }
    for (Vec p : b.vertices) { double d=dot(p,axis); bmin=std::min(bmin,d); bmax=std::max(bmax,d); }
    double ab=bmin-amax, ba=amin-bmax, gap=std::max(ab,ba);
    ++result.axes;
    if (gap>result.gap) { result.gap=gap; result.normal=(axis*(ba>=ab ? 1 : -1)).single(); }
    return gap<=reject;
}
double roundoffGuard(const PosedShape& a, const PosedShape& b) {
    return 2e-7*(1+length(a.p)+length(b.p)+a.shape->boundingRadius()+b.shape->boundingRadius());
}
// Re-evaluate support along a GJK witness direction: the simplex distance alone is an upper bound.
double supportGap(const PosedShape& a, const PosedShape& b, const GjkResult& g) {
    if (g.intersect || g.distance<1e-9f) return 0;
    Vector3 n=normalize(g.pointA-g.pointB);
    return std::max(0.0,double(rf::dot(a.support(-n)-b.support(n),n))-roundoffGuard(a,b));
}
double speedBound(const SweptPose& a, const SweptPose& b) {
    return length((a.p1-a.p0)-(b.p1-b.p0))+a.angularReach()+b.angularReach();
}
void closestSegments(Vector3 a, Vector3 b, Vector3 c, Vector3 d, Vector3& x, Vector3& y) {
    // Ericson, Real-Time Collision Detection, section 5.1.9, with degenerate point segments.
    Vector3 u=b-a,v=d-c,w=a-c;
    float uu=length2(u), vv=length2(v), uv=rf::dot(u,v), uw=rf::dot(u,w), vw=rf::dot(v,w);
    float s=0,t=0;
    if (uu>1e-20f) {
        if (vv<1e-20f) s=clampv(-uw/uu,0.0f,1.0f);
        else {
            float den=uu*vv-uv*uv;
            if (den>1e-20f) s=clampv((uv*vw-uw*vv)/den,0.0f,1.0f);
            t=(uv*s+vw)/vv;
            if (t<0) { t=0; s=clampv(-uw/uu,0.0f,1.0f); }
            if (t>1) { t=1; s=clampv((uv-uw)/uu,0.0f,1.0f); }
        }
    } else if (vv>1e-20f) t=clampv(vw/vv,0.0f,1.0f);
    x=a+u*s; y=c+v*t;
}
std::vector<Vector3> supportVertices(const PosedShape& pose, Vector3 direction) {
    Geometry g;
    worldGeometry(pose,g);
    double best=-1e100;
    for (Vec p : g.vertices) best=std::max(best,dot(p,Vec(direction)));
    std::vector<Vector3> result;
    for (Vec p : g.vertices) if (best-dot(p,Vec(direction))<1e-7) result.push_back(p.single());
    return result;
}
void singleContact(const PosedShape& a, const PosedShape& b, Vector3 n, float depth, ContactManifold& m) {
    const auto fa=supportVertices(a,-n), fb=supportVertices(b,n);
    Vector3 pa=a.support(-n),pb=b.support(n);
    double best=1e100;
    // Segments between support vertices stay on their supporting convex face. For an edge or
    // vertex this is the entire feature, including its interior rather than only its endpoints.
    for (size_t i=0; i<fa.size(); ++i) for (size_t j=i; j<fa.size(); ++j)
        for (size_t k=0; k<fb.size(); ++k) for (size_t l=k; l<fb.size(); ++l) {
            Vector3 x,y; closestSegments(fa[i],fa[j],fb[k],fb[l],x,y);
            double d=length2(x-y);
            if (d<best) { best=d; pa=x; pb=y; }
        }
    ++counts.single; m.add((pa+pb)*.5f,n,depth);
}
void parts(const SweptPose& body, std::vector<SweptPose>& out) {
    if (body.shape->type()!=ShapeType::Compound) { out.push_back(body); return; }
    for (const auto& child : static_cast<const CompoundShape&>(*body.shape).children()) {
        SweptPose p=body; p.shape=child.shape.get(); p.partR=child.R; p.partT=child.t; out.push_back(p);
    }
}
struct Candidate { size_t a=0,b=0; float distance=0; double lowerTime=0; };
std::vector<Candidate> candidates(const std::vector<SweptPose>& a, const std::vector<SweptPose>& b, float tol) {
    std::vector<Candidate> out;
    for (size_t i=0; i<a.size(); ++i) for (size_t j=0; j<b.size(); ++j) {
        ++counts.candidates;
        const PosedShape pa=a[i].at(0), pb=b[j].at(0);
        const GjkResult g=gjk(pa,pb);
        double lower=supportGap(pa,pb,g), speed=speedBound(a[i],b[j]);
        if (lower>speed+tol) { ++counts.culled; continue; }
        double time=speed>1e-12 ? std::max(0.0,(lower-tol)/speed) : 0;
        out.push_back({i,j,g.intersect ? 0 : g.distance,time});
    }
    // Closest current children first. Other potential impacts remain in the list.
    std::stable_sort(out.begin(),out.end(),[](const Candidate& x,const Candidate& y) { return x.distance<y.distance; });
    return out;
}
} // namespace

bool isPolytope(const ConvexShape* s) { return s && (s->type()==ShapeType::Box || s->type()==ShapeType::ConvexHull); }
bool narrowEnabled() { static const bool on=std::getenv("RF_EXPERIMENT_SAT")!=nullptr; return on; }
bool ccdEnabled() { static const bool on=std::getenv("RF_EXPERIMENT_SAT_CCD")!=nullptr; return on; }

SatResult polytopeSat(const PosedShape& a, const PosedShape& b, double reject, bool full) {
    SatResult result;
    if (!isPolytope(a.shape) || !isPolytope(b.shape)) return result;
    thread_local Geometry wa,wb;
    worldGeometry(a,wa); worldGeometry(b,wb);
    result.valid=true;
    auto done=[&]() { ++counts.sat; counts.axes+=uint64_t(result.axes); return result; };
    for (Vec n : wa.normals) if (!testAxis(wa,wb,n,result,reject)) return done();
    for (Vec n : wb.normals) if (!testAxis(wa,wb,n,result,reject)) return done();
    edgeAxes(wa,wb,result,reject,full);
    return done();
}

bool collideSat(const PosedShape& a, const PosedShape& b, ContactManifold& manifold) {
    ++counts.narrow;
    GjkResult g=gjk(a,b,NarrowPhase::margin);
    if (!g.intersect && g.distance>=NarrowPhase::margin) { ++counts.rejected; return false; }
    SatResult sat=polytopeSat(a,b,NarrowPhase::margin);
    if (!sat.valid || sat.gap>=NarrowPhase::margin) return false;
    Vector3 n=sat.normal;
    float depth=float(-sat.gap);
    // Separated shapes use nearest witnesses for speculation; SAT supplies overlap normal/depth.
    if (sat.gap>1e-7 && !g.intersect && g.distance>1e-7f) {
        n=normalize(g.pointA-g.pointB); depth=-g.distance;
    }
    std::vector<ContactPoint> points;
    if (NarrowPhase::faceManifold(a,b,n,depth,depth,points)) {
        reduceManifold(points,4);
        manifold.points.insert(manifold.points.end(),points.begin(),points.end()); return true;
    }
    singleContact(a,b,n,depth,manifold);
    return true;
}

SweepResult sweptSat(const SweptPose& a, const SweptPose& b, float tol, int budget) {
    SweepResult result;
    if (!isPolytope(a.shape) || !isPolytope(b.shape)) { result.status=SweepStatus::Unsupported; return result; }
    ++counts.sweeps;
    const double speed=speedBound(a,b);
    double s=0,lastChecked=0;
    for (int i=0; i<budget; ++i) {
        ++counts.iterations; result.iterations=i+1;
        const PosedShape pa=a.at(float(s)), pb=b.at(float(s));
        lastChecked=s;
        SatResult sat=polytopeSat(pa,pb,tol); // one separated axis is enough for a safe lower bound
        if (sat.gap<=tol) {
            result.s=float(s); result.status=s==0 ? SweepStatus::InitialContact : SweepStatus::Hit;
            if (s==0) ++counts.initial;
            return result;
        }
        // Any separating support plane proves a gap; its closing rate is at most speed.
        GjkResult g=gjk(pa,pb);
        double lower=std::max(sat.gap-roundoffGuard(pa,pb),supportGap(pa,pb,g));
        if (speed<1e-12 || lower>speed*(1-s)+tol) return result;
        double next=s+(lower-.5*tol)/speed;
        if (next>1) return result;
        if (float(next)<=float(s)) break; // numerical/budget failure must not become a false "clear"
        s=next;
    }
    ++counts.exhausted; result.status=SweepStatus::Unresolved; result.s=float(lastChecked); return result;
}

ToiResult conservativeToi(const SweptPose& a, const SweptPose& b, float tol, int budget) {
    SweepResult r=sweptSat(a,b,tol,budget);
    // Existing world accepts only hit/no-hit. Fail conservatively at the last checked safe pose.
    return {r.status==SweepStatus::Hit || r.status==SweepStatus::Unresolved,r.s,r.iterations};
}

ToiResult compoundToi(const SweptPose& a, const SweptPose& b, float tol) {
    std::vector<SweptPose> pa,pb;
    parts(a,pa); parts(b,pb);
    ToiResult first;
    for (const Candidate& c : candidates(pa,pb,tol)) {
        if (c.lowerTime>=first.s) { ++counts.culled; continue; }
        ToiResult r=(isPolytope(pa[c.a].shape) && isPolytope(pb[c.b].shape))
            ? conservativeToi(pa[c.a],pb[c.b],tol) : timeOfImpact(pa[c.a],pb[c.b],tol);
        if (r.hit && r.s<first.s) first=r;
    }
    return first;
}

void writeStats(std::ostream& out) {
    out << "sat_queries " << counts.sat << "\nsat_axes " << counts.axes << "\nsat_narrow " << counts.narrow
        << "\ngjk_rejected " << counts.rejected << "\nsingle_point_fallback " << counts.single
        << "\nsat_sweeps " << counts.sweeps << "\nsat_iterations " << counts.iterations
        << "\nsat_unresolved " << counts.exhausted << "\nsat_initial_contact " << counts.initial
        << "\nchild_candidates " << counts.candidates << "\nchild_culled " << counts.culled << '\n';
}
} // namespace rf::experimental
