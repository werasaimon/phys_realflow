#include "rigid/Shapes.h"

#include <algorithm>

namespace rf {

AABB orientedBounds(const AABB& local, const Matrix3x3& R, const Vector3& p) {
    Vector3 c = R * local.center() + p, h = local.extent() * 0.5f, e;
    for (int i = 0; i < 3; ++i) e[i] = std::fabs(R.m[i][0]) * h.x + std::fabs(R.m[i][1]) * h.y + std::fabs(R.m[i][2]) * h.z;
    return AABB(c - e, c + e);
}

AABB ConvexShape::boundsAt(const Matrix3x3& R, const Vector3& p) const {
    if (type() == ShapeType::Sphere) {
        float r = boundingRadius();
        return AABB(p - Vector3(r), p + Vector3(r));
    }
    if (type() == ShapeType::Box) return orientedBounds(localBounds(), R, p);
    AABB b;
    const Matrix3x3 Rt = R.transposed();
    for (int a = 0; a < 3; ++a) {
        Vector3 e(0.0f);
        e[a] = 1;
        b.hi[a] = (R * support(Rt * e))[a] + p[a];
        b.lo[a] = (R * support(Rt * -e))[a] + p[a];
    }
    return b;
}

// ---------------------------------------------------------------------------
// Supporting features
// ---------------------------------------------------------------------------
static void orderAround(const Vector3& normal, std::vector<Vector3>& pts) {
    if (pts.size() < 3) return;
    Vector3 c(0.0f);
    for (const Vector3& p : pts) c += p;
    c /= float(pts.size());
    Vector3 n = normalize(normal);
    Vector3 u = normalize(std::fabs(n.x) > 0.57f ? Vector3(n.y, -n.x, 0) : Vector3(0, n.z, -n.y));
    Vector3 v = cross(n, u);
    std::sort(pts.begin(), pts.end(), [&](const Vector3& a, const Vector3& b) {
        return std::atan2(dot(a - c, v), dot(a - c, u)) < std::atan2(dot(b - c, v), dot(b - c, u));
    });
}

std::vector<ConvexShape::Face> buildFaces(const TriMesh& m) {
    std::vector<ConvexShape::Face> faces;
    std::vector<float> offsets;
    float size = maxComp(m.bounds().extent());
    for (size_t t = 0; t < m.triangles.size(); ++t) {
        Vector3 n = m.faceNormal(t);
        if (length2(n) < 0.5f) continue;
        float d = dot(n, m.positions[m.triangles[t][0]]);
        size_t f = 0;
        for (; f < faces.size(); ++f)
            if (dot(faces[f].normal, n) > 0.9999f && std::fabs(offsets[f] - d) < 1e-4f * size) break;
        if (f == faces.size()) {
            faces.push_back({n, {}});
            offsets.push_back(d);
        }
        for (int k = 0; k < 3; ++k) {
            const Vector3& p = m.positions[m.triangles[t][k]];
            bool dup = false;
            for (const Vector3& q : faces[f].verts) dup |= length2(q - p) < 1e-12f;
            if (!dup) faces[f].verts.push_back(p);
        }
    }
    for (auto& f : faces) orderAround(f.normal, f.verts);
    return faces;
}

void ConvexShape::supportFeature(const Vector3& dir, std::vector<Vector3>& out) const {
    const std::vector<Face>* fs = faces();
    if (!fs || fs->empty()) {
        out.assign(1, support(dir));
        return;
    }
    const Face* best = &(*fs)[0];
    float bd = -kInf;
    for (const Face& f : *fs) {
        float d = dot(f.normal, dir);
        if (d > bd) { bd = d; best = &f; }
    }
    out = best->verts;
}

void TriangleShape::supportFeature(const Vector3&, std::vector<Vector3>& out) const {
    out.assign(v_, v_ + 3); // a triangle is its own (two-sided) face
}

// ---------------------------------------------------------------------------
// Ray casts (local frame)
// ---------------------------------------------------------------------------
bool SphereShape::raycast(const Vector3& o, const Vector3& d, float maxT, float& t, Vector3& n) const {
    float b = dot(o, d), c = dot(o, o) - r_ * r_;
    if (c > 0 && b > 0) return false;
    float disc = b * b - c;
    if (disc < 0) return false;
    t = std::max(-b - std::sqrt(disc), 0.0f);
    if (t > maxT) return false;
    n = normalize(o + d * t);
    return true;
}

bool BoxShape::raycast(const Vector3& o, const Vector3& d, float maxT, float& t, Vector3& n) const {
    float t0 = 0, t1 = maxT;
    int axis = -1;
    float sign = 1;
    for (int a = 0; a < 3; ++a) {
        if (std::fabs(d[a]) < 1e-12f) {
            if (o[a] < -h_[a] || o[a] > h_[a]) return false;
            continue;
        }
        float inv = 1.0f / d[a];
        float tn = (-h_[a] - o[a]) * inv, tf = (h_[a] - o[a]) * inv;
        float s = -1;
        if (tn > tf) { std::swap(tn, tf); s = 1; }
        if (tn > t0) { t0 = tn; axis = a; sign = s; }
        t1 = std::min(t1, tf);
        if (t0 > t1) return false;
    }
    t = t0;
    n = Vector3(0.0f);
    if (axis >= 0) n[axis] = sign;
    else n = normalize(-d); // origin inside
    return true;
}

// Cyrus-Beck clipping of the ray against the face planes.
bool ConvexHullShape::raycast(const Vector3& o, const Vector3& d, float maxT, float& t, Vector3& n) const {
    float tEnter = 0, tExit = maxT;
    Vector3 nEnter = normalize(-d);
    for (const Plane& pl : planes_) {
        float denom = dot(pl.n, d);
        float dist = pl.d - dot(pl.n, o); // >= 0 when the origin is inside this plane
        if (std::fabs(denom) < 1e-12f) {
            if (dist < 0) return false;
            continue;
        }
        float tt = dist / denom;
        if (denom < 0) {
            if (tt > tEnter) { tEnter = tt; nEnter = pl.n; }
        } else {
            tExit = std::min(tExit, tt);
        }
        if (tEnter > tExit) return false;
    }
    t = tEnter;
    n = nEnter;
    return true;
}

// ---------------------------------------------------------------------------
// Compound
// ---------------------------------------------------------------------------
CompoundShape::CompoundShape(const std::vector<TriMesh>& hulls, const TriMesh& visual) {
    struct Tmp { std::shared_ptr<ConvexHullShape> s; Matrix3x3 R; Vector3 t; float m; };
    std::vector<Tmp> tmp;
    Vector3 com(0.0f);
    float V = 0;
    for (const TriMesh& h : hulls) {
        auto s = std::make_shared<ConvexHullShape>(h);
        if (s->volume() <= 1e-12f) continue;
        tmp.push_back({s, s->principalRotation(), s->centerOfMass(), s->volume()});
        com += s->centerOfMass() * s->volume();
        V += s->volume();
    }
    volume_ = std::max(V, 1e-12f);
    com /= volume_;
    // Inertia of the union (unit density) about its centre of mass: rotated part tensors plus
    // parallel-axis terms, then diagonalised -> principal frame of the compound.
    Matrix3x3 I = Matrix3x3::zero();
    for (const Tmp& c : tmp) {
        Vector3 Ip = c.s->unitInertia() * c.m;
        Matrix3x3 Ic = c.R * Matrix3x3::diag(Ip) * c.R.transposed();
        Vector3 d = c.t - com;
        float d2 = dot(d, d);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) I.m[i][j] += Ic.m[i][j] + c.m * ((i == j ? d2 : 0.0f) - d[i] * d[j]);
    }
    Vector3 eig;
    Matrix3x3 P;
    symmetricEigen(I, eig, P);
    if (dot(cross(P.col(0), P.col(1)), P.col(2)) < 0)
        for (int i = 0; i < 3; ++i) P.m[i][2] = -P.m[i][2];
    inertia_ = Vector3(std::max(eig.x, 1e-12f), std::max(eig.y, 1e-12f), std::max(eig.z, 1e-12f)) / volume_;
    com_ = com;
    principal_ = P;
    Matrix3x3 Pt = P.transposed();
    for (const Tmp& c : tmp) {
        Child ch{c.s, Pt * c.R, Pt * (c.t - com), c.m};
        children_.push_back(ch);
        TriMesh pm = *c.s->mesh();
        for (Vector3& v : pm.positions) v = ch.R * v + ch.t;
        partMeshes_.push_back(std::make_shared<TriMesh>(std::move(pm)));
    }
    TriMesh vis = visual;
    for (Vector3& v : vis.positions) v = Pt * (v - com);
    visual_ = std::make_shared<TriMesh>(std::move(vis));
    radius_ = 0;
    for (const Child& c : children_) {
        radius_ = std::max(radius_, length(c.t) + c.shape->boundingRadius());
        for (const Vector3& v : c.shape->vertices()) bounds_.expand(c.R * v + c.t);
    }
}

AABB CompoundShape::boundsAt(const Matrix3x3& R, const Vector3& p) const {
    AABB b;
    for (const Child& c : children_) b.expand(orientedBounds(c.shape->localBounds(), R * c.R, p + R * c.t));
    return b;
}

Vector3 CompoundShape::support(const Vector3& d) const {
    Vector3 best;
    float bd = -kInf;
    for (const Child& c : children_) {
        Vector3 p = c.R * c.shape->support(c.R.transposed() * d) + c.t;
        float v = dot(p, d);
        if (v > bd) { bd = v; best = p; }
    }
    return best;
}

bool CompoundShape::contains(const Vector3& p) const {
    for (const Child& c : children_) {
        Vector3 q = c.R.transposed() * (p - c.t);
        if (c.shape->localBounds().contains(q) && c.shape->contains(q)) return true;
    }
    return false;
}

float CompoundShape::signedDistance(const Vector3& p, Vector3& n) const {
    float best = kInf;
    for (const Child& c : children_) {
        Vector3 nl;
        float d = c.shape->signedDistance(c.R.transposed() * (p - c.t), nl);
        if (d < best) { best = d; n = c.R * nl; }
    }
    return best;
}

bool CompoundShape::raycast(const Vector3& o, const Vector3& d, float maxT, float& t, Vector3& normal) const {
    bool hit = false;
    for (const Child& c : children_) {
        Matrix3x3 Rt = c.R.transposed();
        float tc;
        Vector3 nc;
        if (c.shape->raycast(Rt * (o - c.t), Rt * d, maxT, tc, nc) && tc < maxT) {
            maxT = tc;
            t = tc;
            normal = c.R * nc;
            hit = true;
        }
    }
    return hit;
}

float SphereShape::signedDistance(const Vector3& p, Vector3& n) const {
    float l = length(p);
    n = l > 1e-9f ? p / l : Vector3(0, 1, 0);
    return l - r_;
}

float BoxShape::signedDistance(const Vector3& q, Vector3& n) const {
    Vector3 d = vabs(q) - h_;
    Vector3 outside = vmax(d, Vector3(0.0f));
    float outLen = length(outside);
    if (outLen > 0) {
        n = Vector3(q.x < 0 ? -outside.x : outside.x, q.y < 0 ? -outside.y : outside.y, q.z < 0 ? -outside.z : outside.z) / outLen;
        return outLen;
    }
    int ax = (d.x > d.y) ? (d.x > d.z ? 0 : 2) : (d.y > d.z ? 1 : 2);
    n = Vector3(0.0f);
    n[ax] = q[ax] < 0 ? -1.0f : 1.0f;
    return d[ax];
}

// ---------------------------------------------------------------------------
// Mass properties (D. Eberly, "Polyhedral Mass Properties (Revisited)")
// ---------------------------------------------------------------------------
static void subexpressions(double w0, double w1, double w2, double& f1, double& f2, double& f3, double& g0,
                           double& g1, double& g2) {
    double temp0 = w0 + w1;
    f1 = temp0 + w2;
    double temp1 = w0 * w0;
    double temp2 = temp1 + w1 * temp0;
    f2 = temp2 + w2 * f1;
    f3 = w0 * temp1 + w1 * temp2 + w2 * f2;
    g0 = f2 + w0 * (f1 + w0);
    g1 = f2 + w1 * (f1 + w1);
    g2 = f2 + w2 * (f1 + w2);
}

MassProperties computeMassProperties(const TriMesh& m) {
    const double mult[10] = {1 / 6.0, 1 / 24.0, 1 / 24.0, 1 / 24.0, 1 / 60.0, 1 / 60.0, 1 / 60.0, 1 / 120.0, 1 / 120.0, 1 / 120.0};
    double in[10] = {0};
    for (const auto& t : m.triangles) {
        const Vector3 &p0 = m.positions[t[0]], &p1 = m.positions[t[1]], &p2 = m.positions[t[2]];
        double x0 = p0.x, y0 = p0.y, z0 = p0.z, x1 = p1.x, y1 = p1.y, z1 = p1.z, x2 = p2.x, y2 = p2.y, z2 = p2.z;
        double a1 = x1 - x0, b1 = y1 - y0, c1 = z1 - z0, a2 = x2 - x0, b2 = y2 - y0, c2 = z2 - z0;
        double d0 = b1 * c2 - b2 * c1, d1 = a2 * c1 - a1 * c2, d2 = a1 * b2 - a2 * b1;
        double f1x, f2x, f3x, g0x, g1x, g2x, f1y, f2y, f3y, g0y, g1y, g2y, f1z, f2z, f3z, g0z, g1z, g2z;
        subexpressions(x0, x1, x2, f1x, f2x, f3x, g0x, g1x, g2x);
        subexpressions(y0, y1, y2, f1y, f2y, f3y, g0y, g1y, g2y);
        subexpressions(z0, z1, z2, f1z, f2z, f3z, g0z, g1z, g2z);
        in[0] += d0 * f1x;
        in[1] += d0 * f2x;
        in[2] += d1 * f2y;
        in[3] += d2 * f2z;
        in[4] += d0 * f3x;
        in[5] += d1 * f3y;
        in[6] += d2 * f3z;
        in[7] += d0 * (y0 * g0x + y1 * g1x + y2 * g2x);
        in[8] += d1 * (z0 * g0y + z1 * g1y + z2 * g2y);
        in[9] += d2 * (x0 * g0z + x1 * g1z + x2 * g2z);
    }
    for (int i = 0; i < 10; ++i) in[i] *= mult[i];
    MassProperties mp;
    double mass = in[0];
    mp.volume = float(mass);
    if (mass <= 0) return mp;
    double cx = in[1] / mass, cy = in[2] / mass, cz = in[3] / mass;
    mp.com = Vector3(float(cx), float(cy), float(cz));
    double Ixx = in[5] + in[6] - mass * (cy * cy + cz * cz);
    double Iyy = in[4] + in[6] - mass * (cz * cz + cx * cx);
    double Izz = in[4] + in[5] - mass * (cx * cx + cy * cy);
    double Ixy = -(in[7] - mass * cx * cy);
    double Iyz = -(in[8] - mass * cy * cz);
    double Ixz = -(in[9] - mass * cz * cx);
    Matrix3x3& I = mp.inertia;
    I.m[0][0] = float(Ixx); I.m[0][1] = float(Ixy); I.m[0][2] = float(Ixz);
    I.m[1][0] = float(Ixy); I.m[1][1] = float(Iyy); I.m[1][2] = float(Iyz);
    I.m[2][0] = float(Ixz); I.m[2][1] = float(Iyz); I.m[2][2] = float(Izz);
    return mp;
}

// ---------------------------------------------------------------------------
// Convex hull
// ---------------------------------------------------------------------------
ConvexHullShape::ConvexHullShape(const TriMesh& input) {
    TriMesh m = input;
    m.orientOutward();
    MassProperties mp = computeMassProperties(m);
    volume_ = std::max(mp.volume, 1e-12f);
    com_ = mp.com;
    Vector3 eig;
    symmetricEigen(mp.inertia, eig, principal_);
    // Proper rotation (det = +1).
    Vector3 c0 = principal_.col(0), c1 = principal_.col(1), c2 = principal_.col(2);
    if (dot(cross(c0, c1), c2) < 0)
        for (int i = 0; i < 3; ++i) principal_.m[i][2] = -principal_.m[i][2];
    inertia_ = Vector3(std::max(eig.x, 1e-12f), std::max(eig.y, 1e-12f), std::max(eig.z, 1e-12f)) / volume_;

    Matrix3x3 Rt = principal_.transposed();
    for (Vector3& p : m.positions) p = Rt * (p - com_);
    bounds_ = m.bounds();
    radius_ = 0;
    for (const Vector3& p : m.positions) radius_ = std::max(radius_, length(p));
    for (size_t t = 0; t < m.triangles.size(); ++t) {
        Vector3 n = m.faceNormal(t);
        if (length2(n) < 0.5f) continue;
        planes_.push_back({n, dot(n, m.positions[m.triangles[t][0]])});
    }
    faces_ = buildFaces(m);
    mesh_ = std::make_shared<TriMesh>(std::move(m));
}

Vector3 ConvexHullShape::support(const Vector3& d) const {
    const auto& P = mesh_->positions;
    size_t best = 0;
    float bestDot = -kInf;
    for (size_t i = 0; i < P.size(); ++i) {
        float v = dot(P[i], d);
        if (v > bestDot) { bestDot = v; best = i; }
    }
    return P[best];
}

bool ConvexHullShape::contains(const Vector3& p) const {
    for (const Plane& pl : planes_)
        if (dot(pl.n, p) > pl.d) return false; // outside one face plane: outside the hull
    return true;
}

float ConvexHullShape::signedDistance(const Vector3& p, Vector3& n) const {
    // Max over face planes: exact inside, a lower bound of the distance outside.
    float best = -kInf;
    for (const Plane& pl : planes_) {
        float d = dot(pl.n, p) - pl.d;
        if (d > best) { best = d; n = pl.n; }
    }
    return best;
}

} // namespace rf
