#pragma once
// Convex collision shapes. Every shape is described by its support mapping
//     support(d) = argmax_{x in shape} dot(x, d)       (local frame, centre of mass at the origin)
// which is all GJK / EPA need; shapes also provide mass properties and an approximate signed
// distance (used by the SPH <-> rigid coupling).

#include "math/Math.h"
#include "core/Mesh.h"

#include <memory>
#include <vector>

namespace rf {

enum class ShapeType { Sphere, Box, ConvexHull, Triangle, Compound };

class ConvexShape {
public:
    virtual ~ConvexShape() = default;
    virtual ShapeType type() const = 0;
    virtual Vector3 support(const Vector3& dir) const = 0;
    virtual AABB localBounds() const = 0;
    virtual float volume() const = 0;
    // Principal moments of inertia per unit mass (local axes are principal axes).
    virtual Vector3 unitInertia() const = 0;
    // Signed distance from a local point to the surface (negative inside) and outward normal.
    virtual float signedDistance(const Vector3& p, Vector3& normal) const = 0;
    // Is the local point inside? (Cheaper than signedDistance where shapes can exit early.)
    virtual bool contains(const Vector3& p) const {
        Vector3 n;
        return signedDistance(p, n) < 0;
    }
    virtual float boundingRadius() const = 0;
    // World AABB of the shape posed at (R, p). Default: exact, from the support mapping along the
    // world axes; shapes override it with something cheaper.
    virtual AABB boundsAt(const Matrix3x3& R, const Vector3& p) const;

    // Supporting face ("boundary simplex") in direction `dir`, as in Jolt's GetSupportingFace: the
    // polygon of the face whose outward normal is most aligned with dir (ordered counter-clockwise
    // seen from outside). Shapes without faces return their single support point.
    virtual void supportFeature(const Vector3& dir, std::vector<Vector3>& out) const;

    // Ray in local coordinates (origin o, unit direction d): first entry point with t in [0, maxT].
    virtual bool raycast(const Vector3& o, const Vector3& d, float maxT, float& t, Vector3& normal) const {
        (void)o; (void)d; (void)maxT; (void)t; (void)normal;
        return false;
    }

    struct Face {
        Vector3 normal;              // outward, unit
        std::vector<Vector3> verts;  // convex polygon, CCW around normal
    };

protected:
    virtual const std::vector<Face>* faces() const { return nullptr; }
};

// AABB of the oriented box (local AABB `local` posed at R, p): centre R c + p, half extents |R| h.
AABB orientedBounds(const AABB& local, const Matrix3x3& R, const Vector3& p);

// Groups coplanar triangles of a convex mesh into polygon faces.
std::vector<ConvexShape::Face> buildFaces(const TriMesh& convexMesh);

class SphereShape final : public ConvexShape {
public:
    explicit SphereShape(float r) : r_(r) {}
    ShapeType type() const override { return ShapeType::Sphere; }
    Vector3 support(const Vector3& d) const override { return normalize(d) * r_; }
    AABB localBounds() const override { return AABB(Vector3(-r_), Vector3(r_)); }
    float volume() const override { return 4.0f / 3.0f * kPi * r_ * r_ * r_; }
    Vector3 unitInertia() const override { return Vector3(0.4f * r_ * r_); }
    float signedDistance(const Vector3& p, Vector3& n) const override;
    float boundingRadius() const override { return r_; }
    float radius() const { return r_; }
    bool raycast(const Vector3& o, const Vector3& d, float maxT, float& t, Vector3& normal) const override;

private:
    float r_;
};

class BoxShape final : public ConvexShape {
public:
    explicit BoxShape(const Vector3& h) : h_(h), faces_(buildFaces(primitives::box(h))) {}
    ShapeType type() const override { return ShapeType::Box; }
    Vector3 support(const Vector3& d) const override {
        return {d.x >= 0 ? h_.x : -h_.x, d.y >= 0 ? h_.y : -h_.y, d.z >= 0 ? h_.z : -h_.z};
    }
    AABB localBounds() const override { return AABB(-h_, h_); }
    float volume() const override { return 8.0f * h_.x * h_.y * h_.z; }
    Vector3 unitInertia() const override {
        Vector3 s = h_ * 2.0f;
        return Vector3(s.y * s.y + s.z * s.z, s.x * s.x + s.z * s.z, s.x * s.x + s.y * s.y) / 12.0f;
    }
    float signedDistance(const Vector3& p, Vector3& n) const override;
    float boundingRadius() const override { return length(h_); }
    const Vector3& halfExtents() const { return h_; }
    bool raycast(const Vector3& o, const Vector3& d, float maxT, float& t, Vector3& normal) const override;

protected:
    const std::vector<Face>* faces() const override { return &faces_; }

private:
    Vector3 h_;
    std::vector<Face> faces_;
};

// Convex polyhedron given by a closed, convex, outward-oriented triangle mesh. The constructor
// moves the centre of mass to the origin and rotates the vertices into the principal axes of
// inertia (so the body frame diagonalises the inertia tensor). No hull computation is done:
// the input must already be convex (all procedural primitives are).
class ConvexHullShape final : public ConvexShape {
public:
    explicit ConvexHullShape(const TriMesh& convexMesh);
    ShapeType type() const override { return ShapeType::ConvexHull; }
    Vector3 support(const Vector3& d) const override;
    AABB localBounds() const override { return bounds_; }
    float volume() const override { return volume_; }
    Vector3 unitInertia() const override { return inertia_; }
    float signedDistance(const Vector3& p, Vector3& n) const override;
    bool contains(const Vector3& p) const override;
    float boundingRadius() const override { return radius_; }
    const std::vector<Vector3>& vertices() const { return mesh_->positions; }
    const std::shared_ptr<const TriMesh>& mesh() const { return mesh_; }
    // Rotation applied to the input to reach the principal frame (world = R * local + offset).
    const Matrix3x3& principalRotation() const { return principal_; }
    const Vector3& centerOfMass() const { return com_; }
    bool raycast(const Vector3& o, const Vector3& d, float maxT, float& t, Vector3& normal) const override;

protected:
    const std::vector<Face>* faces() const override { return &faces_; }

private:
    std::vector<Face> faces_;
    struct Plane { Vector3 n; float d; };
    std::shared_ptr<const TriMesh> mesh_;
    std::vector<Plane> planes_;
    AABB bounds_;
    float volume_ = 0, radius_ = 0;
    Vector3 inertia_, com_;
    Matrix3x3 principal_;
};

// World-space triangle (static mesh collisions); zero volume, used only by GJK/EPA.
class TriangleShape final : public ConvexShape {
public:
    TriangleShape(const Vector3& a, const Vector3& b, const Vector3& c) : v_{a, b, c} {}
    ShapeType type() const override { return ShapeType::Triangle; }
    Vector3 support(const Vector3& d) const override {
        float d0 = dot(v_[0], d), d1 = dot(v_[1], d), d2 = dot(v_[2], d);
        return d0 >= d1 ? (d0 >= d2 ? v_[0] : v_[2]) : (d1 >= d2 ? v_[1] : v_[2]);
    }
    AABB localBounds() const override { AABB b; for (auto& v : v_) b.expand(v); return b; }
    float volume() const override { return 0; }
    Vector3 unitInertia() const override { return Vector3(0.0f); }
    float signedDistance(const Vector3&, Vector3& n) const override { n = normal(); return kInf; }
    float boundingRadius() const override { return 0; }
    Vector3 normal() const { return normalize(cross(v_[1] - v_[0], v_[2] - v_[0])); }
    const Vector3& vertex(int i) const { return v_[i]; }
    void supportFeature(const Vector3& dir, std::vector<Vector3>& out) const override;

private:
    Vector3 v_[3];
};

// Non-convex body as a set of convex hull parts (from a convex decomposition), each with its
// own pose in the compound frame (centre of mass at the origin, principal axes). Collisions go
// part against part in the narrow phase; support() is that of the union's hull (bounds, CCD).
class CompoundShape final : public ConvexShape {
public:
    struct Child {
        std::shared_ptr<const ConvexHullShape> shape;
        Matrix3x3 R;   // part frame -> compound frame
        Vector3 t;
        float mass = 0;
    };
    // `hulls`: convex parts in model coordinates; `visual`: the render mesh in model coordinates.
    CompoundShape(const std::vector<TriMesh>& hulls, const TriMesh& visual);
    ShapeType type() const override { return ShapeType::Compound; }
    Vector3 support(const Vector3& d) const override;
    AABB localBounds() const override { return bounds_; }
    float volume() const override { return volume_; }
    Vector3 unitInertia() const override { return inertia_; }
    float signedDistance(const Vector3& p, Vector3& n) const override;
    bool contains(const Vector3& p) const override;
    float boundingRadius() const override { return radius_; }
    // Union of the parts' oriented boxes: O(parts) instead of 6 support queries over all vertices.
    AABB boundsAt(const Matrix3x3& R, const Vector3& p) const override;
    bool raycast(const Vector3& o, const Vector3& d, float maxT, float& t, Vector3& normal) const override;
    const std::vector<Child>& children() const { return children_; }
    const std::shared_ptr<const TriMesh>& visualMesh() const { return visual_; }
    // Hull meshes of the parts in the compound frame (for debug drawing).
    const std::vector<std::shared_ptr<const TriMesh>>& partMeshes() const { return partMeshes_; }
    // Model -> compound frame: x_compound = principal^T (x_model - com).
    const Vector3& centerOfMass() const { return com_; }
    const Matrix3x3& principalRotation() const { return principal_; }

private:
    std::vector<Child> children_;
    std::shared_ptr<const TriMesh> visual_;
    std::vector<std::shared_ptr<const TriMesh>> partMeshes_;
    Vector3 com_;
    Matrix3x3 principal_;
    AABB bounds_;
    float volume_ = 0, radius_ = 0;
    Vector3 inertia_;
};

// Mass properties of a closed triangle mesh with unit density (Eberly, "Polyhedral Mass Properties").
struct MassProperties {
    float volume = 0;
    Vector3 com;
    Matrix3x3 inertia; // about the centre of mass, unit density
};
MassProperties computeMassProperties(const TriMesh& m);

} // namespace rf
