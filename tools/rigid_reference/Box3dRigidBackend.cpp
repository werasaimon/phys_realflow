// Optional Box3D stepping backend: one SDK body maps to one native body, with full collision refresh.
// Input/output stay in SI; numerical length scaling scales geometry, inertia, gravity and loads together.
#include "Box3dRigidBackend.h"
#include "core/Probe.h"
#include <box3d/box3d.h>
#include <box3d/collision.h>
#include <algorithm>
#include <array>
#include <chrono>
#include <memory>
#include <set>
#include <stdexcept>
#include <unordered_map>

namespace rf::reference {
namespace {
b3Vec3 vector(const rf::Vector3& v) { return {v.x,v.y,v.z}; }
rf::Vector3 vector(b3Vec3 v) { return {v.x,v.y,v.z}; }
b3Quat quaternion(const rf::Quaternion& q) { return {{q.x,q.y,q.z},q.w}; }
}

struct Box3dRigidBackend::Impl {
    b3WorldId world = b3_nullWorldId;
    std::vector<b3BodyId> bodies;
    std::vector<b3ContactData> contacts;
    float supportError = 0, deepest = 0;
    float scale = 1;
    size_t pointCount = 0;
    int importedParts = 0;
    std::vector<RigidBody> imported;
    Timings lastTimings;
    AABB domain;
    bool walls=false;
    using Hull = std::unique_ptr<b3HullData,decltype(&b3DestroyHull)>;
    std::unordered_map<const rf::ConvexShape*,std::vector<Hull>> hulls;
    ~Impl() { if (B3_IS_NON_NULL(world)) b3DestroyWorld(world); }

    // A convex solid is the union of tetrahedra from an interior point to its boundary triangles.
    // Hulls of groups of those tetrahedra stay inside the original solid and still cover it.
    // Limit each group to 30 boundary vertices + the interior point: <=174 half-edges.
    const std::vector<Hull>& splitHull(const rf::ConvexHullShape& source) {
        const auto found = hulls.find(&source);
        if (found != hulls.end()) return found->second;
        const auto& vertices = source.vertices();
        std::vector<b3Vec3> all;
        rf::Vector3 center(0);
        for (const auto& v : vertices) { all.push_back(vector(v*scale)); center += v; }
        center /= float(vertices.size());
        std::vector<Hull> result;
        Hull whole(b3CreateHull(all.data(),int(all.size()),int(all.size())),b3DestroyHull);
        if (whole) result.push_back(std::move(whole));
        else {
            std::vector<uint32_t> group;
            auto flush = [&] {
                std::vector<b3Vec3> points{vector(center*scale)};
                for (uint32_t index : group) points.push_back(vector(vertices[index]*scale));
                Hull hull(b3CreateHull(points.data(),int(points.size()),int(points.size())),b3DestroyHull);
                if (!hull) throw std::runtime_error("Box3D tetrahedron group hull failed");
                result.push_back(std::move(hull));
                group.clear();
            };
            for (size_t t = 0; t < source.mesh()->triangles.size(); ++t) {
                const auto& triangle = source.mesh()->triangles[t];
                const auto n = source.mesh()->faceNormal(t);
                for (const auto& v : vertices)
                    if (rf::dot(n,v-vertices[triangle[0]]) > 5e-6f)
                        throw std::runtime_error("source mesh is not a convex boundary for exact partition");
                int fresh = 0;
                for (uint32_t index : triangle) fresh += std::find(group.begin(),group.end(),index) == group.end();
                if (group.size()+size_t(fresh) > 30) flush();
                for (uint32_t index : triangle)
                    if (std::find(group.begin(),group.end(),index) == group.end()) group.push_back(index);
            }
            if (!group.empty()) flush();
        }
        return hulls.emplace(&source,std::move(result)).first->second;
    }

    void validateHull(const rf::CompoundShape::Child& child, const std::vector<const b3HullData*>& pieces) {
        const auto& source = static_cast<const rf::ConvexHullShape&>(*child.shape);
        const auto& mesh = *source.mesh();
        for (size_t i = 0; i < mesh.triangles.size(); ++i) {
            const auto local = mesh.faceNormal(i);
            const auto n = child.R * local;
            float support = -rf::kInf;
            for (const auto* hull : pieces) {
                const auto* points = b3GetHullPoints(hull);
                for (int k = 0; k < hull->vertexCount; ++k) support = std::max(support,rf::dot(n,vector(points[k])/scale));
            }
            const float expected = rf::dot(n,child.t + child.R * source.support(local));
            supportError = std::max(supportError,std::fabs(support-expected));
        }
        if (supportError > 5e-6f) throw std::runtime_error("Box3D hull support import exceeds 5 micrometers");
    }

    void addBody(const rf::RigidBody& source, const rf::RigidParams& params) {
        b3BodyDef body = b3DefaultBodyDef();
        body.type = source.invMass > 0 ? b3_dynamicBody : b3_staticBody;
        body.position = vector(source.pos*scale);
        body.rotation = quaternion(source.rot);
        body.linearVelocity = vector(source.vel*scale);
        body.angularVelocity = vector(source.angVel);
        body.linearDamping = params.linearDamping;
        body.angularDamping = params.angularDamping;
        body.enableSleep = false;
        const b3BodyId id = b3CreateBody(world,&body);
        if (B3_IS_NULL(id)) throw std::runtime_error("Box3D body creation failed");
        b3ShapeDef shape = b3DefaultShapeDef();
        shape.baseMaterial.friction = source.friction;
        shape.baseMaterial.restitution = source.restitution;
        shape.updateBodyMass = false;
        addShapes(id,shape,source);
        if (source.invMass > 0) setMass(id,source);
        bodies.push_back(id);
    }

    void addShapes(b3BodyId id, const b3ShapeDef& shape, const RigidBody& source) {
        if (source.type()==ShapeType::Box) {
            const auto h=static_cast<const BoxShape&>(*source.shape).halfExtents()*scale;
            const auto hull=b3MakeBoxHull(h.x,h.y,h.z);
            if (B3_IS_NULL(b3CreateHullShape(id,&shape,&hull.base))) throw std::runtime_error("Box3D box import failed");
            ++importedParts; return;
        }
        if (source.type()==ShapeType::Sphere) {
            b3Sphere sphere{}; sphere.radius=static_cast<const SphereShape&>(*source.shape).radius()*scale;
            if (B3_IS_NULL(b3CreateSphereShape(id,&shape,&sphere))) throw std::runtime_error("Box3D sphere import failed");
            ++importedParts; return;
        }
        std::vector<CompoundShape::Child> parts;
        if (source.type()==ShapeType::Compound) parts=static_cast<const CompoundShape&>(*source.shape).children();
        else if (source.type()==ShapeType::ConvexHull) {
            CompoundShape::Child c; c.shape=std::static_pointer_cast<const ConvexHullShape>(source.shape); parts.push_back(c);
        } else throw std::runtime_error("Box3D backend: unsupported shape");
        for (const auto& child : parts) {
            const auto& pieces = splitHull(static_cast<const rf::ConvexHullShape&>(*child.shape));
            std::vector<const b3HullData*> cloned;
            const b3Transform transform{vector(child.t*scale),quaternion(rf::Quaternion::fromMatrix3x3(child.R))};
            for (const auto& hull : pieces) {
                const auto shapeId = b3CreateTransformedHullShape(id,&shape,hull.get(),transform,{1,1,1});
                if (B3_IS_NULL(shapeId)) throw std::runtime_error("Box3D shape creation failed");
                cloned.push_back(b3Shape_GetHull(shapeId));
                ++importedParts;
            }
            validateHull(child,cloned);
        }
    }

    void setMass(b3BodyId id, const RigidBody& source) {
        b3MassData mass{};
        mass.mass = source.mass;
        // With x'=s*x and unchanged mass/time: I'=s^2*I, F'=s*F, torque'=s^2*torque.
        const float s2=scale*scale;
        mass.inertia = {{s2/source.invInertiaLocal.x,0,0},{0,s2/source.invInertiaLocal.y,0},{0,0,s2/source.invInertiaLocal.z}};
        b3Body_SetMassData(id,mass);
        const auto readback = b3Body_GetMassData(id);
        if (std::fabs(readback.mass-source.mass) > 1e-5f*source.mass || rf::length(vector(readback.center)) > 1e-7f ||
            rf::length(vector(readback.inertia.cx)-vector(mass.inertia.cx)) > 1e-4f*rf::length(vector(mass.inertia.cx)) ||
            rf::length(vector(readback.inertia.cy)-vector(mass.inertia.cy)) > 1e-4f*rf::length(vector(mass.inertia.cy)) ||
            rf::length(vector(readback.inertia.cz)-vector(mass.inertia.cz)) > 1e-4f*rf::length(vector(mass.inertia.cz)))
            throw std::runtime_error("Box3D mass/inertia import mismatch");
    }

    void load(const rf::RigidWorld& source) {
        if (!source.joints().empty() || (source.staticMesh() && source.staticMesh()->triangleCount()>0)) throw std::runtime_error("Box3D backend: SDK joints/static mesh are unsupported");
        if (B3_IS_NON_NULL(world)) b3DestroyWorld(world);
        world=b3_nullWorldId; bodies.clear(); hulls.clear(); contacts.clear();
        supportError=deepest=0; importedParts=0; pointCount=0;
        imported=source.bodies();
        domain=source.domain(); walls=source.params.collideWithDomain;
        b3WorldDef settings = b3DefaultWorldDef();
        settings.gravity = vector(source.params.gravity*scale);
        settings.restitutionThreshold*=scale; settings.hitEventThreshold*=scale;
        settings.contactSpeed*=scale; settings.maximumLinearSpeed*=scale;
        settings.enableSleep = false;
        settings.enableContinuous = true;
        settings.workerCount = 1;
        world = b3CreateWorld(&settings);
        if (B3_IS_NULL(world)) throw std::runtime_error("Box3D world creation failed");
        for (const auto& body : source.bodies()) addBody(body,source.params);
        if (!source.params.collideWithDomain) return;
        const auto box = source.domain();
        const auto center = (box.lo + box.hi)*0.5f;
        const auto half = (box.hi - box.lo)*0.5f;
        for (int axis = 0; axis < 3; ++axis) for (int side : {-1,1}) {
            auto position = center;
            auto extent = half + rf::Vector3(1);
            position[axis] += float(side)*(half[axis]+0.5f);
            extent[axis] = 0.5f;
            auto body = b3DefaultBodyDef();
            body.position = vector(position*scale);
            const auto id = b3CreateBody(world,&body);
            auto shape = b3DefaultShapeDef();
            shape.baseMaterial.friction = 0.6f;
            shape.baseMaterial.restitution = 0;
            const auto hull = b3MakeBoxHull(extent.x*scale,extent.y*scale,extent.z*scale);
            if (B3_IS_NULL(b3CreateHullShape(id,&shape,&hull.base))) throw std::runtime_error("Box3D wall creation failed");
        }
    }

    void observe(rf::RigidWorld& target) {
        pointCount=0;
        std::set<std::array<uint32_t,3>> seen;
        for (size_t i = 0; i < bodies.size(); ++i) {
            auto& body = target.bodies()[i];
            body.pos = vector(b3Body_GetPosition(bodies[i]))/scale;
            const auto q = b3Body_GetRotation(bodies[i]);
            body.rot = {q.s,q.v.x,q.v.y,q.v.z};
            body.vel = vector(b3Body_GetLinearVelocity(bodies[i]))/scale;
            body.angVel = vector(b3Body_GetAngularVelocity(bodies[i]));
            body.sleeping=false; body.sleepTimer=0; body.sleepIsland=-1;
            body.updateInertia();
            if (!std::isfinite(rf::length2(body.pos)) || !std::isfinite(rf::length2(body.vel)) ||
                !std::isfinite(rf::length2(body.angVel)) || !std::isfinite(body.rot.w) ||
                !std::isfinite(length2(Vector3(body.rot.x,body.rot.y,body.rot.z))))
                throw std::runtime_error("non-finite Box3D state");
            contacts.resize(size_t(b3Body_GetContactCapacity(bodies[i])));
            const int count = b3Body_GetContactData(bodies[i],contacts.data(),int(contacts.size()));
            for (int k = 0; k < count; ++k) {
                std::array<uint32_t,3> key; b3StoreContactId(contacts[k].contactId,key.data());
                if (!seen.insert(key).second) continue;
                for (int m = 0; m < contacts[k].manifoldCount; ++m) {
                    const auto& manifold = contacts[k].manifolds[m];
                    pointCount+=size_t(manifold.pointCount);
                    for (int p = 0; p < manifold.pointCount; ++p) deepest = std::max(deepest,-manifold.points[p].separation/scale);
                }
            }
        }
    }

    void input(RigidWorld& target) {
        if (target.bodies().size()!=bodies.size()) throw std::runtime_error("reload backend after changing the body list");
        if (!target.joints().empty() || (target.staticMesh() && target.staticMesh()->triangleCount()>0)) throw std::runtime_error("Box3D backend: SDK joints/static mesh are unsupported");
        if (walls!=target.params.collideWithDomain || length2(domain.lo-target.domain().lo)>0 || length2(domain.hi-target.domain().hi)>0)
            throw std::runtime_error("reload backend after changing domain walls");
        // Validate all immutable inputs before applying any loads to the native world.
        for (size_t i=0; i<bodies.size(); ++i) {
            const auto& b=target.bodies()[i]; const auto& old=imported[i];
            if (b.shape!=old.shape || b.mass!=old.mass || b.invMass!=old.invMass || b.alive!=old.alive ||
                length2(b.invInertiaLocal-old.invInertiaLocal)>0 || b.friction!=old.friction || b.restitution!=old.restitution)
                throw std::runtime_error("reload backend after changing shapes, mass, inertia, type or material");
        }
        b3World_SetGravity(world,vector(target.params.gravity*scale));
        for (size_t i=0; i<bodies.size(); ++i) {
            auto& body=target.bodies()[i];
            const auto p=vector(b3Body_GetPosition(bodies[i]))/scale;
            const auto q=b3Body_GetRotation(bodies[i]);
            if (length2(p-body.pos)>1e-14f || length2(Vector3(q.v.x,q.v.y,q.v.z)-Vector3(body.rot.x,body.rot.y,body.rot.z))>1e-14f || q.s!=body.rot.w)
                b3Body_SetTransform(bodies[i],vector(body.pos*scale),quaternion(body.rot));
            b3Body_SetLinearVelocity(bodies[i],vector(body.vel*scale));
            b3Body_SetAngularVelocity(bodies[i],vector(body.angVel));
            b3Body_SetLinearDamping(bodies[i],target.params.linearDamping);
            b3Body_SetAngularDamping(bodies[i],target.params.angularDamping);
            b3Body_ApplyForceToCenter(bodies[i],vector(body.force*scale),true);
            b3Body_ApplyTorque(bodies[i],vector(body.torque*(scale*scale)),true);
            body.prevPos=body.pos; body.prevRot=body.rot;
            body.force=body.torque=body.biasVel=body.biasAngVel=Vector3(0);
        }
    }
};

Box3dRigidBackend::Box3dRigidBackend(float scale) : impl_(std::make_unique<Impl>()) {
    if (!std::isfinite(scale) || scale<1 || scale>100) throw std::runtime_error("numerical length scale must be in [1,100]");
    impl_->scale=scale;
}
Box3dRigidBackend::~Box3dRigidBackend()=default;
void Box3dRigidBackend::load(RigidWorld& source) { impl_->load(source); impl_->observe(source); }
void Box3dRigidBackend::step(RigidWorld& target,float dt,int substeps) {
    if (!B3_IS_NON_NULL(impl_->world) || !std::isfinite(dt) || dt<=0 || substeps<1) throw std::runtime_error("invalid backend step");
    using Clock=std::chrono::steady_clock;
    const auto start=Clock::now(); impl_->input(target); const auto ready=Clock::now();
    b3World_Step(impl_->world,dt,substeps);
    const auto done=Clock::now(); impl_->observe(target); const auto exported=Clock::now();
    impl_->lastTimings={std::chrono::duration<double,std::milli>(ready-start).count(),
        std::chrono::duration<double,std::milli>(done-ready).count(),std::chrono::duration<double,std::milli>(exported-done).count()};
    Probe::set("box3d/contact points",double(impl_->pointCount));
}
float Box3dRigidBackend::deepestContact() const { return impl_->deepest; }
float Box3dRigidBackend::supportImportError() const { return impl_->supportError; }
int Box3dRigidBackend::importedPartCount() const { return impl_->importedParts; }
size_t Box3dRigidBackend::contactPointCount() const { return impl_->pointCount; }
float Box3dRigidBackend::numericalLengthScale() const { return impl_->scale; }
Box3dRigidBackend::Timings Box3dRigidBackend::timings() const { return impl_->lastTimings; }
Box3dRigidBackend::Profile Box3dRigidBackend::profile() const {
    auto p=b3World_GetProfile(impl_->world); return {p.pairs,p.collide,p.solve,p.bullets};
}


} // namespace rf::reference
