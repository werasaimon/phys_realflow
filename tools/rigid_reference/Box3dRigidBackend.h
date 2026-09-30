// Optional native Box3D backend for SDK rigid scenes. The main SDK keeps no engine dependency.
// A step imports external velocities/loads, executes native collision/response and exports SI poses.
// Supported: boxes, spheres, hulls and compounds, domain walls; no SDK joints or static triangle mesh.
#pragma once
#include "rigid/RigidWorld.h"
#include <memory>

namespace rf::reference {
class Box3dRigidBackend {
public:
    struct Profile { double pairs=0,collide=0,solve=0,bullets=0; };
    struct Timings { double input=0,engine=0,output=0; };
    explicit Box3dRigidBackend(float numericalLengthScale=1);
    ~Box3dRigidBackend();
    Box3dRigidBackend(const Box3dRigidBackend&)=delete;
    Box3dRigidBackend& operator=(const Box3dRigidBackend&)=delete;
    void load(RigidWorld& source);
    // substeps=1 refreshes native collisions on every call. For chains call once per SDK substep.
    void step(RigidWorld& target,float dt,int substeps=1);
    float deepestContact() const;
    float supportImportError() const;
    int importedPartCount() const;
    size_t contactPointCount() const;
    Profile profile() const;
    Timings timings() const;
    float numericalLengthScale() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
} // namespace rf::reference
