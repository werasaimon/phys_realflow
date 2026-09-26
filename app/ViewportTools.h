#pragma once
// Mouse interaction strategies for the Viewport (Strategy pattern): the viewport forwards mouse
// events to the active tool, so new tools (probe, brush, picking...) plug in without touching it.

#include "OrbitCamera.h"

#include <QElapsedTimer>
#include <QPoint>

class QMouseEvent;
class Viewport;

class ViewportTool {
public:
    virtual ~ViewportTool() = default;
    virtual void press(Viewport&, QMouseEvent*) {}
    virtual void move(Viewport&, QMouseEvent*) {}
    virtual void release(Viewport&, QMouseEvent*) {}
};

// Middle button (or LMB where no brush applies) orbits; Shift+middle or right button pans.
class OrbitTool : public ViewportTool {
public:
    void press(Viewport& v, QMouseEvent* e) override;
    void move(Viewport& v, QMouseEvent* e) override;

private:
    QPoint last_;
};

// LMB grabs a rigid body with a mouse joint: the cursor ray picks the body (any convex shape),
// the hit point becomes the anchor, and dragging moves the target on the plane through the anchor
// facing the camera.
class GrabTool : public ViewportTool {
public:
    void press(Viewport& v, QMouseEvent* e) override;
    void move(Viewport& v, QMouseEvent* e) override;
    void release(Viewport& v, QMouseEvent* e) override;

private:
    bool grabbing_ = false;
    rf::Vector3 planePoint_, planeNormal_;
    OrbitTool fallback_; // empty space / other buttons move the camera
};

// LMB drag injects a disturbance into the gas: the cursor ray is unprojected into 3D, intersected
// with the slice plane (or the middle of the domain) and the hit point's motion gives the velocity.
class DisturbTool : public ViewportTool {
public:
    void press(Viewport& v, QMouseEvent* e) override;
    void move(Viewport& v, QMouseEvent* e) override;
    void release(Viewport& v, QMouseEvent* e) override;

private:
    void stroke(Viewport& v, QMouseEvent* e, bool first);

    bool active_ = false;
    rf::Vector3 prevHit_;
    QElapsedTimer clock_;
    OrbitTool fallback_; // wheel/right button still move the camera during a stroke
};
