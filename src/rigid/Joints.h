#pragma once
// Joints on a generic Jacobian-row solver.
//
// Every joint describes its constraints as rows of the Jacobian J = [lin_A, ang_A, lin_B, ang_B]:
//   * velocity stage: sequential impulses per row, J v -> target, with accumulated impulse clamped
//     to [lo, hi] (limits, motors, ropes) and warm starting;
//   * position stage: nonlinear Gauss-Seidel (pseudo impulses applied directly to the poses) that
//     removes the drift C(x, q) through the effective mass K = J M^-1 J^T - no Baumgarte term in
//     the velocities, so the joints do not inject energy. Angular errors use log on SO(3).
//
// b = -1 attaches a joint to the static world (its anchor/axis are then given in world space).

#include "rigid/RigidBody.h"

#include <limits>
#include <vector>

namespace rf {

enum class JointType { Ball, Hinge, Slider, Fixed, Distance };

struct JacobianRow {
    Vector3 linA, angA, linB, angB;
    float effMass = 0;   // 1 / (J M^-1 J^T)
    float bias = 0;      // target: J v = -bias
    float lo = -std::numeric_limits<float>::infinity();
    float hi = std::numeric_limits<float>::infinity();
    float softness = 0;  // gamma of a soft constraint (springs), 0 = rigid
    float lambda = 0;    // accumulated impulse
    bool motor = false; // marks actuator work separately from passive constraint rows
};

class Joint {
    friend class RigidStepCheckpoint;
public:
    Joint(int a, int b) : a(a), b(b) {}
    virtual ~Joint() = default;
    virtual JointType type() const = 0;
    virtual const char* name() const = 0;

    // Velocity stage.
    void prepare(std::vector<RigidBody>& bodies, float h, bool warmStart, bool measureMotorWork = false);
    void solveVelocity(std::vector<RigidBody>& bodies);
    // Position stage; returns the remaining error (m or rad) for convergence checks.
    virtual float solvePosition(std::vector<RigidBody>& bodies) = 0;

    Vector3 worldAnchorA(const std::vector<RigidBody>& bodies) const;
    Vector3 worldAnchorB(const std::vector<RigidBody>& bodies) const;
    Vector3 worldAxis(const std::vector<RigidBody>& bodies) const; // axis of A (hinge / slider)
    float appliedImpulse() const;
    double motorWork() const { return motorWork_; } // signed discrete work of the last prepared trial [J]
    void captureMotorMotion(const std::vector<RigidBody>& bodies);
    void finishMotorWork(const std::vector<RigidBody>& bodies);

    int a, b;
    Vector3 localAnchorA, localAnchorB; // B's anchor is in world space when b < 0
    Vector3 localAxisA, localAxisB;     // unit axis in each body frame (world for b < 0)
    Quaternion refRel;                     // qB^-1 qA at creation (fixed / slider)

protected:
    virtual void buildRows(const std::vector<RigidBody>& bodies, float h) = 0;
    void addRow(const std::vector<RigidBody>& bodies, JacobianRow r);
    void applyRow(std::vector<RigidBody>& bodies, const JacobianRow& row, float impulse);
    // Helpers for the position stage.
    float correctPoint(std::vector<RigidBody>& bodies, const Vector3& C, const Vector3& pA, const Vector3& pB, float beta);
    float correctAngle(std::vector<RigidBody>& bodies, Vector3 e, float beta);
    float correctAlong(std::vector<RigidBody>& bodies, const Vector3& n, float C, const Vector3& pA, const Vector3& pB, float beta);
    Vector3 rotationError(const std::vector<RigidBody>& bodies) const; // log(qRel * refRel^-1), world

    std::vector<JacobianRow> rows_;
    std::vector<float> warm_;
    bool measureMotorWork_ = false;
    double motorWork_ = 0;
    Vector3 motorStartLinA_, motorStartAngA_, motorStartLinB_, motorStartAngB_;
};

// 3 linear rows: a point of A coincides with a point of B.
class BallJoint : public Joint {
public:
    BallJoint(int a, int b) : Joint(a, b) {}
    JointType type() const override { return JointType::Ball; }
    const char* name() const override { return "Шаровое"; }
    float solvePosition(std::vector<RigidBody>& bodies) override;

protected:
    void buildRows(const std::vector<RigidBody>& bodies, float h) override;
};

// Ball + 2 angular rows keeping the axes aligned; optional angle limits and motor.
class HingeJoint : public BallJoint {
public:
    HingeJoint(int a, int b) : BallJoint(a, b) {}
    JointType type() const override { return JointType::Hinge; }
    const char* name() const override { return "Шарнир"; }
    float solvePosition(std::vector<RigidBody>& bodies) override;
    float angle(const std::vector<RigidBody>& bodies) const;

    bool limitEnabled = false;
    float lower = 0, upper = 0; // [rad]
    bool motorEnabled = false;
    float motorSpeed = 0;       // [rad/s], relative A-B about the axis
    float maxMotorTorque = 0;   // [N m]
    Vector3 localRefA, localRefB;  // perpendicular reference vectors for the angle

protected:
    void buildRows(const std::vector<RigidBody>& bodies, float h) override;
};

// 3 angular rows + 2 linear rows perpendicular to the axis; optional travel limits.
class SliderJoint : public Joint {
public:
    SliderJoint(int a, int b) : Joint(a, b) {}
    JointType type() const override { return JointType::Slider; }
    const char* name() const override { return "Призматическое"; }
    float solvePosition(std::vector<RigidBody>& bodies) override;
    float translation(const std::vector<RigidBody>& bodies) const;

    bool limitEnabled = false;
    float lower = 0, upper = 0; // [m] along the axis

protected:
    void buildRows(const std::vector<RigidBody>& bodies, float h) override;
};

// Ball + 3 angular rows: the relative pose is frozen.
class FixedJoint : public BallJoint {
public:
    FixedJoint(int a, int b) : BallJoint(a, b) {}
    JointType type() const override { return JointType::Fixed; }
    const char* name() const override { return "Сварка"; }
    float solvePosition(std::vector<RigidBody>& bodies) override;

protected:
    void buildRows(const std::vector<RigidBody>& bodies, float h) override;
};

// One row along the line between the anchors: rigid rod, rope (pull only) or soft spring.
class DistanceJoint : public Joint {
public:
    DistanceJoint(int a, int b) : Joint(a, b) {}
    JointType type() const override { return JointType::Distance; }
    const char* name() const override { return "Дистанция"; }
    float solvePosition(std::vector<RigidBody>& bodies) override;

    float length = 1;
    bool rope = false;      // only resists stretching
    float frequency = 0;    // > 0: spring (Hz), solved as a soft constraint
    float dampingRatio = 0.5f;

protected:
    void buildRows(const std::vector<RigidBody>& bodies, float h) override;
};

} // namespace rf
