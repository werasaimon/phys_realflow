#include "grid/Tokamak.h"

#include <cmath>
#include <complex>

namespace rf {

static constexpr float kMu0 = MagneticField::kMu0;

float Tokamak::plasmaCurrent() const {
    // q_a = a B0 / (R0 B_theta(a)) with B_theta(a) = mu0 I_p / (2 pi a).
    return 2.0f * kPi * minorRadius * minorRadius * toroidalField / (kMu0 * majorRadius * std::max(safetyFactorEdge, 1e-3f));
}

float Tokamak::poloidalField(float r) const {
    const float a = minorRadius, Ip = plasmaCurrent();
    if (r >= a) return kMu0 * Ip / (2.0f * kPi * std::max(r, 1e-6f));
    return kMu0 * Ip * r / (2.0f * kPi * a * a);
}

float Tokamak::safetyFactor(float r) const {
    if (r <= minorRadius) return safetyFactorEdge; // constant current: constant twist
    return safetyFactorEdge * r * r / (minorRadius * minorRadius);
}

float Tokamak::equilibriumShift() const {
    const float betaPoloidal = 0.0f; // a cold plasma: no pressure to hold
    const float lambda = betaPoloidal + 0.5f * internalInductance() - 1.0f;
    const float a = minorRadius, b = vesselRadius;
    return b * b / (2.0f * majorRadius) * (std::log(b / a) + (lambda + 0.5f) * (1.0f - a * a / (b * b)));
}

float Tokamak::verticalFieldStrength() const {
    if (!verticalField) return 0.0f;
    return kMu0 * plasmaCurrent() * equilibriumShift() / (2.0f * kPi * vesselRadius * vesselRadius);
}

float Tokamak::wallLimit() const {
    const float a2 = minorRadius * minorRadius, b2 = vesselRadius * vesselRadius;
    return 2.0f * a2 / (a2 + b2);
}

float Tokamak::kinkGrowthRate(float density) const {
    if (!kinkUnstable()) return 0.0f;
    const float a = minorRadius, a2 = a * a, b2 = vesselRadius * vesselRadius, q = safetyFactorEdge;
    const float vA = poloidalField(a) / std::sqrt(kMu0 * std::max(density, 1e-12f)); // on the poloidal field at the edge
    return vA / a * std::sqrt((1.0f - q) * (q * (b2 + a2) - 2.0f * a2) / b2);
}

bool Tokamak::inside(const Vector3& x) const {
    const Vector3 p = x - centre;
    const float R = std::hypot(p.x, p.z);
    return std::hypot(R - majorRadius, p.y) < vesselRadius;
}

void Tokamak::channelCentre(float phi, float& R, float& y) const {
    // Helix m = 1, n = 1: the offset turns once in the poloidal plane per turn around the torus,
    // in the sense of the field lines' twist (theta = phi for q > 0).
    const float d = seedDisplacement * minorRadius;
    R = majorRadius + d * std::cos(phi);
    y = d * std::sin(phi);
}

float Tokamak::channelRadius(const Vector3& x) const {
    const Vector3 p = x - centre;
    const float R = std::hypot(p.x, p.z), phi = std::atan2(p.z, p.x);
    float cR, cy;
    channelCentre(phi, cR, cy);
    return std::hypot(R - cR, p.y - cy);
}

void Tokamak::controlGains(float density, float& gain, float& damping) const {
    // The shell's restoring stiffness per unit length k = mu0 I^2 / (2 pi b^2); the control adds
    // 4 k (a change dB_v of the vertical field is a force -I dB_v); the ring plus the fluid it
    // pushes around (added mass with the wall) weighs m = rho pi a^2 2 b^2 / (b^2 - a^2) per
    // unit length; critical damping 2 sqrt(5 k m).
    const float I = plasmaCurrent(), a2 = sqr(minorRadius), b2 = sqr(vesselRadius);
    const float k = kMu0 * I * I / (2.0f * kPi * b2);
    const float m = density * kPi * a2 * 2.0f * b2 / (b2 - a2);
    gain = 4.0f * k / I;
    damping = 2.0f * std::sqrt(5.0f * k * m) / I;
}

Vector3 Tokamak::coilPotential(const Vector3& x, float verticalFieldT) const {
    // A = A_y(R) y_hat gives B = (-dA_y/dz, 0, dA_y/dx) = (dA_y/dR) phi_hat: A_y = B0 R0 ln(R / R0)
    // makes B_phi = B0 R0 / R along +phi. A = -B_v x z_hat gives B = B_v y_hat (the vertical field).
    // Inside the inner wall (R < R0 - b, where the coils' legs and the central solenoid are) the
    // 1 / R field is capped at its wall value: no plasma is there, and its Alfven speed would
    // otherwise set the time step of the whole grid.
    const Vector3 p = x - centre;
    const float R = std::max(std::hypot(p.x, p.z), majorRadius - vesselRadius);
    return {0.0f, toroidalField * majorRadius * std::log(R / majorRadius), -verticalFieldT * p.x};
}

// Complete elliptic integrals K(k) and E(k) by the arithmetic-geometric mean.
static void ellipticKE(double k, double& K, double& E) {
    double a = 1.0, b = std::sqrt(1.0 - k * k), c = k, sum = 0.5 * c * c, pow2 = 0.5;
    for (int n = 0; n < 40 && std::fabs(a - b) > 1e-14 * a; ++n) {
        const double an = 0.5 * (a + b);
        c = 0.5 * (a - b);
        b = std::sqrt(a * b);
        a = an;
        pow2 *= 2.0;
        sum += pow2 * c * c;
    }
    K = kPi / (2.0 * a);
    E = K * (1.0 - sum);
}

double Tokamak::loopPotential(double R, double y, double Rl, double yl, double I) {
    // Jackson eq. 5.37: A_phi = mu0 I / (pi k) sqrt(Rl / R) [(1 - k^2/2) K(k) - E(k)],
    // k^2 = 4 Rl R / ((Rl + R)^2 + (y - yl)^2). Near the axis (k -> 0) the bracket is pi k^4 / 32
    // and the difference of two near-equal numbers: its series is used there.
    if (R < 1e-9) return 0.0;
    const double dy = y - yl;
    const double k2 = std::min(4.0 * Rl * R / ((Rl + R) * (Rl + R) + dy * dy), 1.0 - 1e-12), k = std::sqrt(k2);
    const double front = double(kMu0) * I * std::sqrt(Rl / R);
    if (k2 < 1e-4) return front * k * k2 / 32.0;
    double K, E;
    ellipticKE(k, K, E);
    return front / (kPi * k) * ((1.0 - 0.5 * k2) * K - E);
}

Vector3 Tokamak::plasmaPotential(const Vector3& x) const {
    // The channel as a lattice of current loops (spacing a / 12, each with an equal share of
    // I_p: the constant profile); their potentials add. The helical seed displaces the whole
    // configuration: A(x) = A0(x - xi), xi = d (cos phi, sin phi) in the (R, y) plane.
    const int M = 12;
    const float a = minorRadius, s = a / M, Ip = plasmaCurrent();
    const Vector3 p = x - centre;
    const float R = std::hypot(p.x, p.z), phi = std::atan2(p.z, p.x);
    float cR, cy;
    channelCentre(phi, cR, cy);
    const double Rq = R - (cR - majorRadius), yq = p.y - cy; // the point, shifted back by the seed
    int count = 0;
    for (int j = -M; j <= M; ++j)
        for (int i = -M; i <= M; ++i) count += (i * i + j * j) * (s * s) <= a * a;
    double A = 0;
    for (int j = -M; j <= M; ++j)
        for (int i = -M; i <= M; ++i) {
            if ((i * i + j * j) * (s * s) > a * a) continue;
            A += loopPotential(Rq, yq, majorRadius + i * s, j * s, Ip / count);
        }
    const Vector3 phiHat(-std::sin(phi), 0.0f, std::cos(phi));
    return phiHat * float(A);
}

float Tokamak::tracer(const Vector3& x) const {
    const float a = minorRadius;
    return clampv((1.05f * a - channelRadius(x)) / (0.2f * a), 0.0f, 1.0f);
}

float Tokamak::resistivityAt(const Vector3& x, float plasmaResistivity) const {
    // The plasma's resistivity out to 1.2 a (room for the column to move), the vacuum's from
    // 1.4 a, a linear ramp between.
    const float f = clampv((channelRadius(x) - 1.2f * minorRadius) / (0.2f * minorRadius), 0.0f, 1.0f);
    return plasmaResistivity + f * (vacuumResistivity - plasmaResistivity);
}

float Tokamak::measuredCurrent(const MagneticField& m, float dx) const {
    // Ampere: I = (1 / mu0) closed integral of B . dl around a circle of radius r_c about the magnetic
    // axis in the plane phi = 0 (the (x, y) plane) - r_c inside the wall, outside the channel.
    const float rc = vesselRadius - 1.5f * dx;
    const int N = 96;
    double circulation = 0;
    for (int n = 0; n < N; ++n) {
        const float th = 2.0f * kPi * (n + 0.5f) / N;
        const Vector3 p = centre + Vector3(majorRadius + rc * std::cos(th), rc * std::sin(th), 0.0f);
        const Vector3 thetaHat(-std::sin(th), std::cos(th), 0.0f);
        circulation += dot(m.fieldAt(p), thetaHat) * (2.0f * kPi * rc / N);
    }
    return float(circulation / kMu0);
}

void Tokamak::currentCentroids(const MagneticField& m, float dx, std::vector<Vector2>& out) const {
    const int N = 24;
    const float b = vesselRadius - dx, sign = plasmaCurrent() >= 0 ? 1.0f : -1.0f;
    out.clear();
    for (int n = 0; n < N; ++n) {
        const float phi = 2.0f * kPi * float(n) / float(N);
        const Vector3 Rhat(std::cos(phi), 0.0f, std::sin(phi)), phiHat(-std::sin(phi), 0.0f, std::cos(phi));
        double w = 0, sR = 0, sy = 0;
        for (float y = -b; y <= b; y += dx)
            for (float R = majorRadius - b; R <= majorRadius + b; R += dx) {
                if (std::hypot(R - majorRadius, y) >= b) continue;
                const float j = std::max(0.0f, sign * dot(m.currentAt(centre + Rhat * R + Vector3(0.0f, y, 0.0f)), phiHat));
                w += j;
                sR += j * (R - majorRadius);
                sy += j * y;
            }
        out.push_back(w > 0 ? Vector2(float(sR / w), float(sy / w)) : Vector2(0.0f, 0.0f));
    }
}

float Tokamak::measuredShift(const MagneticField& m, float dx) const {
    std::vector<Vector2> c;
    currentCentroids(m, dx, c);
    float s = 0;
    for (const Vector2& d : c) s += d.x;
    return c.empty() ? 0.0f : s / float(c.size());
}

float Tokamak::kinkAmplitude(const MagneticField& m, float dx) const {
    // The n = 1 Fourier amplitude of (dR + i dy) around the torus, for either sense of the helix.
    std::vector<Vector2> c;
    currentCentroids(m, dx, c);
    std::complex<double> plus = 0, minus = 0;
    for (size_t n = 0; n < c.size(); ++n) {
        const float phi = 2.0f * kPi * float(n) / float(c.size());
        const std::complex<double> d(c[n].x, c[n].y), e(std::cos(phi), std::sin(phi));
        plus += d * std::conj(e);
        minus += d * e;
    }
    return c.empty() ? 0.0f : float(std::max(std::abs(plus), std::abs(minus)) / double(c.size()));
}

} // namespace rf
