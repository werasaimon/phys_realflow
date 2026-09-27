#pragma once
// A spacetime is anything that answers one question: "what is the metric g_mn at this point?"
// The metric is the ruler of the world: ds^2 = g_mn dx^m dx^n turns a small step dx into a length
// (or a time) - and everything else, the Christoffel symbols, the curvature, the tides, the
// geodesics, follows from it by differentiation (Curvature.h). Give your own metric by writing a
// class with dimension() and metric(); the built-in ones below are the textbook cases, each with
// the known answer its test checks:
//
//   MinkowskiCartesian    flat, ds^2 = -dt^2 + dx^2 + dy^2 + dz^2          curvature 0
//   MinkowskiSpherical    the same flat world in (t, r, theta, phi)        Gamma != 0, curvature 0
//   TwoSphere             a sphere's surface of radius R, (theta, phi)     scalar curvature 2 / R^2
//   SchwarzschildSpacetime  a mass M                                       Ricci 0, K = 48 M^2 / r^6
//   KerrSpacetime         a spinning mass (wraps Metric, bit-identical)    Ricci 0, K of Henry 2000
//   ReissnerNordstrom     a charged mass M, Q                              R = 0 but Ricci != 0
//   DeSitterStatic        empty space with a cosmological constant Lambda  R = 4 Lambda
//   FlatFLRW              an expanding flat universe, a(t)                 G_tt = 3 (a'/a)^2 (Friedmann)
//   EllisWormhole         Morris-Thorne with b(r) = b0^2 / r               rho = -b0^2 / (8 pi r^4) < 0
//
// Geometric units G = c = 1 (MTW sign convention, signature -+++). The derivatives of the metric
// are taken by 4th-order central differences with differenceStep() unless a class gives them
// exactly (metricDerivatives); see Curvature.h for the error budget.
#include "relativity/Metric.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <string>

namespace rf {

class Spacetime {
public:
    virtual ~Spacetime() = default;
    virtual int dimension() const = 0;           // 2, 3 or 4
    virtual std::string name() const = 0;
    virtual std::string coordinates() const = 0; // the names of x^0 .. x^(n-1), e.g. "t r theta phi"
    // g_mn at x into g[m * n + n'] (row-major, symmetric).
    virtual void metric(const double* x, double* g) const = 0;
    // d_k g_mn into dg[(k * dim + m) * dim + n], if the class knows them exactly; false: the
    // curvature code takes finite differences instead.
    virtual bool metricDerivatives(const double* x, double* dg) const { return false; }
    // The finite-difference step in coordinate k at x: 1e-3 of the coordinate's size (at least 1),
    // the measured optimum for the curvature: truncation (~h^4) and rounding (~1e-16 / h^2) of the
    // nested second derivatives meet there at ~2e-10 relative (the budget is in Curvature.h).
    virtual double differenceStep(int k, const double* x) const { return 1e-3 * std::max(1.0, std::fabs(x[k])); }
};

class MinkowskiCartesian : public Spacetime {
public:
    int dimension() const override { return 4; }
    std::string name() const override { return "Minkowski (t, x, y, z)"; }
    std::string coordinates() const override { return "t x y z"; }
    void metric(const double* x, double* g) const override;
};

class MinkowskiSpherical : public Spacetime {
public:
    int dimension() const override { return 4; }
    std::string name() const override { return "Minkowski (t, r, theta, phi)"; }
    std::string coordinates() const override { return "t r theta phi"; }
    void metric(const double* x, double* g) const override;
};

class TwoSphere : public Spacetime {
public:
    explicit TwoSphere(double radius) : R(radius) {}
    double R;
    int dimension() const override { return 2; }
    std::string name() const override { return "2-sphere"; }
    std::string coordinates() const override { return "theta phi"; }
    void metric(const double* x, double* g) const override;
};

class SchwarzschildSpacetime : public Spacetime {
public:
    explicit SchwarzschildSpacetime(double mass) : M(mass) {}
    double M;
    int dimension() const override { return 4; }
    std::string name() const override { return "Schwarzschild"; }
    std::string coordinates() const override { return "t r theta phi"; }
    void metric(const double* x, double* g) const override;
};

// The Kerr metric of relativity/Metric (Boyer-Lindquist or Eddington-Finkelstein): the very same
// numbers the geodesic integrator and the ray tracer use.
class KerrSpacetime : public Spacetime {
public:
    explicit KerrSpacetime(const Metric& m) : kerr(m) {}
    Metric kerr;
    int dimension() const override { return 4; }
    std::string name() const override { return "Kerr"; }
    std::string coordinates() const override { return "t r theta phi"; }
    void metric(const double* x, double* g) const override;
};

class ReissnerNordstrom : public Spacetime {
public:
    ReissnerNordstrom(double mass, double charge) : M(mass), Q(charge) {}
    double M, Q;
    int dimension() const override { return 4; }
    std::string name() const override { return "Reissner-Nordstrom"; }
    std::string coordinates() const override { return "t r theta phi"; }
    void metric(const double* x, double* g) const override;
};

// The static patch of de Sitter space: f = 1 - Lambda r^2 / 3 (inside the horizon r < sqrt(3 / Lambda)).
class DeSitterStatic : public Spacetime {
public:
    explicit DeSitterStatic(double lambda) : Lambda(lambda) {}
    double Lambda;
    int dimension() const override { return 4; }
    std::string name() const override { return "de Sitter (static patch)"; }
    std::string coordinates() const override { return "t r theta phi"; }
    void metric(const double* x, double* g) const override;
};

// ds^2 = -dt^2 + a(t)^2 (dx^2 + dy^2 + dz^2): a flat expanding universe. The scale factor is any
// function; a(t) = t^(2/3) is a matter-dominated universe, t^(1/2) a radiation-dominated one.
class FlatFLRW : public Spacetime {
public:
    explicit FlatFLRW(std::function<double(double)> scaleFactor) : a(std::move(scaleFactor)) {}
    std::function<double(double)> a;
    int dimension() const override { return 4; }
    std::string name() const override { return "flat FLRW"; }
    std::string coordinates() const override { return "t x y z"; }
    void metric(const double* x, double* g) const override;
};

// The Morris-Thorne wormhole with shape function b(r) = b0^2 / r (Ellis 1973; Morris & Thorne
// 1988, Am. J. Phys. 56, 395), written in the proper radial distance l (r^2 = l^2 + b0^2), which
// is regular at the throat l = 0:  ds^2 = -dt^2 + dl^2 + (l^2 + b0^2) dOmega^2.
class EllisWormhole : public Spacetime {
public:
    explicit EllisWormhole(double throat) : b0(throat) {}
    double b0;
    int dimension() const override { return 4; }
    std::string name() const override { return "Morris-Thorne (Ellis) wormhole"; }
    std::string coordinates() const override { return "t l theta phi"; }
    void metric(const double* x, double* g) const override;
};

} // namespace rf
