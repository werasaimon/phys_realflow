#pragma once
// Input uncertainty by Monte Carlo: the inputs of a model that are not known exactly (a friction
// coefficient to 5 %, a density to 1 %, the size of an experiment's water column) are drawn from
// their distributions, the model runs once per draw, and the spread of the outputs is the input
// uncertainty of the answer (ASME V&V 20-2009, u_input; JCGM 101:2008, "Propagation of
// distributions using a Monte Carlo method").
//
// The draws are a Latin hypercube (McKay, Beckman & Conover 1979): each input's probability
// range is cut into n equal strata, every stratum is sampled once, and the strata of different
// inputs are paired at random. A few dozen runs then cover the range much more evenly than plain
// random draws. The seed is fixed, so the same study gives the same numbers.

#include <cstdint>
#include <functional>
#include <string>
#include <vector>

namespace rf::verify {

struct UncertainInput {
    enum class Distribution { Normal, Uniform };
    std::string name;
    double nominal = 0;
    double spread = 0; // Normal: the standard deviation; Uniform: the half-width
    Distribution distribution = Distribution::Normal;
};

struct UqSummary {
    int samples = 0;
    double mean = 0;
    double stddev = 0;         // the standard input uncertainty u_input
    double lo95 = 0, hi95 = 0; // the 2.5 % and 97.5 % quantiles of the outputs
    std::vector<double> outputs;
};

// Runs model(inputs) for `samples` Latin-hypercube draws of the inputs (in the order given).
UqSummary propagate(const std::vector<UncertainInput>& inputs, const std::function<double(const std::vector<double>&)>& model,
                    int samples, uint32_t seed);

// The standard normal quantile: the x with P(Z < x) = p (Acklam's rational approximation,
// relative error below 1.2e-9).
double normalQuantile(double p);

} // namespace rf::verify
