// Wall-clock measurements around native stepping, excluding topology/import observers.
#pragma once
#include "core/Format.h"
#include <chrono>
#include <string>

struct ReferenceTiming {
    using Clock = std::chrono::steady_clock;
    Clock::time_point previous;
    double engineMs = 0, observerMs = 0;
    void start() { previous = Clock::now(); }
    double lap() {
        const auto now = Clock::now();
        const double ms = std::chrono::duration<double, std::milli>(now - previous).count();
        previous = now;
        return ms;
    }
    void finishEngine() { engineMs += lap(); }
    void finishObserver() { observerMs += lap(); }
    std::string csv() const { return "," + rf::numberText(engineMs) + "," + rf::numberText(observerMs); }
};
