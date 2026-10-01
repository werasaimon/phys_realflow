#pragma once
// Same physical benchmark inputs; only the stepping method is changed explicitly.
#include "BouncingSdk.h"

namespace rf::benchmark {
class BouncingVariational : public BouncingSdk {
public:
    BouncingVariational(int iterations, bool) : BouncingSdk(iterations, true, true) {}
};
} // namespace rf::benchmark
