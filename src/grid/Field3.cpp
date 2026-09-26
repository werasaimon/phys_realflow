#include "grid/Field3.h"

#include <algorithm>

namespace rf {

// Cell index and weight along one axis of a trilinear lookup (clamped to the samples).
static inline void axisCoord(float f, int n, int& i0, int& i1, float& t) {
    if (n <= 1) { i0 = i1 = 0; t = 0; return; }
    f = clampv(f, 0.0f, float(n - 1));
    i0 = std::min(int(f), n - 2);
    i1 = i0 + 1;
    t = f - float(i0);
}

float Field3::sample(const Vector3& gp) const {
    int i0, i1, j0, j1, k0, k1;
    float tx, ty, tz;
    axisCoord(gp.x - offset.x, nx, i0, i1, tx);
    axisCoord(gp.y - offset.y, ny, j0, j1, ty);
    axisCoord(gp.z - offset.z, nz, k0, k1, tz);
    float c00 = at(i0, j0, k0) * (1 - tx) + at(i1, j0, k0) * tx;
    float c10 = at(i0, j1, k0) * (1 - tx) + at(i1, j1, k0) * tx;
    float c01 = at(i0, j0, k1) * (1 - tx) + at(i1, j0, k1) * tx;
    float c11 = at(i0, j1, k1) * (1 - tx) + at(i1, j1, k1) * tx;
    float c0 = c00 * (1 - ty) + c10 * ty;
    float c1 = c01 * (1 - ty) + c11 * ty;
    return c0 * (1 - tz) + c1 * tz;
}

void Field3::sampleMinMax(const Vector3& gp, float& lo, float& hi) const {
    int i0, i1, j0, j1, k0, k1;
    float tx, ty, tz;
    axisCoord(gp.x - offset.x, nx, i0, i1, tx);
    axisCoord(gp.y - offset.y, ny, j0, j1, ty);
    axisCoord(gp.z - offset.z, nz, k0, k1, tz);
    float vals[8] = {at(i0, j0, k0), at(i1, j0, k0), at(i0, j1, k0), at(i1, j1, k0),
                     at(i0, j0, k1), at(i1, j0, k1), at(i0, j1, k1), at(i1, j1, k1)};
    lo = hi = vals[0];
    for (float v : vals) { lo = std::min(lo, v); hi = std::max(hi, v); }
}

} // namespace rf
