#pragma once
// A scalar field on a regular 3D grid - one component of a staggered (MAC) grid field or a
// cell-centred quantity. Sample (i, j, k) sits at (i, j, k) + offset in cell units.

#include "math/Math.h"

#include <vector>

namespace rf {

struct Field3 {
    int nx = 0, ny = 0, nz = 0;
    Vector3 offset;  // sample location of (i,j,k) in cell units is (i,j,k) + offset
    std::vector<float> d;

    void init(int x, int y, int z, const Vector3& off, float v = 0) {
        nx = x; ny = y; nz = z; offset = off;
        d.assign(size_t(x) * y * z, v);
    }
    size_t idx(int i, int j, int k) const { return size_t(i) + size_t(nx) * (size_t(j) + size_t(ny) * size_t(k)); }
    float& at(int i, int j, int k) { return d[idx(i, j, k)]; }
    float at(int i, int j, int k) const { return d[idx(i, j, k)]; }
    // Trilinear value at gp: position in cell units relative to the grid origin (clamped inside).
    float sample(const Vector3& gp) const;
    // Smallest and largest of the 8 samples around gp (limiter of MacCormack advection).
    void sampleMinMax(const Vector3& gp, float& lo, float& hi) const;
};

} // namespace rf
