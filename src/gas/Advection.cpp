// Advection of GasSolver: semi-Lagrangian back-tracing (RK2) of the face velocities and the
// cell scalars, with the MacCormack second-order correction and its limiter (Selle et al. 2008).
#include "gas/GasSolver.h"
#include "gas/GasSolverInternal.h"

#include "core/Parallel.h"

namespace rf {

// ---------------------------------------------------------------------------
// Advection
// ---------------------------------------------------------------------------
Vector3 GasSolver::sampleVelGrid(const Vector3& gp) const { return sampleVel(u0_, v0_, w0_, gp) / dx_; }

Vector3 GasSolver::backtrace(const Vector3& gp, float dt) const {
    Vector3 v1 = sampleVelGrid(gp);
    Vector3 mid = gp - v1 * (0.5f * dt);
    Vector3 v2 = sampleVelGrid(mid);
    return gp - v2 * dt;
}

void GasSolver::advect(const Field3& src, Field3& dst, float dt, bool scalar) const {
    dst.nx = src.nx; dst.ny = src.ny; dst.nz = src.nz; dst.offset = src.offset;
    dst.d.resize(src.d.size());
    const int sx = src.nx, sy = src.ny, sz = src.nz;
    const float size[3] = {float(nx_), float(ny_), float(nz_)};
    // Did the backtraced point come in from outside through an open side of the domain?
    auto fromOutside = [&](const Vector3& b) {
        for (int a = 0; a < 3; ++a) {
            if (b[a] < 0 && params.bc[2 * a] != BoundaryType::Wall) return true;
            if (b[a] > size[a] && params.bc[2 * a + 1] != BoundaryType::Wall) return true;
        }
        return false;
    };
    parallelFor(int(sz), [&](int k_) {
        int k = k_;
        for (int j = 0; j < sy; ++j)
            for (int i = 0; i < sx; ++i) {
                Vector3 gp = Vector3(float(i), float(j), float(k)) + src.offset;
                Vector3 b = backtrace(gp, dt);
                dst.at(i, j, k) = scalar && fromOutside(b) ? 0.0f : src.sample(b);
            }
    }, 1);
}

void GasSolver::advectField(Field3& f, Field3& t1, Field3& t2, float dt, bool scalar) {
    advect(f, t1, dt, scalar);
    if (!params.maccormack) {
        f.d.swap(t1.d);
        return;
    }
    advect(t1, t2, -dt, scalar);
    const int sx = f.nx, sy = f.ny, sz = f.nz;
    parallelFor(int(sz), [&](int k_) {
        int k = k_;
        for (int j = 0; j < sy; ++j)
            for (int i = 0; i < sx; ++i) {
                size_t id = f.idx(i, j, k);
                float val = t1.d[id] + 0.5f * (f.d[id] - t2.d[id]);
                float lo, hi;
                f.sampleMinMax(backtrace(Vector3(float(i), float(j), float(k)) + f.offset, dt), lo, hi);
                t2.d[id] = (val < lo || val > hi) ? t1.d[id] : val;
            }
    }, 1);
    f.d.swap(t2.d);
}

void GasSolver::advectScalars(const std::vector<Field3*>& fields, float dt) {
    // All cell-centred scalars ride the same flow: the departure points (and the arrival points
    // for MacCormack's backward step) are traced once for all of them. What arrives from outside
    // through an open side is clean ambient gas (0).
    const size_t n = size_t(nx_) * ny_ * nz_;
    departure_.resize(n);
    arrival_.resize(n);
    fromOutside_.resize(n);
    const float size[3] = {float(nx_), float(ny_), float(nz_)};
    auto outside = [&](const Vector3& b) {
        for (int a = 0; a < 3; ++a) {
            if (b[a] < 0 && params.bc[2 * a] != BoundaryType::Wall) return true;
            if (b[a] > size[a] && params.bc[2 * a + 1] != BoundaryType::Wall) return true;
        }
        return false;
    };
    parallelFor(nz_, [&](int k) {
        for (int j = 0; j < ny_; ++j)
            for (int i = 0; i < nx_; ++i) {
                const size_t c = cidx(i, j, k);
                const Vector3 gp(i + 0.5f, j + 0.5f, k + 0.5f);
                departure_[c] = backtrace(gp, dt);
                arrival_[c] = backtrace(gp, -dt);
                fromOutside_[c] = uint8_t(outside(departure_[c]) | (outside(arrival_[c]) << 1));
            }
    }, 1);
    for (Field3* field : fields) {
        Field3& f = *field;
        t1_.init(nx_, ny_, nz_, f.offset);
        t2_.init(nx_, ny_, nz_, f.offset);
        parallelFor(int(n), [&](int c) { t1_.d[c] = (fromOutside_[c] & 1) ? 0.0f : f.sample(departure_[c]); }, 4096);
        if (params.maccormack) {
            // MacCormack: trace the result back, correct by half the round-trip error, fall back to
            // the plain value where the correction leaves the range of the departure cell (limiter).
            parallelFor(int(n), [&](int c) { t2_.d[c] = (fromOutside_[c] & 2) ? 0.0f : t1_.sample(arrival_[c]); }, 4096);
            parallelFor(int(n), [&](int c) {
                if (fromOutside_[c] & 1) return;
                const float val = t1_.d[c] + 0.5f * (f.d[c] - t2_.d[c]);
                float lo, hi;
                f.sampleMinMax(departure_[c], lo, hi);
                if (val >= lo && val <= hi) t1_.d[c] = val;
            }, 4096);
        }
        f.d.swap(t1_.d);
    }
}

} // namespace rf
