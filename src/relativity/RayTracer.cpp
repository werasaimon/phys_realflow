// The CPU ray tracer of a black hole: every pixel's photon is followed backwards along its
// geodesic to the sky, the accretion disk or the hole. The reference for the editor's GPU
// version, and the source of the shadow test. See RayTracer.h.
#include "relativity/RayTracer.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace rf {

static const double kPiD = 3.14159265358979323846;

void RayTracer::skyColor(double theta, double phi, float rgb[3]) {
    // 15-degree checkerboard in (theta, phi): the lensing shows as the squares bend.
    double ph = std::fmod(phi, 2.0 * kPiD);
    if (ph < 0) ph += 2.0 * kPiD;
    const int a = int(std::floor(theta / (kPiD / 12.0))), b = int(std::floor(ph / (kPiD / 12.0)));
    const bool light = ((a + b) & 1) == 0;
    rgb[0] = light ? 0.92f : 0.16f;
    rgb[1] = light ? 0.92f : 0.28f;
    rgb[2] = light ? 0.96f : 0.55f;
}

// The photon arriving at the camera through pixel (sx, sy): where it is (the camera) and its
// 4-momentum, from the pixel's direction in the observer's local frame. Forward is towards the
// hole (-e_r), right is +e_phi, up is -e_theta (theta grows downwards). The photon moves along
// -n; its path is followed backwards in the affine parameter to where it came from.
void RayTracer::arrivingPhoton(double sx, double sy, GeodesicState& s) const {
    const double tanHalf = std::tan(0.5 * camera.fovDeg * kPiD / 180.0);
    const double aspect = double(camera.height) / double(camera.width);
    double n[3] = {-1.0, -sy * tanHalf * aspect, sx * tanHalf};
    const double len = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    for (double& v : n) v /= len;
    const double arriving[3] = {-n[0], -n[1], -n[2]};
    s.x[0] = 0.0;
    s.x[1] = camera.r;
    s.x[2] = camera.theta;
    s.x[3] = camera.phi;
    Geodesic::momentumFromLocal(metric_, s.x, arriving, 1.0, 1.0, s.p);
}

// The colour of the disk where the ray crossed it. Redshift factor of light from the Keplerian
// disk at radius r seen far away: g = E / (u^t (E - Omega L)) with
// u^t = 1 / sqrt(-(g_tt + 2 Omega g_tphi + Omega^2 g_phiphi)) (Cunningham 1975); the observed
// intensity is g^4 times the emitted one (I / nu^3 is invariant: g^3 per unit frequency, one
// more g for the frequency width). The disk's own emissivity here is the toy law (r_in / r)^3 -
// a real one is Novikov & Thorne 1973.
void RayTracer::shadeDisk(const GeodesicState& s, double diskR, double rIn, float rgb[3]) const {
    double x[4] = {0.0, diskR, 0.5 * kPiD, 0.0}, g[4][4];
    metric_.covariant(x, g);
    const double Omega = metric_.keplerianOmega(diskR, true);
    const double ut = 1.0 / std::sqrt(-(g[0][0] + 2.0 * Omega * g[0][3] + Omega * Omega * g[3][3]));
    const double E = -s.p[0], L = s.p[3];
    const double gf = E / (ut * (E - Omega * L));
    const double brightness = std::pow(std::max(gf, 0.0), 4.0) * std::pow(rIn / diskR, 3.0);
    rgb[0] = float(std::min(1.0, 1.0 * brightness));
    rgb[1] = float(std::min(1.0, 0.75 * brightness));
    rgb[2] = float(std::min(1.0, 0.4 * brightness));
}

// One pixel: the arriving photon is followed back until it came from the sky, from the disk or
// out of the hole, and the pixel takes the colour of what it came from.
RayTracer::Hit RayTracer::trace(double sx, double sy, float rgb[3]) const {
    GeodesicState s;
    arrivingPhoton(sx, sy, s);
    Geodesic::Options opt = options;
    opt.backward = true;
    opt.mu = 0;
    opt.rMax = std::max(opt.rMax, 2.0 * camera.r);
    // A ray from outside that is still falling below the innermost (prograde) photon orbit has
    // no turning point left: it is in the hole. Ending it there instead of at the horizon saves
    // the spiral of ever smaller steps that Boyer-Lindquist coordinates force near r+.
    if (opt.rMin <= 0) opt.rMin = std::max(metric_.horizonRadius() * (1.0 + 1e-3), 0.98 * metric_.photonSphereRadius(true));
    const double rIn = disk.rIn > 0 ? disk.rIn : metric_.iscoRadius(true);
    double diskR = 0;
    bool diskHit = false;
    auto crossesDisk = [&](const GeodesicState& before, const GeodesicState& after) {
        if (!disk.enabled) return false;
        const double a = before.x[2] - 0.5 * kPiD, b = after.x[2] - 0.5 * kPiD;
        if (a * b > 0) return false;
        const double f = a / (a - b); // where between the two states the plane was crossed
        const double r = before.x[1] + f * (after.x[1] - before.x[1]);
        if (r < rIn || r > disk.rOut) return false;
        diskR = r;
        diskHit = true;
        return true;
    };
    const GeodesicResult res = geodesic_.integrate(s, opt, crossesDisk);
    if (diskHit) {
        shadeDisk(s, diskR, rIn, rgb);
        return Hit::Disk;
    }
    switch (res.end) {
    case GeodesicEnd::Escaped:
        skyColor(s.x[2], s.x[3], rgb);
        return Hit::Sky;
    case GeodesicEnd::Captured:
        rgb[0] = rgb[1] = rgb[2] = 0.0f;
        return Hit::Hole;
    default:
        rgb[0] = rgb[1] = rgb[2] = 0.2f;
        return Hit::Lost;
    }
}

Image RayTracer::render() const {
    Image img;
    img.width = camera.width;
    img.height = camera.height;
    img.rgb.assign(size_t(3) * img.width * img.height, 0.0f);
    for (int y = 0; y < img.height; ++y) {
        const double sy = 1.0 - 2.0 * (y + 0.5) / img.height;
        for (int x = 0; x < img.width; ++x) {
            const double sx = 2.0 * (x + 0.5) / img.width - 1.0;
            trace(sx, sy, img.pixel(x, y));
        }
    }
    return img;
}

void RayTracer::writePpm(const Image& img, const std::string& path) {
    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) return;
    std::fprintf(f, "P6\n%d %d\n255\n", img.width, img.height);
    for (size_t i = 0; i < img.rgb.size(); ++i) {
        const unsigned char b = (unsigned char)std::lround(255.0 * std::min(1.0f, std::max(0.0f, img.rgb[i])));
        std::fwrite(&b, 1, 1, f);
    }
    std::fclose(f);
}

} // namespace rf
