#pragma once
// What an observer near a black hole sees: for every pixel of a camera the null geodesic is
// integrated backwards (towards where the light came from) until it escapes to the sky, falls
// into the hole, or crosses a thin equatorial disk. The reference renderer on the CPU - the way
// DNGR rendered Interstellar (James, von Tunzelmann, Franklin & Thorne 2015, Class. Quantum Grav.
// 32, 065001) and GYOTO (Vincent et al. 2011) trace images; a GPU port of the same integrator is
// the editor's business and is checked against this pixel by pixel. Output is an array of RGB
// floats; nothing here writes files except the PPM helper the tests use.
#include "relativity/Geodesic.h"

#include <string>
#include <vector>

namespace rf {

struct Image {
    int width = 0, height = 0;
    std::vector<float> rgb; // width * height * 3, row by row from the top
    float* pixel(int x, int y) { return &rgb[size_t(3) * (size_t(y) * width + x)]; }
    const float* pixel(int x, int y) const { return &rgb[size_t(3) * (size_t(y) * width + x)]; }
};

class RayTracer {
public:
    struct Camera {
        double r = 30, theta = 1.4, phi = 0; // Boyer-Lindquist position of the observer
        double fovDeg = 20;                  // horizontal field of view
        int width = 96, height = 96;
    };
    struct Disk {
        bool enabled = false;
        double rIn = 0;   // 0: the ISCO
        double rOut = 20; // in units of M
    };
    enum class Hit { Sky, Hole, Disk, Lost };

    explicit RayTracer(const Metric& metric) : metric_(metric), geodesic_(metric) {}
    Camera camera;
    Disk disk;
    Geodesic::Options options; // tol, rMax and the limits of every ray

    Image render() const;
    // One ray: the pixel's direction in the camera plane (sx, sy in [-1, 1]), returns what it hit
    // and the colour. Exposed for the tests (the shadow edge) and the GPU comparison.
    Hit trace(double sx, double sy, float rgb[3]) const;

    // Sky and disk colours, for the tests to reason about pixels.
    static void skyColor(double theta, double phi, float rgb[3]); // a 15-degree checkerboard
    static void writePpm(const Image& img, const std::string& path);

private:
    const Metric& metric_;
    Geodesic geodesic_;
};

} // namespace rf
