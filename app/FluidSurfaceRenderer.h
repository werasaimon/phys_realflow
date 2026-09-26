#pragma once
// Screen-space rendering of a particle liquid as a smooth water surface
// (Green 2010, "Screen Space Fluid Rendering with Curvature Flow", NVIDIA; van der Laan, Green,
// Sainz 2009). Per frame, in screen space - no mesh is built:
//  1. depth:     every particle as a sphere, the nearest surface per pixel (eye-space distance);
//  2. thickness: the spheres added up along the view ray - how much water the light crosses;
//  3. smoothing: the bumpy sphere depth is filtered into a smooth surface with the narrow-range
//                filter of Truong & Yuksel 2018 ("A Narrow-Range Filter for Screen-Space Fluid
//                Rendering"): a separable Gaussian that keeps only samples within a depth range
//                around the pixel - it smooths the surface but not across its silhouettes;
//  4. shading:   normals from the smoothed depth; Fresnel reflection of the sky (Schlick, water
//                n = 1.33), refraction of the scene behind (offset along the normal), Beer-Lambert
//                absorption through the thickness (red is absorbed first - deep water turns blue),
//                a sun highlight. The surface writes its depth, so bodies in front of it hide it.
// OpenGL 3.0 / GLSL 1.30, so it also runs on the CPU (llvmpipe) - there at half resolution.

#include "math/Vector3.h"

#include <QMatrix4x4>
#include <QOpenGLShaderProgram>
#include <QString>

#include <vector>

class QOpenGLExtraFunctions;

class FluidSurfaceRenderer {
public:
    struct Settings {
        float sphereScale = 1.6f;                 // sprite radius / particle radius (the spheres must overlap)
        float smoothing = 5.0f;                   // filter width in particle radii
        // [1/m] per colour (R, G, B): Beer-Lambert extinction. Pure water ~(0.45, 0.07, 0.02) - too
        // faint to see in a tank this size; these are of slightly turbid (pool / sea) water.
        QVector3D absorption{3.0f, 0.9f, 0.45f};
        QVector3D scatterColor{0.05f, 0.25f, 0.35f}; // colour of light scattered back inside the water
        float refraction = 0.04f;                 // screen offset of what is seen through the water
    };
    Settings settings;

    // Needs the current OpenGL context; glslHeader is the #version line of the context.
    void initialize(const QString& glslHeader);
    void release();

    // Draws the liquid (particle centres, radius) into sceneFbo - the framebuffer that already holds
    // the opaque scene (its colour is refracted through the water, its depth hides the water).
    // lowQuality: half-resolution buffers and fewer filter taps (the CPU renderer).
    void render(const std::vector<rf::Vector3>& particles, float radius, const QMatrix4x4& view, const QMatrix4x4& proj,
                int width, int height, unsigned sceneFbo, bool lowQuality);

private:
    // An offscreen target: one colour texture and, optionally, a depth buffer.
    struct Target {
        unsigned fbo = 0, color = 0, depth = 0;
        int w = 0, h = 0;
    };
    void resize(Target& t, int w, int h, int internalFormat, int format, int type, bool withDepth);
    void destroy(Target& t);
    void drawFullScreen();

    QOpenGLExtraFunctions* gl_ = nullptr;
    QOpenGLShaderProgram depthProg_, thickProg_, smoothProg_, shadeProg_;
    unsigned vao_ = 0, vbo_ = 0, emptyVao_ = 0;
    Target depth_, thickness_, smoothA_, smoothB_, scene_;
};
