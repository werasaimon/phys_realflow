#include "FluidSurfaceRenderer.h"

#include <QOpenGLContext>
#include <QOpenGLExtraFunctions>

#include <algorithm>

// ---------------------------------------------------------------------------
// Shaders (the #version line is added for the context: 330 core or 130)
// ---------------------------------------------------------------------------

// Particles as point sprites the size of their sphere on screen.
static const char* kSpriteVS = R"(
in vec3 aPos;
uniform mat4 uView, uProj;
uniform float uRadius, uPointScale;
out vec3 vCenter;
void main() {
    vec4 vp = uView * vec4(aPos, 1.0);
    vCenter = vp.xyz;
    gl_Position = uProj * vp;
    gl_PointSize = 2.0 * uRadius * uPointScale / max(-vp.z, 1e-3);
})";

// 1. Nearest sphere surface: eye-space distance in the colour, true depth for the depth test.
static const char* kDepthFS = R"(
in vec3 vCenter;
uniform mat4 uProj;
uniform float uRadius;
out vec4 o;
void main() {
    vec2 c = gl_PointCoord * 2.0 - 1.0;
    c.y = -c.y;
    float r2 = dot(c, c);
    if (r2 > 1.0) discard;
    vec3 p = vCenter + vec3(c, sqrt(1.0 - r2)) * uRadius;
    vec4 clip = uProj * vec4(p, 1.0);
    gl_FragDepth = clip.z / clip.w * 0.5 + 0.5;
    o = vec4(-p.z, 0.0, 0.0, 1.0);
})";

// 2. Thickness: the chord of every sphere along the ray, added up (additive blending). The sprites
//    overlap, so every chord is scaled by (particle volume / sprite volume): the sum is then the
//    length of water the ray really crosses.
static const char* kThicknessFS = R"(
in vec3 vCenter;
uniform float uRadius, uVolumeScale;
out vec4 o;
void main() {
    vec2 c = gl_PointCoord * 2.0 - 1.0;
    float r2 = dot(c, c);
    if (r2 > 1.0) discard;
    o = vec4(2.0 * uRadius * sqrt(1.0 - r2) * uVolumeScale, 0.0, 0.0, 1.0);
})";

// Full-screen triangle without vertex data.
static const char* kFullScreenVS = R"(
out vec2 uv;
void main() {
    vec2 p = vec2((gl_VertexID << 1) & 2, gl_VertexID & 2);
    uv = p;
    gl_Position = vec4(p * 2.0 - 1.0, 0.0, 1.0);
})";

// 3. Narrow-range filter, one direction (Truong & Yuksel 2018): Gaussian over the neighbours whose
//    depth lies within [z - range, z + range] of the pixel; nearer samples (another surface in
//    front) are dropped, farther ones are clamped to z + range, so silhouettes stay sharp.
static const char* kSmoothFS = R"(
in vec2 uv;
uniform sampler2D uDepth;
uniform vec2 uStep;          // one texel along the filter direction
uniform float uWorldRadius;  // filter radius [m]
uniform float uProjScale;    // pixels per metre at unit distance
uniform float uRange;        // depth range [m]
uniform int uMaxTaps;
out vec4 o;
void main() {
    float zc = texture(uDepth, uv).r;
    if (zc <= 0.0) { o = vec4(0.0); return; }
    float radiusPx = clamp(uWorldRadius * uProjScale / zc, 1.0, float(uMaxTaps));
    float sigma = 0.5 * radiusPx;
    float lo = zc - uRange, hi = zc + uRange;
    float sum = zc, wsum = 1.0;
    for (int i = 1; i <= uMaxTaps; ++i) {
        if (float(i) > radiusPx) break;
        float w = exp(-float(i * i) / (2.0 * sigma * sigma));
        for (int s = -1; s <= 1; s += 2) {
            float z = texture(uDepth, uv + float(s * i) * uStep).r;
            if (z <= 0.0 || z < lo) continue;
            sum += w * min(z, hi);
            wsum += w;
        }
    }
    o = vec4(sum / wsum, 0.0, 0.0, 1.0);
})";

// 4. Shading of the smoothed surface.
static const char* kShadeFS = R"(
in vec2 uv;
uniform sampler2D uDepth, uThickness, uScene;
uniform mat4 uProj, uInvView;
uniform vec2 uTexel;
uniform vec3 uAbsorption, uScatter;
uniform float uRefraction;
uniform float uDropSize; // [m] one particle diameter of water
out vec4 o;
vec3 eyePosition(vec2 t) {
    float z = texture(uDepth, t).r;
    vec2 ndc = t * 2.0 - 1.0;
    return vec3(ndc.x * z / uProj[0][0], ndc.y * z / uProj[1][1], -z);
}
vec3 sky(vec3 d) { // environment seen in the reflection (world direction)
    vec3 horizon = vec3(0.78, 0.84, 0.9), zenith = vec3(0.35, 0.52, 0.78), ground = vec3(0.12, 0.13, 0.15);
    vec3 c = d.y > 0.0 ? mix(horizon, zenith, sqrt(d.y)) : mix(horizon, ground, sqrt(-d.y));
    vec3 sun = normalize(vec3(0.4, 0.8, 0.45));
    return c + vec3(1.0, 0.95, 0.85) * pow(max(dot(d, sun), 0.0), 600.0) * 6.0;
}
void main() {
    float z = texture(uDepth, uv).r;
    if (z <= 0.0) discard;
    vec3 P = eyePosition(uv);
    // Normal from the smoothed depth: of the two one-sided differences take the smaller one, so
    // a neighbour across a silhouette does not tilt it; a neighbour without water (background)
    // does not count at all. A lone drop with no water around faces the viewer.
    vec2 ex = vec2(uTexel.x, 0.0), ey = vec2(0.0, uTexel.y);
    bool right = texture(uDepth, uv + ex).r > 0.0, left = texture(uDepth, uv - ex).r > 0.0;
    bool up = texture(uDepth, uv + ey).r > 0.0, down = texture(uDepth, uv - ey).r > 0.0;
    vec3 n = vec3(0.0, 0.0, 1.0);
    if ((right || left) && (up || down)) {
        vec3 dx = eyePosition(uv + ex) - P, dx2 = P - eyePosition(uv - ex);
        vec3 dy = eyePosition(uv + ey) - P, dy2 = P - eyePosition(uv - ey);
        if (!right || (left && abs(dx2.z) < abs(dx.z))) dx = dx2;
        if (!up || (down && abs(dy2.z) < abs(dy.z))) dy = dy2;
        n = normalize(cross(dx, dy));
    }
    vec3 v = normalize(-P);
    if (dot(n, v) < 0.0) n = -n;

    float thickness = texture(uThickness, uv).r;
    float cosTheta = max(dot(n, v), 0.0);
    float fresnel = 0.02 + 0.98 * pow(1.0 - cosTheta, 5.0); // Schlick, F0 of water (n = 1.33)
    vec3 reflected = sky(normalize(mat3(uInvView) * reflect(-v, n)));
    // What lies behind, bent by the surface, dimmed by Beer-Lambert through the water.
    vec2 bent = clamp(uv + n.xy * uRefraction * min(thickness / 0.1, 1.0), vec2(0.001), vec2(0.999));
    vec3 behind = texture(uScene, bent).rgb;
    vec3 transmit = exp(-uAbsorption * thickness);
    vec3 water = behind * transmit + uScatter * (vec3(1.0) - transmit);
    vec3 l = normalize(vec3(0.35, 0.7, 0.6));
    float spec = pow(max(dot(n, normalize(l + v)), 0.0), 120.0) * 0.8;
    vec3 c = mix(water, reflected, fresnel) + vec3(spec);
    // Spray: where the water is only a few drops thick the light is scattered by the drops (Mie
    // scattering, nearly white) instead of passing through a clear body of water.
    float spray = 1.0 - smoothstep(0.5, 1.5, thickness / uDropSize);
    c = mix(c, vec3(0.82, 0.86, 0.9) * (0.6 + 0.4 * max(dot(n, l), 0.0)), spray * 0.85);
    o = vec4(c, 1.0);
    vec4 clip = uProj * vec4(P, 1.0);
    gl_FragDepth = clip.z / clip.w * 0.5 + 0.5;
})";

// ---------------------------------------------------------------------------

static void build(QOpenGLShaderProgram& p, const QString& header, const char* vs, const char* fs, const char* name) {
    p.bindAttributeLocation("aPos", 0);
    if (!p.addShaderFromSourceCode(QOpenGLShader::Vertex, header + vs) ||
        !p.addShaderFromSourceCode(QOpenGLShader::Fragment, header + fs) || !p.link())
        qWarning("Fluid shader %s failed: %s", name, qPrintable(p.log()));
}

void FluidSurfaceRenderer::initialize(const QString& glslHeader) {
    gl_ = QOpenGLContext::currentContext()->extraFunctions();
    build(depthProg_, glslHeader, kSpriteVS, kDepthFS, "depth");
    build(thickProg_, glslHeader, kSpriteVS, kThicknessFS, "thickness");
    build(smoothProg_, glslHeader, kFullScreenVS, kSmoothFS, "smooth");
    build(shadeProg_, glslHeader, kFullScreenVS, kShadeFS, "shade");
    gl_->glGenVertexArrays(1, &vao_);
    gl_->glGenBuffers(1, &vbo_);
    gl_->glBindVertexArray(vao_);
    gl_->glBindBuffer(GL_ARRAY_BUFFER, vbo_);
    gl_->glEnableVertexAttribArray(0);
    gl_->glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), nullptr);
    gl_->glGenVertexArrays(1, &emptyVao_);
    gl_->glBindVertexArray(0);
}

void FluidSurfaceRenderer::release() {
    if (!gl_) return;
    for (Target* t : {&depth_, &thickness_, &smoothA_, &smoothB_, &scene_}) destroy(*t);
    gl_->glDeleteBuffers(1, &vbo_);
    gl_->glDeleteVertexArrays(1, &vao_);
    gl_->glDeleteVertexArrays(1, &emptyVao_);
    gl_ = nullptr;
}

void FluidSurfaceRenderer::destroy(Target& t) {
    if (t.fbo) gl_->glDeleteFramebuffers(1, &t.fbo);
    if (t.color) gl_->glDeleteTextures(1, &t.color);
    if (t.depth) gl_->glDeleteRenderbuffers(1, &t.depth);
    t = Target();
}

void FluidSurfaceRenderer::resize(Target& t, int w, int h, int internalFormat, int format, int type, bool withDepth) {
    if (t.fbo && t.w == w && t.h == h) return;
    destroy(t);
    t.w = w;
    t.h = h;
    gl_->glGenTextures(1, &t.color);
    gl_->glBindTexture(GL_TEXTURE_2D, t.color);
    gl_->glTexImage2D(GL_TEXTURE_2D, 0, internalFormat, w, h, 0, GLenum(format), GLenum(type), nullptr);
    // Depth must not be blended between texels (a silhouette would get in-between depths).
    const GLint filter = format == GL_RGBA ? GL_LINEAR : GL_NEAREST;
    gl_->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, filter);
    gl_->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, filter);
    gl_->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    gl_->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    gl_->glGenFramebuffers(1, &t.fbo);
    gl_->glBindFramebuffer(GL_FRAMEBUFFER, t.fbo);
    gl_->glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, t.color, 0);
    if (withDepth) {
        gl_->glGenRenderbuffers(1, &t.depth);
        gl_->glBindRenderbuffer(GL_RENDERBUFFER, t.depth);
        gl_->glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, w, h);
        gl_->glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, t.depth);
    }
    if (gl_->glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) qWarning("Fluid render target incomplete");
}

void FluidSurfaceRenderer::drawFullScreen() {
    gl_->glBindVertexArray(emptyVao_);
    gl_->glDrawArrays(GL_TRIANGLES, 0, 3);
}

void FluidSurfaceRenderer::render(const std::vector<rf::Vector3>& particles, float radius, const QMatrix4x4& view,
                                  const QMatrix4x4& proj, int width, int height, unsigned sceneFbo, bool lowQuality) {
    if (!gl_ || particles.empty() || width <= 0 || height <= 0) return;
    for (QOpenGLShaderProgram* p : {&depthProg_, &thickProg_, &smoothProg_, &shadeProg_})
        if (!p->isLinked()) return; // a shader did not build (see the warning): no water rather than garbage
    const int W = lowQuality ? std::max(1, width / 2) : width, H = lowQuality ? std::max(1, height / 2) : height;
    resize(depth_, W, H, GL_R32F, GL_RED, GL_FLOAT, true);
    resize(thickness_, W, H, GL_R32F, GL_RED, GL_FLOAT, false);
    resize(smoothA_, W, H, GL_R32F, GL_RED, GL_FLOAT, false);
    resize(smoothB_, W, H, GL_R32F, GL_RED, GL_FLOAT, false);
    resize(scene_, width, height, GL_RGBA8, GL_RGBA, GL_UNSIGNED_BYTE, false);

    // The opaque scene as a texture, to be seen through the water.
    gl_->glBindFramebuffer(GL_READ_FRAMEBUFFER, sceneFbo);
    gl_->glBindFramebuffer(GL_DRAW_FRAMEBUFFER, scene_.fbo);
    gl_->glBlitFramebuffer(0, 0, width, height, 0, 0, width, height, GL_COLOR_BUFFER_BIT, GL_NEAREST);

    gl_->glBindBuffer(GL_ARRAY_BUFFER, vbo_);
    gl_->glBufferData(GL_ARRAY_BUFFER, GLsizeiptr(particles.size() * sizeof(rf::Vector3)), particles.data(), GL_STREAM_DRAW);
    const float sprite = radius * settings.sphereScale;
    const float projScale = proj(1, 1) * 0.5f * float(H); // pixels per metre at unit distance
    auto sprites = [&](QOpenGLShaderProgram& p) {
        p.bind();
        p.setUniformValue("uView", view);
        p.setUniformValue("uProj", proj);
        p.setUniformValue("uRadius", sprite);
        p.setUniformValue("uPointScale", projScale);
        gl_->glBindVertexArray(vao_);
        gl_->glDrawArrays(GL_POINTS, 0, GLsizei(particles.size()));
    };
    gl_->glViewport(0, 0, W, H);
    gl_->glDisable(GL_BLEND);
    gl_->glClearColor(0, 0, 0, 0);

    // 1. depth
    gl_->glBindFramebuffer(GL_FRAMEBUFFER, depth_.fbo);
    gl_->glEnable(GL_DEPTH_TEST);
    gl_->glDepthMask(GL_TRUE);
    gl_->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    sprites(depthProg_);

    // 2. thickness
    gl_->glBindFramebuffer(GL_FRAMEBUFFER, thickness_.fbo);
    gl_->glDisable(GL_DEPTH_TEST);
    gl_->glClear(GL_COLOR_BUFFER_BIT);
    gl_->glEnable(GL_BLEND);
    gl_->glBlendFunc(GL_ONE, GL_ONE);
    thickProg_.bind();
    const float cell = 2.0f * radius; // particle spacing: every particle stands for a cube (2r)^3
    thickProg_.setUniformValue("uVolumeScale", cell * cell * cell / (4.18879f * sprite * sprite * sprite));
    sprites(thickProg_);
    gl_->glDisable(GL_BLEND);

    // 3. smoothing: (horizontal, vertical) x 2
    smoothProg_.bind();
    smoothProg_.setUniformValue("uDepth", 0);
    smoothProg_.setUniformValue("uWorldRadius", settings.smoothing * radius);
    smoothProg_.setUniformValue("uProjScale", projScale);
    smoothProg_.setUniformValue("uRange", 2.0f * sprite);
    smoothProg_.setUniformValue("uMaxTaps", lowQuality ? 8 : 24);
    gl_->glActiveTexture(GL_TEXTURE0);
    unsigned source = depth_.color;
    for (int pass = 0; pass < 4; ++pass) {
        Target& out = pass % 2 == 0 ? smoothA_ : smoothB_;
        gl_->glBindFramebuffer(GL_FRAMEBUFFER, out.fbo);
        gl_->glBindTexture(GL_TEXTURE_2D, source);
        smoothProg_.setUniformValue("uStep", pass % 2 == 0 ? QVector2D(1.0f / W, 0.0f) : QVector2D(0.0f, 1.0f / H));
        drawFullScreen();
        source = out.color;
    }

    // 4. shading into the scene, depth-tested against it
    gl_->glBindFramebuffer(GL_FRAMEBUFFER, sceneFbo);
    gl_->glViewport(0, 0, width, height);
    gl_->glEnable(GL_DEPTH_TEST);
    gl_->glDepthFunc(GL_LESS);
    shadeProg_.bind();
    shadeProg_.setUniformValue("uDepth", 0);
    shadeProg_.setUniformValue("uThickness", 1);
    shadeProg_.setUniformValue("uScene", 2);
    shadeProg_.setUniformValue("uProj", proj);
    shadeProg_.setUniformValue("uInvView", view.inverted());
    shadeProg_.setUniformValue("uTexel", QVector2D(1.0f / W, 1.0f / H));
    shadeProg_.setUniformValue("uAbsorption", settings.absorption);
    shadeProg_.setUniformValue("uScatter", settings.scatterColor);
    shadeProg_.setUniformValue("uRefraction", settings.refraction);
    shadeProg_.setUniformValue("uDropSize", 2.0f * radius);
    gl_->glActiveTexture(GL_TEXTURE0);
    gl_->glBindTexture(GL_TEXTURE_2D, smoothB_.color);
    gl_->glActiveTexture(GL_TEXTURE1);
    gl_->glBindTexture(GL_TEXTURE_2D, thickness_.color);
    gl_->glActiveTexture(GL_TEXTURE2);
    gl_->glBindTexture(GL_TEXTURE_2D, scene_.color);
    drawFullScreen();
    gl_->glActiveTexture(GL_TEXTURE0);
    gl_->glBindVertexArray(0);
}
