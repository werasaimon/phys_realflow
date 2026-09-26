#include "Viewport.h"

#include <QColor>
#include <QDateTime>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QPainter>
#include <QPainterPath>
#include <QTimer>
#include <QWheelEvent>

#include <cmath>

using rf::Vector3;

// ---------------------------------------------------------------------------
// Shaders
// ---------------------------------------------------------------------------
static const char* kColormapGLSL = R"(
vec3 turbo(float x) {
    const vec4 kR4 = vec4(0.13572138, 4.61539260, -42.66032258, 132.13108234);
    const vec4 kG4 = vec4(0.09140261, 2.19418839, 4.84296658, -14.18503333);
    const vec4 kB4 = vec4(0.10667330, 12.64194608, -60.58204836, 110.36276771);
    const vec2 kR2 = vec2(-152.94239396, 59.28637943);
    const vec2 kG2 = vec2(4.27729857, 2.82956604);
    const vec2 kB2 = vec2(-89.90310912, 27.34824973);
    x = clamp(x, 0.0, 1.0);
    vec4 v4 = vec4(1.0, x, x * x, x * x * x);
    vec2 v2 = v4.zw * v4.z;
    return vec3(dot(v4, kR4) + dot(v2, kR2), dot(v4, kG4) + dot(v2, kG2), dot(v4, kB4) + dot(v2, kB2));
}
vec3 viridis(float t) {
    t = clamp(t, 0.0, 1.0);
    const vec3 c0 = vec3(0.2777273272234177, 0.005407344544966578, 0.3340998053353061);
    const vec3 c1 = vec3(0.1050930431085774, 1.404613529898575, 1.384590162594685);
    const vec3 c2 = vec3(-0.3308618287255563, 0.214847559468213, 0.09509516302823659);
    const vec3 c3 = vec3(-4.634230498983486, -5.799100973351585, -19.33244095627987);
    const vec3 c4 = vec3(6.228269936347081, 14.17993336680509, 56.69055260068105);
    const vec3 c5 = vec3(4.776384997670288, -13.74514537774601, -65.35303263337234);
    const vec3 c6 = vec3(-5.435455855934631, 4.645852612178535, 26.3124352495832);
    return c0 + t * (c1 + t * (c2 + t * (c3 + t * (c4 + t * (c5 + t * c6)))));
}
vec3 coolwarm(float t) {
    t = clamp(t, 0.0, 1.0);
    vec3 a = vec3(0.230, 0.299, 0.754), m = vec3(0.865, 0.865, 0.865), b = vec3(0.706, 0.016, 0.150);
    return t < 0.5 ? mix(a, m, t * 2.0) : mix(m, b, t * 2.0 - 1.0);
}
vec3 cmap(int m, float t) {
    if (m == 1) return viridis(t);
    if (m == 2) return coolwarm(t);
    if (m == 3) return vec3(clamp(t, 0.0, 1.0));
    return turbo(t);
}
)";

static const char* kBgVS = R"(out vec2 uv;
void main() {
    vec2 p = vec2((gl_VertexID << 1) & 2, gl_VertexID & 2);
    uv = p;
    gl_Position = vec4(p * 2.0 - 1.0, 0.9999, 1.0);
})";
static const char* kBgFS = R"(in vec2 uv; out vec4 o;
void main() {
    vec3 top = vec3(0.17, 0.19, 0.23), bot = vec3(0.06, 0.07, 0.09);
    o = vec4(mix(bot, top, uv.y), 1.0);
})";

static const char* kMeshVS = R"(in vec3 aPos;
in vec3 aNormal;
in float aScalar;
uniform mat4 uModel, uView, uProj;
out vec3 vN; out vec3 vPosV; out float vS;
void main() {
    vec4 vp = uView * uModel * vec4(aPos, 1.0);
    vPosV = vp.xyz;
    vN = mat3(uView) * mat3(uModel) * aNormal;
    vS = aScalar;
    gl_Position = uProj * vp;
})";
static const char* kMeshFS = R"(
in vec3 vN; in vec3 vPosV; in float vS;
uniform vec3 uColor; uniform int uUseScalar; uniform float uMin, uMax; uniform int uCmap;
out vec4 o;
void main() {
    vec3 n = normalize(vN);
    vec3 v = normalize(-vPosV);
    if (dot(n, v) < 0.0) n = -n;
    vec3 base = uUseScalar == 1 ? cmap(uCmap, (vS - uMin) / max(uMax - uMin, 1e-6)) : uColor;
    // Burning fabric (uUseScalar 2, vS = how far it has burnt, 0..1): it browns and chars black,
    // and the zone that is burning right now - the pyrolysis front - glows orange.
    vec3 ember = vec3(0.0);
    if (uUseScalar == 2) {
        base = mix(uColor, vec3(0.35, 0.22, 0.12), smoothstep(0.0, 0.25, vS));
        base = mix(base, vec3(0.05, 0.045, 0.04), smoothstep(0.3, 0.9, vS));
        ember = vec3(1.0, 0.42, 0.08) * 1.6 * smoothstep(0.02, 0.2, vS) * (1.0 - smoothstep(0.75, 1.0, vS));
    }
    vec3 l = normalize(vec3(0.35, 0.7, 0.6));
    float diff = max(dot(n, l), 0.0);
    float fill = max(dot(n, v), 0.0);
    vec3 h = normalize(l + v);
    float spec = pow(max(dot(n, h), 0.0), 48.0) * 0.3;
    vec3 c = base * (0.18 + 0.55 * diff + 0.35 * fill) + vec3(spec) + ember;
    o = vec4(c, 1.0);
})";

static const char* kSphereVS = R"(in vec3 aPos;
in float aScalar;
in vec3 aColor;
in float aRadius;
uniform mat4 uView, uProj;
uniform float uPointScale, uRadiusScale;
out vec3 vCenter; out float vRad; out float vS; out vec3 vCol;
void main() {
    vec4 vp = uView * vec4(aPos, 1.0);
    vCenter = vp.xyz;
    vRad = aRadius * uRadiusScale;
    vS = aScalar;
    vCol = aColor;
    gl_Position = uProj * vp;
    gl_PointSize = max(2.0, 2.0 * vRad * uPointScale / max(-vp.z, 1e-3));
})";
static const char* kSphereFS = R"(
in vec3 vCenter; in float vRad; in float vS; in vec3 vCol;
uniform mat4 uProj; uniform int uUseScalar; uniform float uMin, uMax; uniform int uCmap;
out vec4 o;
void main() {
    vec2 c = gl_PointCoord * 2.0 - 1.0;
    c.y = -c.y;
    float r2 = dot(c, c);
    if (r2 > 1.0) discard;
    vec3 n = vec3(c, sqrt(1.0 - r2));
    vec3 pv = vCenter + n * vRad;
    vec4 clip = uProj * vec4(pv, 1.0);
    gl_FragDepth = clip.z / clip.w * 0.5 + 0.5;
    vec3 base = uUseScalar == 1 && vS > -1e8 ? cmap(uCmap, (vS - uMin) / max(uMax - uMin, 1e-6)) : vCol;
    vec3 l = normalize(vec3(0.35, 0.7, 0.6));
    float diff = max(dot(n, l), 0.0);
    vec3 h = normalize(l + vec3(0, 0, 1));
    float spec = pow(max(dot(n, h), 0.0), 40.0) * 0.35;
    o = vec4(base * (0.25 + 0.5 * diff + 0.3 * n.z) + vec3(spec), 1.0);
})";

// A magnetised sphere drawn as an Earth-like planet (procedural - no image files): oceans and
// continents from fractal noise on the sphere, ice caps, drifting clouds, the day side facing the
// Sun (upstream of the plasma wind), an atmosphere rim, sun glint on the oceans, and aurora ovals
// around the magnetic poles (the dipole axis is +y) at ~67 deg magnetic latitude - where the
// field lines that reach the magnetosphere's boundary come down to the surface.
static const char* kPlanetVS = R"(in vec3 aPos;
in vec3 aNormal;
in float aScalar;
uniform mat4 uView, uProj;
out vec3 vWorld;
void main() { vWorld = aPos; gl_Position = uProj * uView * vec4(aPos, 1.0); })";
static const char* kPlanetFS = R"(
in vec3 vWorld;
uniform vec3 uCenter, uEye, uSun;
uniform float uTime;
out vec4 o;
float hash(vec3 p) {
    p = fract(p * 0.3183099 + vec3(0.1, 0.2, 0.3));
    p *= 17.0;
    return fract(p.x * p.y * p.z * (p.x + p.y + p.z));
}
float noise(vec3 x) {
    vec3 i = floor(x), f = fract(x);
    f = f * f * (3.0 - 2.0 * f);
    return mix(mix(mix(hash(i), hash(i + vec3(1, 0, 0)), f.x), mix(hash(i + vec3(0, 1, 0)), hash(i + vec3(1, 1, 0)), f.x), f.y),
               mix(mix(hash(i + vec3(0, 0, 1)), hash(i + vec3(1, 0, 1)), f.x), mix(hash(i + vec3(0, 1, 1)), hash(i + vec3(1, 1, 1)), f.x), f.y),
               f.z);
}
float fbm(vec3 p) {
    float s = 0.0, a = 0.5;
    for (int i = 0; i < 6; ++i) { s += a * noise(p); p = p * 2.03 + vec3(1.7, 9.2, 3.1); a *= 0.5; }
    return s;
}
void main() {
    vec3 n = normalize(vWorld - uCenter);           // direction on the planet = its normal
    vec3 v = normalize(uEye - vWorld);
    float lat = n.y;                                  // sine of the (magnetic) latitude
    // Surface: land where the fractal height is above sea level, ice towards the poles.
    float h = fbm(n * 2.2) - 0.52;
    float ice = smoothstep(0.80, 0.86, abs(lat) + 0.08 * (fbm(n * 6.0) - 0.5));
    vec3 ocean = mix(vec3(0.02, 0.07, 0.20), vec3(0.04, 0.22, 0.38), smoothstep(-0.12, 0.0, h));
    float dry = smoothstep(0.35, 0.65, fbm(n * 3.0 + vec3(5.0)) + 0.3 * (1.0 - abs(lat)));
    vec3 land = mix(vec3(0.13, 0.30, 0.10), vec3(0.55, 0.45, 0.28), dry) * (0.8 + 1.5 * clamp(h, 0.0, 0.2));
    bool isLand = h > 0.0;
    vec3 albedo = mix(isLand ? land : ocean, vec3(0.92, 0.95, 1.0), ice);
    // Clouds drifting eastwards.
    float c = uTime * 0.05;
    vec3 nc = vec3(n.x * cos(c) - n.z * sin(c), n.y, n.x * sin(c) + n.z * cos(c));
    float cloud = smoothstep(0.52, 0.75, fbm(nc * 4.0 + vec3(0.0, 3.0, 0.0)));
    albedo = mix(albedo, vec3(1.0), cloud * 0.85);
    // Sunlight: day and night, glint on open water, atmosphere at the limb.
    float day = max(dot(n, uSun), 0.0);
    float twilight = smoothstep(-0.1, 0.15, dot(n, uSun));
    vec3 col = albedo * (0.03 + 1.1 * day);
    if (!isLand && ice < 0.5) col += vec3(1.0, 0.95, 0.85) * pow(max(dot(n, normalize(uSun + v)), 0.0), 80.0) * 0.6 * (1.0 - cloud) * twilight;
    float rim = pow(1.0 - max(dot(n, v), 0.0), 3.0);
    col += vec3(0.30, 0.55, 1.0) * rim * (0.15 + 0.85 * twilight);
    // Aurora ovals around both magnetic poles, shimmering.
    float oval = exp(-pow((abs(lat) - 0.92) / 0.025, 2.0));
    float curtain = 0.5 + 0.5 * sin(atan(n.z, n.x) * 14.0 + uTime * 1.3 + 6.0 * noise(n * 8.0 + vec3(uTime * 0.3)));
    col += vec3(0.25, 1.0, 0.45) * oval * curtain * (1.0 - 0.7 * day) * 1.4;
    o = vec4(col, 1.0);
})";

static const char* kLineVS = R"(in vec3 aPos;
in float aScalar;
uniform mat4 uVP;
out float vS;
void main() { vS = aScalar; gl_Position = uVP * vec4(aPos, 1.0); })";
static const char* kLineFS = R"(
in float vS;
uniform vec4 uColor; uniform int uUseScalar; uniform float uMin, uMax; uniform int uCmap;
out vec4 o;
void main() {
    o = uUseScalar == 1 ? vec4(cmap(uCmap, (vS - uMin) / max(uMax - uMin, 1e-6)), uColor.a) : uColor;
})";

static const char* kSliceVS = R"(in vec3 aPos;
in vec2 aUV;
uniform mat4 uVP;
out vec2 vUV;
void main() { vUV = aUV; gl_Position = uVP * vec4(aPos, 1.0); })";
static const char* kSliceFS = R"(
in vec2 vUV;
uniform sampler2D uTex; uniform float uMin, uMax, uAlpha; uniform int uCmap;
out vec4 o;
void main() {
    vec2 t = texture(uTex, vUV).rg;
    if (t.g > 0.5) { o = vec4(0.22, 0.23, 0.26, 1.0); return; }
    o = vec4(cmap(uCmap, (t.r - uMin) / max(uMax - uMin, 1e-6)), uAlpha);
})";

static const char* kVolVS = R"(in vec3 aPos;
uniform mat4 uVP; uniform vec3 uLo, uSize;
out vec3 vWorld;
void main() { vWorld = uLo + aPos * uSize; gl_Position = uVP * vec4(vWorld, 1.0); })";
// Ray marching through the gas volume: smoke scatters light (self-shadowed); in fire mode the hot
// gas also emits - blackbody colour of its temperature, brightness ~ T^4 (Stefan-Boltzmann) from
// the Draper point (~800 K, the onset of visible glow) up - and the soot darkens what lies behind.
static const char* kVolFS = R"(in vec3 vWorld;
uniform sampler3D uVol; uniform vec3 uEye, uLo, uSize; uniform float uDensity;
uniform int uSteps;          // samples along the ray (fewer on the CPU renderer)
uniform int uFire;           // 1: channel g is the temperature
uniform float uTempScale;    // kelvin above ambient for g = 1
uniform float uAmbient;      // [K]
uniform float uFlame;        // brightness of the flame
uniform sampler2D uDepth; uniform int uUseDepth; uniform mat4 uInvVP; uniform vec2 uViewport; // scene depth
out vec4 o;
// Colour of a black body at temperature T [K], normalised (T. Helland's fit of the Planck locus).
vec3 blackbody(float T) {
    float t = T / 100.0;
    vec3 c;
    c.r = t <= 66.0 ? 1.0 : 1.292936 * pow(t - 60.0, -0.1332047);
    c.g = t <= 66.0 ? 0.3900816 * log(t) - 0.6318414 : 1.1298909 * pow(t - 60.0, -0.0755148);
    c.b = t >= 66.0 ? 1.0 : (t <= 19.0 ? 0.0 : 0.5432068 * log(t - 10.0) - 1.1962541);
    return clamp(c, 0.0, 1.0);
}
void main() {
    vec3 dir = normalize(vWorld - uEye);
    vec3 inv = 1.0 / (dir + vec3(1e-7));
    vec3 t0 = (uLo - uEye) * inv, t1 = (uLo + uSize - uEye) * inv;
    vec3 tmin = min(t0, t1), tmax = max(t0, t1);
    float tn = max(max(max(tmin.x, tmin.y), tmin.z), 0.0);
    float tf = min(min(tmax.x, tmax.y), tmax.z);
    if (uUseDepth == 1) { // stop at the first solid surface of the scene
        vec2 uv = gl_FragCoord.xy / uViewport;
        float z = texture(uDepth, uv).r;
        if (z < 1.0) {
            vec4 w = uInvVP * vec4(uv * 2.0 - 1.0, z * 2.0 - 1.0, 1.0);
            tf = min(tf, dot(w.xyz / w.w - uEye, dir));
        }
    }
    if (tf <= tn) discard;
    float dt = (tf - tn) / float(uSteps);
    vec4 acc = vec4(0.0);
    vec3 glow = vec3(0.0); // emitted light that reaches the eye
    for (int i = 0; i < uSteps; ++i) {
        float t = tn + (float(i) + 0.5) * dt;
        vec3 p = (uEye + dir * t - uLo) / uSize;
        vec2 s = texture(uVol, p).rg;
        float d = s.r;
        if (uFire == 1) {
            float T = uAmbient + s.g * uTempScale;
            if (T > 800.0) {
                float I = uFlame * pow(T / 1500.0, 4.0) * smoothstep(800.0, 1100.0, T);
                glow += (1.0 - acc.a) * blackbody(T) * I * dt;
            }
        }
        if (uFire == 2 && d > 0.004) {
            // Plasma: an optically thin glowing gas - it emits (line radiation of ionised gas,
            // violet-pink as argon / hydrogen discharges) and hardly absorbs.
            glow += vec3(0.62, 0.32, 1.0) * (d * uFlame * dt);
            continue;
        }
        if (d > 0.004) {
            float a = 1.0 - exp(-d * uDensity * dt);
            vec3 col = mix(vec3(0.55, 0.6, 0.68), vec3(0.97, 0.97, 1.0), clamp(d * 1.5, 0.0, 1.0));
            if (uFire == 1) col *= 0.3; // soot: dark grey
            // Self-shadowing: optical depth towards the light over a few samples; the inside of a
            // dense plume darkens and the volume reads as a volume instead of a flat white blob.
            const vec3 L = vec3(0.37, 0.86, 0.35);
            const float ls = 0.04;
            float od = 0.0;
            for (int k = 1; k <= 6; ++k) od += texture(uVol, p + L * (float(k) * ls) / uSize).r;
            col *= 0.3 + 0.7 * exp(-od * uDensity * ls);
            acc.rgb += (1.0 - acc.a) * a * col;
            acc.a += (1.0 - acc.a) * a;
            if (acc.a > 0.98) break;
        }
    }
    // The flame's light is added over what lies behind (premultiplied blending), tone-mapped so
    // that the hot core saturates to white-yellow instead of clipping.
    o = vec4(acc.rgb + (vec3(1.0) - exp(-glow)), acc.a);
})";


// ---------------------------------------------------------------------------
// Colormaps (CPU side, for the legend)
// ---------------------------------------------------------------------------
QColor Viewport::colormapColor(int map, float t) {
    t = std::clamp(t, 0.0f, 1.0f);
    float r, g, b;
    if (map == Viridis) {
        const float c[7][3] = {{0.2777273f, 0.0054073f, 0.3340998f}, {0.1050930f, 1.4046135f, 1.3845902f},
                               {-0.3308618f, 0.2148476f, 0.0950952f}, {-4.6342305f, -5.7991010f, -19.3324410f},
                               {6.2282699f, 14.1799334f, 56.6905526f}, {4.7763850f, -13.7451454f, -65.3530326f},
                               {-5.4354559f, 4.6458526f, 26.3124352f}};
        float v[3];
        for (int k = 0; k < 3; ++k) {
            float acc = c[6][k];
            for (int i = 5; i >= 0; --i) acc = c[i][k] + t * acc;
            v[k] = acc;
        }
        r = v[0]; g = v[1]; b = v[2];
    } else if (map == CoolWarm) {
        float a[3] = {0.230f, 0.299f, 0.754f}, m[3] = {0.865f, 0.865f, 0.865f}, bb[3] = {0.706f, 0.016f, 0.150f};
        float v[3];
        for (int k = 0; k < 3; ++k) v[k] = t < 0.5f ? a[k] + (m[k] - a[k]) * t * 2 : m[k] + (bb[k] - m[k]) * (t * 2 - 1);
        r = v[0]; g = v[1]; b = v[2];
    } else if (map == Gray) {
        r = g = b = t;
    } else {
        float x = t, x2 = x * x, x3 = x2 * x, x4 = x2 * x2, x5 = x4 * x;
        r = 0.13572138f + 4.61539260f * x - 42.66032258f * x2 + 132.13108234f * x3 - 152.94239396f * x4 + 59.28637943f * x5;
        g = 0.09140261f + 2.19418839f * x + 4.84296658f * x2 - 14.18503333f * x3 + 4.27729857f * x4 + 2.82956604f * x5;
        b = 0.10667330f + 12.64194608f * x - 60.58204836f * x2 + 110.36276771f * x3 - 89.90310912f * x4 + 27.34824973f * x5;
    }
    auto cl = [](float v) { return int(std::clamp(v, 0.0f, 1.0f) * 255.0f + 0.5f); };
    return QColor(cl(r), cl(g), cl(b));
}

// ---------------------------------------------------------------------------
// Setup
// ---------------------------------------------------------------------------
Viewport::Viewport(QWidget* parent) : QOpenGLWidget(parent) {
    setFocusPolicy(Qt::StrongFocus);
    setMinimumSize(320, 240);
}

Viewport::~Viewport() {
    makeCurrent();
    if (sliceTex_) glDeleteTextures(1, &sliceTex_);
    if (volumeTex_) glDeleteTextures(1, &volumeTex_);
    if (depthTex_) glDeleteTextures(1, &depthTex_);
    if (depthFbo_) glDeleteFramebuffers(1, &depthFbo_);
    hullCache_.clear(); // GL objects must die while the context is current
    fluidSurface_.release();
    doneCurrent();
}

// Vertex attributes are bound by name (attribute k = attributes[k]) rather than with
// layout(location), which GLSL 1.30 does not have.
static void buildProgram(QOpenGLShaderProgram& p, const QString& vs, const QString& fs, const char* name,
                         std::initializer_list<const char*> attributes) {
    int location = 0;
    for (const char* a : attributes) p.bindAttributeLocation(a, location++);
    if (!p.addShaderFromSourceCode(QOpenGLShader::Vertex, vs) || !p.addShaderFromSourceCode(QOpenGLShader::Fragment, fs) ||
        !p.link())
        qWarning("Shader %s failed: %s", name, qPrintable(p.log()));
}

void Viewport::buildShaders() {
    const QString h = glslHeader_;
    const QString hc = h + kColormapGLSL; // fragment shaders that colour by a field
    buildProgram(bgProg_, h + kBgVS, h + kBgFS, "background", {});
    buildProgram(meshProg_, h + kMeshVS, hc + kMeshFS, "mesh", {"aPos", "aNormal", "aScalar"});
    buildProgram(sphereProg_, h + kSphereVS, hc + kSphereFS, "sphere", {"aPos", "aScalar", "aColor", "aRadius"});
    buildProgram(lineProg_, h + kLineVS, hc + kLineFS, "line", {"aPos", "aScalar"});
    buildProgram(sliceProg_, h + kSliceVS, hc + kSliceFS, "slice", {"aPos", "aUV"});
    buildProgram(volumeProg_, h + kVolVS, h + kVolFS, "volume", {"aPos"});
    buildProgram(planetProg_, h + kPlanetVS, h + kPlanetFS, "planet", {"aPos", "aNormal", "aScalar"});
}

void Viewport::initializeGL() {
    initializeOpenGLFunctions();
    renderer_ = QString::fromLatin1(reinterpret_cast<const char*>(glGetString(GL_RENDERER)));
    const QSurfaceFormat f = context()->format();
    const int version = f.majorVersion() * 10 + f.minorVersion();
    if (version < 30) {
        glOk_ = false; // nothing below would work: say so, offer the software renderer
        QTimer::singleShot(0, this, [this] { emit openGLUnsupported(renderer_); });
        return;
    }
    // Qt gives a 3.3 core context on normal drivers; old drivers and Qt's llvmpipe give 3.0.
    const bool core = version >= 33 && f.profile() == QSurfaceFormat::CoreProfile;
    glslHeader_ = core ? "#version 330 core\n" : "#version 130\n";
    buildShaders();
    fluidSurface_.initialize(glslHeader_);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_PROGRAM_POINT_SIZE);
    glEnable(GL_MULTISAMPLE);
    if (!core) glEnable(0x8861); // GL_POINT_SPRITE: gl_PointCoord needs it outside core profiles

    bg_.vao.create();

    // Unit cube with face normals (for boxes), positions in [-1, 1].
    {
        std::vector<float> v;
        const int faces[6][4][3] = {
            {{-1, -1, -1}, {-1, -1, 1}, {-1, 1, 1}, {-1, 1, -1}}, {{1, -1, -1}, {1, 1, -1}, {1, 1, 1}, {1, -1, 1}},
            {{-1, -1, -1}, {1, -1, -1}, {1, -1, 1}, {-1, -1, 1}}, {{-1, 1, -1}, {-1, 1, 1}, {1, 1, 1}, {1, 1, -1}},
            {{-1, -1, -1}, {-1, 1, -1}, {1, 1, -1}, {1, -1, -1}}, {{-1, -1, 1}, {1, -1, 1}, {1, 1, 1}, {-1, 1, 1}}};
        const float normals[6][3] = {{-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1}};
        const int order[6] = {0, 1, 2, 0, 2, 3};
        for (int f = 0; f < 6; ++f)
            for (int k : order) {
                for (int c = 0; c < 3; ++c) v.push_back(float(faces[f][k][c]));
                for (int c = 0; c < 3; ++c) v.push_back(normals[f][c]);
                v.push_back(0.0f);
            }
        cube_.vao.create();
        cube_.vao.bind();
        cube_.vbo.create();
        cube_.vbo.bind();
        cube_.vbo.allocate(v.data(), int(v.size() * sizeof(float)));
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 7 * sizeof(float), nullptr);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 7 * sizeof(float), reinterpret_cast<void*>(3 * sizeof(float)));
        glEnableVertexAttribArray(2);
        glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, 7 * sizeof(float), reinterpret_cast<void*>(6 * sizeof(float)));
        cube_.count = 36;
        cube_.vao.release();
    }
    // Unit cube [0,1]^3 as triangles for the volume pass.
    {
        std::vector<float> v;
        const int faces[6][4][3] = {
            {{0, 0, 0}, {0, 0, 1}, {0, 1, 1}, {0, 1, 0}}, {{1, 0, 0}, {1, 1, 0}, {1, 1, 1}, {1, 0, 1}},
            {{0, 0, 0}, {1, 0, 0}, {1, 0, 1}, {0, 0, 1}}, {{0, 1, 0}, {0, 1, 1}, {1, 1, 1}, {1, 1, 0}},
            {{0, 0, 0}, {0, 1, 0}, {1, 1, 0}, {1, 0, 0}}, {{0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}}};
        const int order[6] = {0, 1, 2, 0, 2, 3};
        for (auto& f : faces)
            for (int k : order)
                for (int c = 0; c < 3; ++c) v.push_back(float(f[k][c]));
        volumeBox_.vao.create();
        volumeBox_.vao.bind();
        volumeBox_.vbo.create();
        volumeBox_.vbo.bind();
        volumeBox_.vbo.allocate(v.data(), int(v.size() * sizeof(float)));
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), nullptr);
        volumeBox_.count = 36;
        volumeBox_.vao.release();
    }
    // Dynamic buffers.
    auto makeDyn = [&](Buffer& b) {
        b.vao.create();
        b.vbo.create();
        b.vbo.setUsagePattern(QOpenGLBuffer::DynamicDraw);
    };
    makeDyn(obstacle_);
    makeDyn(particles_);
    makeDyn(lines_);
    makeDyn(slice_);
    obstacleScalar_.create();
    obstacleScalar_.setUsagePattern(QOpenGLBuffer::DynamicDraw);

    lines_.vao.bind();
    lines_.vbo.bind();
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 4 * sizeof(float), nullptr);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 4 * sizeof(float), reinterpret_cast<void*>(3 * sizeof(float)));
    lines_.vao.release();

    particles_.vao.bind();
    particles_.vbo.bind();
    const int stride = 8 * sizeof(float);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, nullptr);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(3 * sizeof(float)));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(4 * sizeof(float)));
    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(7 * sizeof(float)));
    particles_.vao.release();

    slice_.vao.bind();
    slice_.vbo.bind();
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), nullptr);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), reinterpret_cast<void*>(3 * sizeof(float)));
    slice_.vao.release();

    glGenTextures(1, &sliceTex_);
    glBindTexture(GL_TEXTURE_2D, sliceTex_);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glGenTextures(1, &volumeTex_);
    glBindTexture(GL_TEXTURE_3D, volumeTex_);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
}

void Viewport::resizeGL(int, int) {}

// ---------------------------------------------------------------------------
// Camera & interaction
// ---------------------------------------------------------------------------
void Viewport::frameScene() {
    if (!snap_ || !snap_->domain.valid()) return;
    // Voxel/vector scenes look straight at the Z slice; others use a 3/4 view.
    if (snap_->vis.gridDisplay > 0 || snap_->vis.vectorDisplay > 0) camera_.frame(snap_->domain, 90.0f, 12.0f);
    else camera_.frame(snap_->domain, snap_->mode == rf::SimMode::WindTunnel ? 60.0f : 55.0f, 22.0f);
    update();
}

ViewportTool& Viewport::toolFor(Qt::MouseButton button, const QPointF& pos) {
    // LMB is the disturbance ray wherever there is gas to disturb - unless it hits a rigid body,
    // which it then grabs; the camera lives on the wheel.
    bool gas = snap_ && snap_->mode == rf::SimMode::WindTunnel;
    if (button == Qt::LeftButton && gas) {
        int body;
        Vector3 hit;
        const Ray ray = camera_.screenRay(pos, size());
        if (!snap_->bodies.empty() && pickBody(ray, body, hit)) return grabTool_;
        if (pickParticle(ray, hit)) return grabTool_; // cloth / soft body
        return disturbTool_;
    }
    if (button == Qt::LeftButton && snap_ && (!snap_->bodies.empty() || !snap_->cloths.empty() || !snap_->softMeshes.empty()))
        return grabTool_;
    // Liquid / soft bodies / cloth: LMB grabs a particle (empty space still orbits).
    if (button == Qt::LeftButton && snap_ && snap_->mode == rf::SimMode::Fluid) return grabTool_;
    return orbitTool_;
}

void Viewport::mousePressEvent(QMouseEvent* e) {
    if (!grabbed_) grabbed_ = &toolFor(e->button(), e->position());
    grabbed_->press(*this, e);
}

void Viewport::mouseMoveEvent(QMouseEvent* e) {
    if (grabbed_) grabbed_->move(*this, e);
}

void Viewport::mouseReleaseEvent(QMouseEvent* e) {
    if (grabbed_) grabbed_->release(*this, e);
    if (e->buttons() == Qt::NoButton) grabbed_ = nullptr;
}

void Viewport::wheelEvent(QWheelEvent* e) {
    camera_.zoom(e->angleDelta().y() / 120.0f);
    update();
}

void Viewport::mouseDoubleClickEvent(QMouseEvent*) { frameScene(); }

void Viewport::keyPressEvent(QKeyEvent* e) {
    if (e->key() == Qt::Key_F) frameScene();
    else QOpenGLWidget::keyPressEvent(e);
}

bool Viewport::pickPoint(const Ray& ray, Vector3& hit) const {
    if (!snap_ || snap_->mode != rf::SimMode::WindTunnel || !snap_->domain.valid()) return false;
    const rf::AABB& d = snap_->domain;
    if (snap_->gridNx > 0) {
        int axis = std::clamp(snap_->vis.sliceAxis, 0, 2);
        float coord = snap_->gridOrigin[axis] + (snap_->sliceLayer + 0.5f) * snap_->gridDx;
        float t = intersectAxisPlane(ray, axis, coord);
        if (t > 0) {
            Vector3 p = ray.at(t);
            if (d.contains(p)) { hit = p; return true; }
        }
    }
    float tn, tf;
    if (!intersectBox(ray, d, tn, tf)) return false;
    hit = ray.at(0.5f * (tn + tf));
    return true;
}

bool Viewport::pickParticle(const Ray& ray, Vector3& hit) const {
    if (!snap_) return false;
    const float r = snap_->particleRadius * 1.5f; // a little generous: particles are small
    float best = rf::kInf;
    auto test = [&](const Vector3& c) {
        Vector3 oc = c - ray.origin;
        float t = rf::dot(oc, ray.dir);
        if (t <= 0 || t >= best) return;
        if (rf::length2(oc - ray.dir * t) > r * r) return;
        best = t;
        hit = c;
    };
    for (const Vector3& p : snap_->particles) test(p);
    for (const auto& c : snap_->cloths)
        for (const Vector3& p : c.positions) test(p);
    for (const auto& m : snap_->softMeshes)
        for (const Vector3& p : m.positions) test(p);
    return best < rf::kInf;
}

bool Viewport::pickBody(const Ray& ray, int& body, Vector3& hit) const {
    body = -1;
    if (!snap_) return false;
    float best = 1e9f;
    for (int i = 0; i < int(snap_->bodies.size()); ++i) {
        const auto& b = snap_->bodies[i];
        if (!b.collisionShape) continue;
        rf::Matrix3x3 Rt = b.rot.toMatrix3x3().transposed();
        float t;
        Vector3 n;
        if (b.collisionShape->raycast(Rt * (ray.origin - b.pos), Rt * ray.dir, best, t, n) && t < best) {
            best = t;
            body = i;
        }
    }
    if (body >= 0) {
        hit = ray.at(best);
        return true;
    }
    // Sticky picking: small or distant bodies are only a few pixels wide, so when the exact ray
    // misses, take the body whose centre passes closest to the ray within ~15 pixels on screen.
    const float pixel = 2.0f * std::tan(qDegreesToRadians(camera_.fovDeg() * 0.5f)) / std::max(1, height());
    float bestRatio = 1.0f;
    for (int i = 0; i < int(snap_->bodies.size()); ++i) {
        const auto& b = snap_->bodies[i];
        float t = rf::dot(b.pos - ray.origin, ray.dir);
        if (t <= 0) continue;
        float perp = rf::length(b.pos - ray.at(t));
        float tol = b.radius + 15.0f * pixel * t;
        if (perp / tol < bestRatio) {
            bestRatio = perp / tol;
            body = i;
            Vector3 off = ray.at(t) - b.pos;
            float l = rf::length(off), r = 0.8f * b.radius;
            hit = b.pos + (l > r ? off * (r / l) : off); // anchor inside the body, towards the cursor
        }
    }
    return body >= 0;
}

void Viewport::drawJoints(const QMatrix4x4& vp) {
    if (!snap_ || snap_->joints.empty()) return;
    const float c = 0.05f;
    auto cross3 = [&](std::vector<float>& d, const Vector3& p) {
        for (int k = 0; k < 3; ++k) {
            Vector3 e(0.0f);
            e[k] = c;
            Vector3 p0 = p - e, p1 = p + e;
            d.insert(d.end(), {p0.x, p0.y, p0.z, 0, p1.x, p1.y, p1.z, 0});
        }
    };
    const QVector4D colors[5] = {{1.0f, 0.55f, 0.2f, 1}, {0.3f, 0.8f, 1.0f, 1}, {0.4f, 1.0f, 0.45f, 1},
                                 {0.9f, 0.45f, 1.0f, 1}, {1.0f, 0.9f, 0.35f, 1}};
    glDisable(GL_DEPTH_TEST);
    for (int type = 0; type < 5; ++type) {
        std::vector<float> d;
        for (const auto& j : snap_->joints) {
            if (int(j.type) != type) continue;
            Vector3 a = j.anchorA, b = j.anchorB;
            d.insert(d.end(), {a.x, a.y, a.z, 0, b.x, b.y, b.z, 0});
            cross3(d, b);
            if (j.type == rf::JointType::Hinge || j.type == rf::JointType::Slider) {
                float L = j.type == rf::JointType::Slider ? 1.0f : 0.25f;
                Vector3 p0 = a - j.axis * L, p1 = a + j.axis * L;
                d.insert(d.end(), {p0.x, p0.y, p0.z, 0, p1.x, p1.y, p1.z, 0});
            }
        }
        drawLines(d, GL_LINES, vp, colors[type], 1);
    }
    glEnable(GL_DEPTH_TEST);
}

void Viewport::drawGrab(const QMatrix4x4& vp) {
    if (!snap_ || !snap_->grabActive) return;
    Vector3 a = snap_->grabAnchor, b = snap_->grabTarget;
    std::vector<float> d = {a.x, a.y, a.z, 0, b.x, b.y, b.z, 0};
    float c = 0.03f * std::max(0.3f, rf::length(snap_->domain.extent()) * 0.1f);
    for (const Vector3& p : {a, b})
        for (int k = 0; k < 3; ++k) {
            Vector3 e(0.0f);
            e[k] = c;
            Vector3 p0 = p - e, p1 = p + e;
            d.insert(d.end(), {p0.x, p0.y, p0.z, 0, p1.x, p1.y, p1.z, 0});
        }
    glDisable(GL_DEPTH_TEST);
    drawLines(d, GL_LINES, vp, QVector4D(1.0f, 0.9f, 0.3f, 1.0f), 1);
    glEnable(GL_DEPTH_TEST);
}

void Viewport::showProbe(const Ray& ray, const Vector3& hit) {
    probeRay_ = ray;
    probeHit_ = hit;
    probeTime_ = QDateTime::currentMSecsSinceEpoch();
    update();
}

void Viewport::drawProbe(const QMatrix4x4& vp) {
    if (!snap_ || QDateTime::currentMSecsSinceEpoch() - probeTime_ > 1500) return;
    std::vector<float> d;
    // Ray segment inside the domain.
    float tn, tf;
    if (intersectBox(probeRay_, snap_->domain, tn, tf)) {
        Vector3 a = probeRay_.at(tn), b = probeRay_.at(tf);
        d.insert(d.end(), {a.x, a.y, a.z, 0, b.x, b.y, b.z, 0});
    }
    // Brush circle facing the camera + a small cross at the hit point.
    Vector3 n = probeRay_.dir;
    Vector3 u = rf::normalize(std::fabs(n.y) < 0.9f ? rf::cross(n, Vector3(0, 1, 0)) : rf::cross(n, Vector3(1, 0, 0)));
    Vector3 w = rf::cross(n, u);
    const int seg = 40;
    for (int i = 0; i < seg; ++i) {
        float t0 = 2 * rf::kPi * i / seg, t1 = 2 * rf::kPi * (i + 1) / seg;
        Vector3 p0 = probeHit_ + (u * std::cos(t0) + w * std::sin(t0)) * brushRadius_;
        Vector3 p1 = probeHit_ + (u * std::cos(t1) + w * std::sin(t1)) * brushRadius_;
        d.insert(d.end(), {p0.x, p0.y, p0.z, 0, p1.x, p1.y, p1.z, 0});
    }
    float c = brushRadius_ * 0.25f;
    for (int a = 0; a < 3; ++a) {
        Vector3 e(0.0f);
        e[a] = c;
        Vector3 p0 = probeHit_ - e, p1 = probeHit_ + e;
        d.insert(d.end(), {p0.x, p0.y, p0.z, 0, p1.x, p1.y, p1.z, 0});
    }
    glDisable(GL_DEPTH_TEST);
    drawLines(d, GL_LINES, vp, QVector4D(1.0f, 0.85f, 0.2f, 1.0f), 1);
    glEnable(GL_DEPTH_TEST);
    QTimer::singleShot(100, this, [this] { update(); }); // fade out
}

// ---------------------------------------------------------------------------
// Data upload
// ---------------------------------------------------------------------------
void Viewport::setSnapshot(std::shared_ptr<const rf::RenderSnapshot> s) {
    snap_ = std::move(s);
    if (snap_ && (!framedOnce_ || snap_->preset != lastPreset_)) {
        framedOnce_ = true;
        lastPreset_ = snap_->preset;
        frameScene();
    }
    update();
}

void Viewport::uploadObstacle(const rf::TriMesh& m) {
    // Expanded, flat-or-smooth shaded triangles (crease angle 35 degrees).
    auto vn = m.vertexNormals();
    std::vector<float> v;
    v.reserve(m.triangles.size() * 3 * 7);
    obstacleVertexMap_.clear();
    const float creaseCos = std::cos(35.0f * rf::kPi / 180.0f);
    for (size_t t = 0; t < m.triangles.size(); ++t) {
        Vector3 fn = m.faceNormal(t);
        for (int k = 0; k < 3; ++k) {
            uint32_t id = m.triangles[t][k];
            const Vector3& p = m.positions[id];
            Vector3 n = rf::dot(vn[id], fn) > creaseCos ? vn[id] : fn;
            v.insert(v.end(), {p.x, p.y, p.z, n.x, n.y, n.z});
            obstacleVertexMap_.push_back(id);
        }
    }
    obstacle_.count = int(obstacleVertexMap_.size());
    obstacle_.vao.bind();
    obstacle_.vbo.bind();
    obstacle_.vbo.allocate(v.data(), int(v.size() * sizeof(float)));
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), nullptr);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), reinterpret_cast<void*>(3 * sizeof(float)));
    obstacleScalar_.bind();
    std::vector<float> zeros(obstacleVertexMap_.size(), 0.0f);
    obstacleScalar_.allocate(zeros.data(), int(zeros.size() * sizeof(float)));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, sizeof(float), nullptr);
    obstacle_.vao.release();
}

// ---------------------------------------------------------------------------
// Drawing
// ---------------------------------------------------------------------------
void Viewport::drawBackground() {
    glDisable(GL_DEPTH_TEST);
    bgProg_.bind();
    bg_.vao.bind();
    glDrawArrays(GL_TRIANGLES, 0, 3);
    bg_.vao.release();
    glEnable(GL_DEPTH_TEST);
}

void Viewport::drawLines(const std::vector<float>& data, GLenum mode, const QMatrix4x4& vp, const QVector4D& color, float) {
    if (data.empty()) return;
    lineProg_.bind();
    lineProg_.setUniformValue("uVP", vp);
    lineProg_.setUniformValue("uColor", color);
    lineProg_.setUniformValue("uUseScalar", 0);
    lines_.vao.bind();
    lines_.vbo.bind();
    lines_.vbo.allocate(data.data(), int(data.size() * sizeof(float)));
    glDrawArrays(mode, 0, GLsizei(data.size() / 4));
    lines_.vao.release();
}

void Viewport::drawObstacle(const QMatrix4x4& view, const QMatrix4x4& proj) {
    if (!snap_->obstacle || snap_->obstacle->empty()) return;
    if (snap_->obstacleVersion != obstacleVersion_) {
        uploadObstacle(*snap_->obstacle);
        obstacleVersion_ = snap_->obstacleVersion;
    }
    bool useScalar = !snap_->obstacleScalar.empty() && snap_->obstacleScalar.size() == snap_->obstacle->positions.size();
    float lo = -1.0f, hi = 1.0f;
    if (useScalar) {
        std::vector<float> s(obstacleVertexMap_.size());
        float mn = 1e9f, mx = -1e9f;
        for (size_t i = 0; i < s.size(); ++i) {
            s[i] = snap_->obstacleScalar[obstacleVertexMap_[i]];
            mn = std::min(mn, s[i]);
            mx = std::max(mx, s[i]);
        }
        lo = std::max(mn, -3.0f);
        hi = std::min(std::max(mx, lo + 0.1f), 1.2f);
        obstacleScalar_.bind();
        obstacleScalar_.write(0, s.data(), int(s.size() * sizeof(float)));
    }
    if (snap_->vis.planetSurface) { // the magnetised sphere as a planet, lit by the Sun upstream
        const rf::AABB b = snap_->obstacle->bounds();
        const Vector3 c = b.center();
        planetProg_.bind();
        planetProg_.setUniformValue("uView", view);
        planetProg_.setUniformValue("uProj", proj);
        planetProg_.setUniformValue("uCenter", QVector3D(c.x, c.y, c.z));
        planetProg_.setUniformValue("uEye", eyePosition());
        planetProg_.setUniformValue("uSun", QVector3D(-0.92f, 0.15f, 0.36f).normalized());
        planetProg_.setUniformValue("uTime", snap_->time);
        obstacle_.vao.bind();
        glDrawArrays(GL_TRIANGLES, 0, obstacle_.count);
        obstacle_.vao.release();
        return;
    }
    meshProg_.bind();
    meshProg_.setUniformValue("uModel", QMatrix4x4());
    meshProg_.setUniformValue("uView", view);
    meshProg_.setUniformValue("uProj", proj);
    meshProg_.setUniformValue("uColor", QVector3D(0.78f, 0.80f, 0.84f));
    meshProg_.setUniformValue("uUseScalar", useScalar ? 1 : 0);
    meshProg_.setUniformValue("uMin", lo);
    meshProg_.setUniformValue("uMax", hi);
    meshProg_.setUniformValue("uCmap", int(CoolWarm));
    obstacle_.vao.bind();
    glDrawArrays(GL_TRIANGLES, 0, obstacle_.count);
    obstacle_.vao.release();
}

Viewport::GpuMesh& Viewport::hullMesh(const std::shared_ptr<const rf::TriMesh>& m, bool smooth) {
    GpuMesh& g = hullCache_[m.get()];
    if (g.vao) return g;
    std::vector<float> v;
    std::vector<Vector3> vn;
    if (smooth) vn = m->vertexNormals();
    for (size_t t = 0; t < m->triangles.size(); ++t) {
        Vector3 fn = m->faceNormal(t); // flat shading shows the facets GJK/EPA works with
        for (int k = 0; k < 3; ++k) {
            const Vector3& p = m->positions[m->triangles[t][k]];
            const Vector3& n = smooth ? vn[m->triangles[t][k]] : fn;
            v.insert(v.end(), {p.x, p.y, p.z, n.x, n.y, n.z, 0.0f});
        }
    }
    g.vao = std::make_unique<QOpenGLVertexArrayObject>();
    g.vao->create();
    g.vao->bind();
    g.vbo = std::make_unique<QOpenGLBuffer>(QOpenGLBuffer::VertexBuffer);
    g.vbo->create();
    g.vbo->bind();
    g.vbo->allocate(v.data(), int(v.size() * sizeof(float)));
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 7 * sizeof(float), nullptr);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 7 * sizeof(float), reinterpret_cast<void*>(3 * sizeof(float)));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, 7 * sizeof(float), reinterpret_cast<void*>(6 * sizeof(float)));
    g.vao->release();
    g.count = int(v.size() / 7);
    return g;
}

void Viewport::drawBodies(const QMatrix4x4& view, const QMatrix4x4& proj) {
    bool anyHull = false;
    auto meshBody = [](const auto& b) { return (b.shape == rf::ShapeType::ConvexHull || b.shape == rf::ShapeType::Compound) && b.mesh; };
    for (const auto& b : snap_->bodies) anyHull |= meshBody(b);
    if (anyHull) {
        meshProg_.bind();
        meshProg_.setUniformValue("uView", view);
        meshProg_.setUniformValue("uProj", proj);
        meshProg_.setUniformValue("uUseScalar", 0);
        for (const auto& b : snap_->bodies) {
            if (!meshBody(b)) continue;
            rf::Matrix3x3 R = b.rot.toMatrix3x3();
            QMatrix4x4 M(R.m[0][0], R.m[0][1], R.m[0][2], b.pos.x, R.m[1][0], R.m[1][1], R.m[1][2], b.pos.y,
                         R.m[2][0], R.m[2][1], R.m[2][2], b.pos.z, 0, 0, 0, 1);
            meshProg_.setUniformValue("uModel", M);
            float dim = b.sleeping ? 0.55f : 1.0f; // sleeping bodies are drawn darker
            meshProg_.setUniformValue("uColor", QVector3D(b.color.x, b.color.y, b.color.z) * dim);
            if (b.shape == rf::ShapeType::Compound && showConvexParts_ && b.collisionShape) {
                // Each convex part in its own shade: shows what the collision pipeline sees.
                const auto& parts = static_cast<const rf::CompoundShape*>(b.collisionShape.get())->partMeshes();
                for (size_t k = 0; k < parts.size(); ++k) {
                    float h = float(k) * 0.618034f;
                    h -= std::floor(h);
                    QColor c = QColor::fromHsvF(h, 0.55f, 0.95f * dim);
                    meshProg_.setUniformValue("uColor", QVector3D(c.redF(), c.greenF(), c.blueF()));
                    GpuMesh& g = hullMesh(parts[k]);
                    g.vao->bind();
                    glDrawArrays(GL_TRIANGLES, 0, g.count);
                    g.vao->release();
                }
                continue;
            }
            GpuMesh& g = hullMesh(b.mesh, b.shape == rf::ShapeType::Compound);
            g.vao->bind();
            glDrawArrays(GL_TRIANGLES, 0, g.count);
            g.vao->release();
        }
    }
    bool anyBox = false;
    for (const auto& b : snap_->bodies) anyBox |= b.shape == rf::ShapeType::Box;
    if (anyBox) {
        meshProg_.bind();
        meshProg_.setUniformValue("uView", view);
        meshProg_.setUniformValue("uProj", proj);
        meshProg_.setUniformValue("uUseScalar", 0);
        cube_.vao.bind();
        for (const auto& b : snap_->bodies) {
            if (b.shape != rf::ShapeType::Box) continue;
            rf::Matrix3x3 R = b.rot.toMatrix3x3();
            QMatrix4x4 M(R.m[0][0] * b.halfExtents.x, R.m[0][1] * b.halfExtents.y, R.m[0][2] * b.halfExtents.z, b.pos.x,
                         R.m[1][0] * b.halfExtents.x, R.m[1][1] * b.halfExtents.y, R.m[1][2] * b.halfExtents.z, b.pos.y,
                         R.m[2][0] * b.halfExtents.x, R.m[2][1] * b.halfExtents.y, R.m[2][2] * b.halfExtents.z, b.pos.z,
                         0, 0, 0, 1);
            meshProg_.setUniformValue("uModel", M);
            float dim = b.sleeping ? 0.55f : 1.0f;
            meshProg_.setUniformValue("uColor", QVector3D(b.color.x, b.color.y, b.color.z) * dim);
            glDrawArrays(GL_TRIANGLES, 0, cube_.count);
        }
        cube_.vao.release();
    }
}

void Viewport::drawCloths(const QMatrix4x4& view, const QMatrix4x4& proj) {
    if (snap_->cloths.empty() && snap_->softMeshes.empty()) return;
    // Both sides of every sheet (each with its own normal), smooth normals from the grid.
    std::vector<float> v;
    std::vector<std::pair<int, int>> ranges; // first vertex, count per cloth
    for (const auto& c : snap_->cloths) {
        const int w = c.width, h = c.height;
        const auto& P = c.positions;
        std::vector<Vector3> N(P.size(), Vector3(0.0f));
        auto id = [w](int x, int y) { return x + w * y; };
        for (int y = 0; y + 1 < h; ++y)
            for (int x = 0; x + 1 < w; ++x) {
                int a = id(x, y), b = id(x + 1, y), cc = id(x, y + 1), d = id(x + 1, y + 1);
                Vector3 n1 = cross(P[b] - P[a], P[d] - P[a]), n2 = cross(P[d] - P[a], P[cc] - P[a]);
                N[a] += n1 + n2; N[b] += n1; N[d] += n1 + n2; N[cc] += n2;
            }
        for (Vector3& n : N) n = rf::normalize(n);
        const int first = int(v.size() / 7);
        auto vert = [&](int i, float side) {
            const float burnt = c.burnt.empty() ? 0.0f : c.burnt[i];
            v.insert(v.end(), {P[i].x, P[i].y, P[i].z, N[i].x * side, N[i].y * side, N[i].z * side, burnt});
        };
        // Cells a crack has cut through are not drawn: the torn edge runs along the threads.
        auto tri = [&](int p, int q, int r) {
            vert(p, 1); vert(q, 1); vert(r, 1);    // front
            vert(p, -1); vert(r, -1); vert(q, -1); // back
        };
        for (int y = 0; y + 1 < h; ++y)
            for (int x = 0; x + 1 < w; ++x) {
                if (!c.cellIntact.empty() && !c.cellIntact[x + (w - 1) * y]) continue;
                int a = id(x, y), b = id(x + 1, y), cc = id(x, y + 1), d = id(x + 1, y + 1);
                tri(a, b, d);
                tri(a, d, cc);
            }
        ranges.push_back({first, int(v.size() / 7) - first});
    }
    // Soft bodies: their skinned surfaces with smooth normals.
    for (const auto& m : snap_->softMeshes) {
        const auto& P = m.positions;
        std::vector<Vector3> N(P.size(), Vector3(0.0f));
        for (const auto& t : m.triangles) {
            Vector3 n = rf::cross(P[t[1]] - P[t[0]], P[t[2]] - P[t[0]]);
            for (uint32_t k : t) N[k] += n;
        }
        for (Vector3& n : N) n = rf::normalize(n);
        const int first = int(v.size() / 7);
        for (const auto& t : m.triangles)
            for (uint32_t k : t) v.insert(v.end(), {P[k].x, P[k].y, P[k].z, N[k].x, N[k].y, N[k].z, 0.0f});
        ranges.push_back({first, int(v.size() / 7) - first});
    }
    GpuMesh& g = clothMesh_;
    if (!g.vao) {
        g.vao = std::make_unique<QOpenGLVertexArrayObject>();
        g.vao->create();
        g.vao->bind();
        g.vbo = std::make_unique<QOpenGLBuffer>(QOpenGLBuffer::VertexBuffer);
        g.vbo->create();
        g.vbo->setUsagePattern(QOpenGLBuffer::DynamicDraw);
        g.vbo->bind();
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 7 * sizeof(float), nullptr);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 7 * sizeof(float), reinterpret_cast<void*>(3 * sizeof(float)));
        glEnableVertexAttribArray(2);
        glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, 7 * sizeof(float), reinterpret_cast<void*>(6 * sizeof(float)));
        g.vao->release();
    }
    g.vbo->bind();
    g.vbo->allocate(v.data(), int(v.size() * sizeof(float)));
    meshProg_.bind();
    meshProg_.setUniformValue("uView", view);
    meshProg_.setUniformValue("uProj", proj);
    meshProg_.setUniformValue("uModel", QMatrix4x4());
    meshProg_.setUniformValue("uUseScalar", 0);
    g.vao->bind();
    for (size_t k = 0; k < ranges.size(); ++k) {
        const bool cloth = k < snap_->cloths.size();
        const Vector3& c = cloth ? snap_->cloths[k].color : snap_->softMeshes[k - snap_->cloths.size()].color;
        meshProg_.setUniformValue("uColor", QVector3D(c.x, c.y, c.z));
        meshProg_.setUniformValue("uUseScalar", cloth && !snap_->cloths[k].burnt.empty() ? 2 : 0); // 2: charring
        glDrawArrays(GL_TRIANGLES, ranges[k].first, ranges[k].second);
    }
    g.vao->release();
}

void Viewport::drawParticles(const QMatrix4x4& view, const QMatrix4x4& proj) {
    // Build one sprite batch: SPH particles (coloured by scalar) + spherical rigid bodies (own colour).
    const auto& P = snap_->particles;
    size_t nSph = 0;
    for (const auto& b : snap_->bodies) nSph += b.shape == rf::ShapeType::Sphere;
    const size_t total = P.size() + nSph;
    if (total == 0) return;
    if (particleFrame_ != snap_->frame || particleCount_ != int(total)) {
        std::vector<float> v;
        v.reserve(total * 8);
        const bool own = snap_->particleColor.size() == P.size();
        for (size_t i = 0; i < P.size(); ++i) {
            float s = i < snap_->particleScalar.size() ? snap_->particleScalar[i] : 0.0f;
            Vector3 c = own && s == rf::RenderSnapshot::kOwnColor ? snap_->particleColor[i] : Vector3(0.3f, 0.55f, 0.95f);
            v.insert(v.end(), {P[i].x, P[i].y, P[i].z, s, c.x, c.y, c.z, snap_->particleRadius});
        }
        for (const auto& b : snap_->bodies)
            if (b.shape == rf::ShapeType::Sphere)
            {
                float dim = b.sleeping ? 0.55f : 1.0f;
                v.insert(v.end(), {b.pos.x, b.pos.y, b.pos.z, -1e9f, b.color.x * dim, b.color.y * dim, b.color.z * dim, b.radius});
            }
        particles_.vbo.bind();
        particles_.vbo.allocate(v.data(), int(v.size() * sizeof(float)));
        particleFrame_ = snap_->frame;
        particleCount_ = int(total);
    }
    const float pointScale = height() / (2.0f * std::tan(qDegreesToRadians(camera_.fovDeg() * 0.5f)));
    sphereProg_.bind();
    sphereProg_.setUniformValue("uView", view);
    sphereProg_.setUniformValue("uProj", proj);
    sphereProg_.setUniformValue("uPointScale", float(pointScale * devicePixelRatioF()));
    sphereProg_.setUniformValue("uMin", snap_->colorMin);
    sphereProg_.setUniformValue("uMax", snap_->colorMax);
    sphereProg_.setUniformValue("uCmap", colormap_);
    particles_.vao.bind();
    if (!P.empty()) {
        bool uniform = snap_->vis.particleColoring == rf::ParticleColoring::Uniform;
        sphereProg_.setUniformValue("uUseScalar", uniform ? 0 : 1);
        sphereProg_.setUniformValue("uRadiusScale", particleScale_);
        glDrawArrays(GL_POINTS, 0, GLsizei(P.size()));
    }
    if (nSph) {
        sphereProg_.setUniformValue("uUseScalar", 0);
        sphereProg_.setUniformValue("uRadiusScale", 1.0f);
        glDrawArrays(GL_POINTS, GLint(P.size()), GLsizei(nSph));
    }
    particles_.vao.release();
}

void Viewport::drawSlice(const QMatrix4x4& vp) {
    if (!snap_->hasSlice || snap_->sliceW <= 0) return;
    if (sliceFrame_ != snap_->frame + snap_->paramsVersion * 1000003ull + uint64_t(snap_->sliceW) * 7 +
                           uint64_t(snap_->vis.sliceField) * 131 + uint64_t(snap_->vis.sliceAxis) * 17 +
                           uint64_t(snap_->vis.slicePosition * 1000)) {
        std::vector<float> tex(size_t(snap_->sliceW) * snap_->sliceH * 2);
        for (size_t i = 0; i < snap_->slice.size(); ++i) {
            tex[2 * i] = snap_->slice[i];
            tex[2 * i + 1] = snap_->sliceSolid[i] ? 1.0f : 0.0f;
        }
        glBindTexture(GL_TEXTURE_2D, sliceTex_);
        glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RG32F, snap_->sliceW, snap_->sliceH, 0, GL_RG, GL_FLOAT, tex.data());
        Vector3 c = snap_->sliceCorner, U = snap_->sliceU, V = snap_->sliceV;
        Vector3 p0 = c, p1 = c + U, p2 = c + U + V, p3 = c + V;
        float q[] = {p0.x, p0.y, p0.z, 0, 0, p1.x, p1.y, p1.z, 1, 0, p2.x, p2.y, p2.z, 1, 1,
                     p0.x, p0.y, p0.z, 0, 0, p2.x, p2.y, p2.z, 1, 1, p3.x, p3.y, p3.z, 0, 1};
        slice_.vbo.bind();
        slice_.vbo.allocate(q, sizeof(q));
        sliceFrame_ = snap_->frame + snap_->paramsVersion * 1000003ull + uint64_t(snap_->sliceW) * 7 +
                      uint64_t(snap_->vis.sliceField) * 131 + uint64_t(snap_->vis.sliceAxis) * 17 +
                      uint64_t(snap_->vis.slicePosition * 1000);
    }
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDisable(GL_CULL_FACE);
    sliceProg_.bind();
    sliceProg_.setUniformValue("uVP", vp);
    sliceProg_.setUniformValue("uMin", snap_->colorMin);
    sliceProg_.setUniformValue("uMax", snap_->colorMax);
    sliceProg_.setUniformValue("uAlpha", sliceOpacity_);
    sliceProg_.setUniformValue("uCmap", colormap_);
    sliceProg_.setUniformValue("uTex", 0);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, sliceTex_);
    slice_.vao.bind();
    glDrawArrays(GL_TRIANGLES, 0, 6);
    slice_.vao.release();
    glDisable(GL_BLEND);
}

void Viewport::drawStreamlines(const QMatrix4x4& vp) {
    const auto& L = snap_->streamlines;
    if (L.empty()) return;
    std::vector<float> data;
    std::vector<GLint> first;
    std::vector<GLsizei> count;
    for (size_t i = 0; i < L.size(); ++i) {
        if (L[i].size() < 2) continue;
        first.push_back(GLint(data.size() / 4));
        count.push_back(GLsizei(L[i].size()));
        for (size_t k = 0; k < L[i].size(); ++k)
            data.insert(data.end(), {L[i][k].x, L[i][k].y, L[i][k].z, snap_->streamlineSpeed[i][k]});
    }
    if (first.empty()) return;
    lineProg_.bind();
    lineProg_.setUniformValue("uVP", vp);
    lineProg_.setUniformValue("uColor", QVector4D(1, 1, 1, 1));
    lineProg_.setUniformValue("uUseScalar", 1);
    lineProg_.setUniformValue("uMin", 0.0f);
    lineProg_.setUniformValue("uMax", snap_->streamlineMax);
    lineProg_.setUniformValue("uCmap", colormap_);
    lines_.vao.bind();
    lines_.vbo.bind();
    lines_.vbo.allocate(data.data(), int(data.size() * sizeof(float)));
    for (size_t k = 0; k < first.size(); ++k) glDrawArrays(GL_LINE_STRIP, first[k], count[k]);
    lines_.vao.release();
}

void Viewport::drawFieldLines(const QMatrix4x4& vp) {
    // Magnetic field lines, brighter where the field is strong (log scale over two decades).
    const auto& L = snap_->fieldLines;
    if (L.empty()) return;
    std::vector<float> data;
    std::vector<GLint> first;
    std::vector<GLsizei> count;
    const float logMax = std::log10(std::max(snap_->fieldLineMax, 1e-12f));
    for (size_t i = 0; i < L.size(); ++i) {
        first.push_back(GLint(data.size() / 4));
        count.push_back(GLsizei(L[i].size()));
        for (size_t k = 0; k < L[i].size(); ++k) {
            const float s = std::clamp((std::log10(std::max(snap_->fieldLineStrength[i][k], 1e-12f)) - logMax) / 2.0f + 1.0f, 0.0f, 1.0f);
            data.insert(data.end(), {L[i][k].x, L[i][k].y, L[i][k].z, s});
        }
    }
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    lineProg_.bind();
    lineProg_.setUniformValue("uVP", vp);
    lineProg_.setUniformValue("uColor", QVector4D(1, 1, 1, 0.85f));
    lineProg_.setUniformValue("uUseScalar", 1);
    lineProg_.setUniformValue("uMin", 0.0f);
    lineProg_.setUniformValue("uMax", 1.0f);
    lineProg_.setUniformValue("uCmap", int(Viridis));
    lines_.vao.bind();
    lines_.vbo.bind();
    lines_.vbo.allocate(data.data(), int(data.size() * sizeof(float)));
    for (size_t k = 0; k < first.size(); ++k) glDrawArrays(GL_LINE_STRIP, first[k], count[k]);
    lines_.vao.release();
    glDisable(GL_BLEND);
}

void Viewport::drawVoxelGrid(const QMatrix4x4& vp) {
    const int mode = snap_->vis.gridDisplay;
    if (mode == 0 || snap_->gridNx <= 0) return;
    const int n[3] = {snap_->gridNx, snap_->gridNy, snap_->gridNz};
    const float dx = snap_->gridDx;
    const Vector3 o = snap_->gridOrigin;
    std::vector<float> g;
    auto line = [&](Vector3 a, Vector3 b) { g.insert(g.end(), {a.x, a.y, a.z, 0, b.x, b.y, b.z, 0}); };
    if (mode == 1) {
        // Cell boundaries of the slice layer, drawn in the plane through the cell centres.
        int axis = std::clamp(snap_->vis.sliceAxis, 0, 2);
        int ua = axis == 0 ? 2 : 0, va = axis == 1 ? 2 : 1;
        float planeCoord = o[axis] + (snap_->sliceLayer + 0.5f) * dx;
        for (int a = 0; a <= n[ua]; ++a) {
            Vector3 p0 = o, p1 = o;
            p0[axis] = p1[axis] = planeCoord;
            p0[ua] = p1[ua] = o[ua] + a * dx;
            p1[va] = o[va] + n[va] * dx;
            line(p0, p1);
        }
        for (int b = 0; b <= n[va]; ++b) {
            Vector3 p0 = o, p1 = o;
            p0[axis] = p1[axis] = planeCoord;
            p0[va] = p1[va] = o[va] + b * dx;
            p1[ua] = o[ua] + n[ua] * dx;
            line(p0, p1);
        }
    } else {
        // Whole lattice: lines along each axis through every grid node.
        for (int axis = 0; axis < 3; ++axis) {
            int ua = (axis + 1) % 3, va = (axis + 2) % 3;
            for (int a = 0; a <= n[ua]; ++a)
                for (int b = 0; b <= n[va]; ++b) {
                    Vector3 p0 = o, p1 = o;
                    p0[ua] = p1[ua] = o[ua] + a * dx;
                    p0[va] = p1[va] = o[va] + b * dx;
                    p1[axis] = o[axis] + n[axis] * dx;
                    line(p0, p1);
                }
        }
    }
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDepthMask(GL_FALSE);
    drawLines(g, GL_LINES, vp, QVector4D(0.75f, 0.8f, 0.9f, mode == 1 ? 0.35f : 0.07f), 1);
    glDepthMask(GL_TRUE);
    glDisable(GL_BLEND);
}

void Viewport::drawVectors(const QMatrix4x4& vp, const QVector3D& eyeQ) {
    const auto& P = snap_->arrowPos;
    const auto& V = snap_->arrowVel;
    if (P.empty()) return;
    const Vector3 eye(eyeQ.x(), eyeQ.y(), eyeQ.z());
    const float vmax = snap_->arrowMax, L = snap_->arrowLength;
    std::vector<float> d;
    d.reserve(P.size() * 6 * 4);
    for (size_t i = 0; i < P.size(); ++i) {
        float sp = rf::length(V[i]);
        if (sp < 1e-6f) continue;
        Vector3 dir = V[i] / sp;
        float len = L * sp / vmax;
        // Arrow centred on the cell centre.
        Vector3 a = P[i] - dir * (0.5f * len), b = P[i] + dir * (0.5f * len);
        Vector3 side = rf::normalize(rf::cross(dir, eye - b));
        Vector3 h1 = b - dir * (0.3f * len) + side * (0.15f * len);
        Vector3 h2 = b - dir * (0.3f * len) - side * (0.15f * len);
        d.insert(d.end(), {a.x, a.y, a.z, sp, b.x, b.y, b.z, sp, b.x, b.y, b.z, sp, h1.x, h1.y, h1.z, sp,
                           b.x, b.y, b.z, sp, h2.x, h2.y, h2.z, sp});
    }
    if (d.empty()) return;
    lineProg_.bind();
    lineProg_.setUniformValue("uVP", vp);
    lineProg_.setUniformValue("uColor", QVector4D(1, 1, 1, 1));
    lineProg_.setUniformValue("uUseScalar", 1);
    lineProg_.setUniformValue("uMin", 0.0f);
    lineProg_.setUniformValue("uMax", vmax);
    lineProg_.setUniformValue("uCmap", colormap_);
    lines_.vao.bind();
    lines_.vbo.bind();
    lines_.vbo.allocate(d.data(), int(d.size() * sizeof(float)));
    glDrawArrays(GL_LINES, 0, GLsizei(d.size() / 4));
    lines_.vao.release();
}

void Viewport::drawHeatSource(const QMatrix4x4& vp) {
    if (snap_->mode != rf::SimMode::WindTunnel || !snap_->heat.enabled) return;
    const Vector3 c = snap_->heat.center;
    const float r = snap_->heat.radius;
    std::vector<float> d;
    const int seg = 48;
    for (int plane = 0; plane < 3; ++plane)
        for (int i = 0; i < seg; ++i)
            for (int k = 0; k < 2; ++k) {
                float t = 2 * rf::kPi * (i + k) / seg;
                Vector3 p = c;
                p[(plane + 1) % 3] += r * std::cos(t);
                p[(plane + 2) % 3] += r * std::sin(t);
                d.insert(d.end(), {p.x, p.y, p.z, 0});
            }
    drawLines(d, GL_LINES, vp, QVector4D(1.0f, 0.6f, 0.15f, 1.0f), 1);
}

// Blits the depth of the frame drawn so far (the widget's framebuffer, possibly multisampled) into
// depthTex_. The formats must match for a depth blit: QOpenGLWidget's buffer is depth 24 + stencil 8.
bool Viewport::copySceneDepth(int w, int h) {
    if (!depthFbo_ || depthW_ != w || depthH_ != h) {
        if (!depthTex_) glGenTextures(1, &depthTex_);
        glBindTexture(GL_TEXTURE_2D, depthTex_);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH24_STENCIL8, w, h, 0, GL_DEPTH_STENCIL, GL_UNSIGNED_INT_24_8, nullptr);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        if (!depthFbo_) glGenFramebuffers(1, &depthFbo_);
        glBindFramebuffer(GL_FRAMEBUFFER, depthFbo_);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_TEXTURE_2D, depthTex_, 0);
        depthW_ = w;
        depthH_ = h;
    }
    glBindFramebuffer(GL_READ_FRAMEBUFFER, defaultFramebufferObject());
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, depthFbo_);
    const bool complete = glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE;
    while (glGetError() != GL_NO_ERROR) {}
    if (complete) glBlitFramebuffer(0, 0, w, h, 0, 0, w, h, GL_DEPTH_BUFFER_BIT, GL_NEAREST);
    const bool ok = complete && glGetError() == GL_NO_ERROR;
    glBindFramebuffer(GL_FRAMEBUFFER, defaultFramebufferObject());
    return ok;
}

void Viewport::drawVolume(const QMatrix4x4& vp, const QVector3D& eye) {
    if (!snap_->hasVolume || snap_->volume.empty()) return;
    const int W = int(width() * devicePixelRatioF()), H = int(height() * devicePixelRatioF());
    const bool depthOk = copySceneDepth(W, H);
    if (volumeFrame_ != snap_->frame + snap_->paramsVersion * 1000003ull) {
        glBindTexture(GL_TEXTURE_3D, volumeTex_);
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        const bool fire = snap_->volumeChannels == 2; // smoke + temperature
        glTexImage3D(GL_TEXTURE_3D, 0, fire ? GL_RG8 : GL_R8, snap_->volX, snap_->volY, snap_->volZ, 0, fire ? GL_RG : GL_RED,
                     GL_UNSIGNED_BYTE, snap_->volume.data());
        glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
        volumeFrame_ = snap_->frame + snap_->paramsVersion * 1000003ull;
    }
    const rf::AABB& d = snap_->domain;
    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_FRONT); // draw back faces so it also works with the camera inside the box
    volumeProg_.bind();
    volumeProg_.setUniformValue("uVP", vp);
    volumeProg_.setUniformValue("uLo", QVector3D(d.lo.x, d.lo.y, d.lo.z));
    Vector3 e = d.extent();
    volumeProg_.setUniformValue("uSize", QVector3D(e.x, e.y, e.z));
    volumeProg_.setUniformValue("uEye", eye);
    volumeProg_.setUniformValue("uDensity", smokeDensity_);
    volumeProg_.setUniformValue("uVol", 0);
    volumeProg_.setUniformValue("uSteps", softwareRenderer() ? 96 : 192); // the CPU renderer: half the samples
    const int mode = snap_->volumePlasma ? 2 : snap_->volumeChannels == 2 ? 1 : 0; // plasma / fire / smoke
    volumeProg_.setUniformValue("uFire", mode);
    volumeProg_.setUniformValue("uTempScale", snap_->volumeTemperatureScale);
    volumeProg_.setUniformValue("uAmbient", snap_->ambientTemperature);
    volumeProg_.setUniformValue("uFlame", mode == 2 ? 3.0f : 12.0f);
    volumeProg_.setUniformValue("uUseDepth", depthOk ? 1 : 0);
    volumeProg_.setUniformValue("uDepth", 1);
    volumeProg_.setUniformValue("uInvVP", vp.inverted());
    volumeProg_.setUniformValue("uViewport", QVector2D(float(W), float(H)));
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, depthTex_);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_3D, volumeTex_);
    volumeBox_.vao.bind();
    glDrawArrays(GL_TRIANGLES, 0, volumeBox_.count);
    volumeBox_.vao.release();
    glDisable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);
}

void Viewport::paintGL() {
    if (!glOk_) {
        QPainter p(this);
        p.fillRect(rect(), QColor(26, 28, 33));
        p.setPen(QColor(230, 200, 120));
        p.drawText(rect(), Qt::AlignCenter,
                   QString("OpenGL 3.0 недоступен (%1).\nЗапустите с ключом --software-gl — рендер на процессоре.").arg(renderer_));
        return;
    }
    glViewport(0, 0, int(width() * devicePixelRatioF()), int(height() * devicePixelRatioF()));
    // QPainter (overlay) may have changed state during the previous frame.
    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);
    glDepthFunc(GL_LESS);
    glEnable(GL_PROGRAM_POINT_SIZE);
    glDisable(GL_BLEND);
    glDisable(GL_CULL_FACE);
    glDisable(GL_SCISSOR_TEST);
    glDisable(GL_STENCIL_TEST);
    glClearColor(0.1f, 0.11f, 0.13f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    drawBackground();
    if (!snap_) return;

    const QMatrix4x4 view = viewMatrix(), proj = projMatrix(), vp = proj * view;
    const rf::AABB& d = snap_->domain;

    // Floor grid and domain box.
    if (showFloor_ && d.valid()) {
        std::vector<float> g;
        float y = d.lo.y;
        Vector3 c = d.center(), e = d.extent();
        float half = std::max(e.x, e.z) * 1.0f;
        float step = std::pow(10.0f, std::floor(std::log10(half))) * 0.5f;
        int n = int(half / step);
        for (int i = -n; i <= n; ++i) {
            float o = i * step;
            g.insert(g.end(), {c.x + o, y, c.z - n * step, 0, c.x + o, y, c.z + n * step, 0});
            g.insert(g.end(), {c.x - n * step, y, c.z + o, 0, c.x + n * step, y, c.z + o, 0});
        }
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        drawLines(g, GL_LINES, vp, QVector4D(0.5f, 0.55f, 0.62f, 0.18f), 1);
        glDisable(GL_BLEND);
    }
    if (showDomain_ && d.valid()) {
        std::vector<float> b;
        Vector3 lo = d.lo, hi = d.hi;
        Vector3 c[8];
        for (int i = 0; i < 8; ++i) c[i] = Vector3((i & 1) ? hi.x : lo.x, (i & 2) ? hi.y : lo.y, (i & 4) ? hi.z : lo.z);
        const int edges[12][2] = {{0, 1}, {2, 3}, {4, 5}, {6, 7}, {0, 2}, {1, 3}, {4, 6}, {5, 7}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
        for (auto& ed : edges)
            b.insert(b.end(), {c[ed[0]].x, c[ed[0]].y, c[ed[0]].z, 0, c[ed[1]].x, c[ed[1]].y, c[ed[1]].z, 0});
        drawLines(b, GL_LINES, vp, QVector4D(0.45f, 0.65f, 0.95f, 1.0f), 1);
    }

    drawObstacle(view, proj);
    drawBodies(view, proj);
    drawParticles(view, proj);
    drawCloths(view, proj);
    // The water surface goes after everything opaque: it refracts what lies behind it.
    if (!snap_->liquid.empty()) {
        const qreal dpr = devicePixelRatioF();
        fluidSurface_.render(snap_->liquid, snap_->liquidRadius, view, proj, int(width() * dpr), int(height() * dpr),
                             defaultFramebufferObject(), softwareRenderer());
        glUseProgram(0);
    }
    drawStreamlines(vp);
    drawFieldLines(vp);
    drawHeatSource(vp);
    drawSlice(vp);
    drawVolume(vp, eyePosition());
    // Diagnostic overlays go on top of the smoke so they stay readable.
    drawVoxelGrid(vp);
    drawVectors(vp, eyePosition());
    drawProbe(vp);
    drawJoints(vp);
    drawGrab(vp);

    // Emitter / heat source markers.
    if (snap_->mode == rf::SimMode::Fluid && snap_->emitter.enabled) {
        const auto& em = snap_->emitter;
        Vector3 a = em.position, b2 = em.position + rf::normalize(em.direction) * 0.15f;
        drawLines({a.x, a.y, a.z, 0, b2.x, b2.y, b2.z, 0}, GL_LINES, vp, QVector4D(1, 0.8f, 0.2f, 1), 1);
    }

    glBindVertexArray(0);
    glUseProgram(0);
    drawOverlay();

    ++fpsFrames_;
    qint64 now = QDateTime::currentMSecsSinceEpoch();
    if (now - fpsLast_ > 1000) {
        fps_ = fpsFrames_ * 1000.0f / float(now - fpsLast_);
        fpsFrames_ = 0;
        fpsLast_ = now;
    }
}

// ---------------------------------------------------------------------------
// 2D overlay
// ---------------------------------------------------------------------------
void Viewport::drawLegend(QPainter& p, const QRect& r, float lo, float hi, const QString& label, int map) {
    QLinearGradient g(r.bottomLeft(), r.topLeft());
    for (int i = 0; i <= 16; ++i) g.setColorAt(i / 16.0, colormapColor(map, i / 16.0f));
    p.setPen(QColor(20, 20, 24));
    p.setBrush(g);
    p.drawRoundedRect(r, 3, 3);
    p.setPen(QColor(225, 228, 235));
    QFont f = p.font();
    f.setPointSizeF(8.5);
    p.setFont(f);
    for (int i = 0; i <= 4; ++i) {
        float t = i / 4.0f;
        int y = r.bottom() - int(t * r.height());
        float v = lo + (hi - lo) * t;
        QString txt = (std::fabs(v) >= 1000 || (std::fabs(v) < 0.01f && v != 0)) ? QString::number(v, 'e', 1)
                                                                                  : QString::number(v, 'f', 2);
        p.drawLine(r.right() + 1, y, r.right() + 5, y);
        p.drawText(r.right() + 8, y + 4, txt);
    }
    f.setBold(true);
    p.setFont(f);
    p.drawText(r.left() - 4, r.top() - 8, label);
}

void Viewport::drawOverlay() {
    QPainter p(this);
    p.setRenderHint(QPainter::Antialiasing);
    const int W = width(), H = height();

    // Title
    QFont f = p.font();
    f.setPointSizeF(10.5);
    f.setBold(true);
    p.setFont(f);
    p.setPen(QColor(235, 238, 245));
    p.drawText(14, 24, QString::fromStdString(rf::presetName(snap_->preset)));
    f.setBold(false);
    f.setPointSizeF(9);
    p.setFont(f);
    p.setPen(QColor(170, 176, 188));
    const bool cpu = softwareRenderer();
    p.drawText(14, 42, QString("t = %1 с   ·   расчёт %2 мс/кадр   ·   рендер %3 fps (%4)")
                           .arg(snap_->time, 0, 'f', 3)
                           .arg(snap_->stepMs, 0, 'f', 1)
                           .arg(fps_, 0, 'f', 0)
                           .arg(cpu ? QString("процессор, llvmpipe") : QString("GPU")));

    // Legends
    bool mainLegend = !snap_->colorLabel.empty() &&
                      (snap_->hasSlice || (!snap_->particles.empty() &&
                                           snap_->vis.particleColoring != rf::ParticleColoring::Uniform));
    int x = W - 110;
    if (mainLegend) {
        drawLegend(p, QRect(x, H - 230, 16, 180), snap_->colorMin, snap_->colorMax,
                   QString::fromStdString(snap_->colorLabel), colormap_);
        x -= 110;
    }
    if (!snap_->arrowPos.empty() && !snap_->hasSlice) {
        drawLegend(p, QRect(x, H - 230, 16, 180), 0.0f, snap_->arrowMax, "|V| векторы, м/с", colormap_);
        x -= 110;
    }
    if (!snap_->obstacleScalar.empty()) {
        float mn = 1e9f, mx = -1e9f;
        for (float v : snap_->obstacleScalar) { mn = std::min(mn, v); mx = std::max(mx, v); }
        float lo = std::max(mn, -3.0f), hi = std::min(std::max(mx, lo + 0.1f), 1.2f);
        drawLegend(p, QRect(x, H - 230, 16, 180), lo, hi, "Cp пов.", CoolWarm);
    }

    // Axis gizmo
    QMatrix4x4 v = viewMatrix();
    QPointF o(52, H - 52);
    const QVector3D axes[3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    const QColor cols[3] = {QColor(235, 90, 90), QColor(110, 220, 110), QColor(100, 150, 255)};
    const char* names[3] = {"X", "Y", "Z"};
    for (int i = 0; i < 3; ++i) {
        QVector3D a = v.mapVector(axes[i]);
        QPointF e = o + QPointF(a.x(), -a.y()) * 32;
        p.setPen(QPen(cols[i], 2.2));
        p.drawLine(o, e);
        p.drawText(e + QPointF(3, 4), names[i]);
    }
    if (snap_->mode == rf::SimMode::WindTunnel) {
        p.setPen(QColor(170, 176, 188));
        p.drawText(14, H - 12, "ЛКМ по телу — схватить, по газу — возмущение  ·  колесо зажато — вращение, Shift+колесо / ПКМ — сдвиг, прокрутка — масштаб  ·  F — показать всё");
    } else {
        p.setPen(QColor(170, 176, 188));
        p.drawText(90, H - 12, "ЛКМ по телу — схватить и тащить, ЛКМ по пустоте / колесо зажато — вращение, Shift+колесо / ПКМ — сдвиг, прокрутка — масштаб");
    }
    p.end();
}
