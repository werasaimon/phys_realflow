#pragma once
// The GLSL sources of the viewport's programs, GLSL 1.30 (OpenGL 3.0, also Mesa's llvmpipe) with
// the version header prepended by Viewport::initializeGL (buildProgram), 3.30 core where available.

// Colormaps (turbo, cool-warm, gray, viridis) shared by the fragment shaders that colour by a field.
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

// Background: a full-screen gradient.
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

// Lit meshes: the obstacle, rigid bodies, cloth; also the glass vessel (uAlpha < 1).
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
uniform float uAlpha = 1.0; // < 1: glass (drawn blended, after the volume)
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
    // Glass: the rim (grazing view) is brighter, the face almost clear (Schlick's Fresnel).
    float alpha = uAlpha < 1.0 ? uAlpha + (1.0 - uAlpha) * pow(1.0 - max(dot(n, v), 0.0), 4.0) : 1.0;
    o = vec4(c, alpha);
})";

// Particles as sphere impostors (point sprites shaded as spheres).
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

// The magnetised sphere of the magnetosphere scene drawn as a planet (procedural Earth, aurora).
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

// Lines: streamlines, field lines, joints, the domain box.
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

// A coloured slice of a grid field.
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

// The smoke / fire / plasma volume: ray marching through the 3D texture against the scene depth.
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
