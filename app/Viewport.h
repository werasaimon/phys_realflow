#pragma once
// OpenGL viewport: orbit camera, particles as ray-traced sphere impostors, meshes, field slice,
// streamlines, ray-marched smoke volume and a colour legend overlay.
// Needs only OpenGL 3.0 (GLSL 1.30): runs on any GPU of the last ~15 years and on the CPU through
// Mesa llvmpipe (--software-gl). A 3.3 core context is used when the driver offers one.

#include "FluidSurfaceRenderer.h"
#include "OrbitCamera.h"
#include "ViewportTools.h"
#include "scene/Simulation.h"

#include <QColor>
#include <QMatrix4x4>
#include <QOpenGLBuffer>
#include <QOpenGLExtraFunctions>
#include <QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLWidget>
#include <QVector3D>

#include <map>
#include <memory>
#include <string>

class Viewport : public QOpenGLWidget, protected QOpenGLExtraFunctions {
    Q_OBJECT
public:
    enum Colormap { Turbo, Viridis, CoolWarm, Gray };

    explicit Viewport(QWidget* parent = nullptr);
    ~Viewport() override;

    void setSnapshot(std::shared_ptr<const rf::RenderSnapshot> s);
    void setColormap(int c) { colormap_ = c; update(); }
    void setShowDomain(bool on) { showDomain_ = on; update(); }
    void setShowFloor(bool on) { showFloor_ = on; update(); }
    // Draw non-convex bodies as their convex decomposition instead of the smooth model.
    void setShowConvexParts(bool on) { showConvexParts_ = on; update(); }
    void setParticleScale(float s) { particleScale_ = s; update(); }
    void setSmokeDensity(float d) { smokeDensity_ = d; update(); }
    void setSliceOpacity(float a) { sliceOpacity_ = a; update(); }
    void frameScene();

    static QColor colormapColor(int map, float t);
    float renderFps() const { return fps_; }

    // Interaction (used by the ViewportTool strategies)
    OrbitCamera& camera() { return camera_; }
    // Point on the cursor ray inside the domain: on the slice plane if hit, else mid-domain.
    bool pickPoint(const Ray& ray, rf::Vector3& hit) const;
    // Nearest rigid body hit by the ray (exact ray cast against its convex shape).
    bool pickBody(const Ray& ray, int& body, rf::Vector3& hit) const;
    // Nearest particle (liquid, soft body) or cloth vertex hit by the ray; hit = its centre.
    bool pickParticle(const Ray& ray, rf::Vector3& hit) const;
    void showProbe(const Ray& ray, const rf::Vector3& hit);
    void setBrushRadius(float r) { brushRadius_ = r; }

signals:
    // World-space point and velocity of a mouse stroke (velocity is zero on the first click).
    void disturbanceRequested(rf::Vector3 position, rf::Vector3 velocity);
    void grabStarted(int body, rf::Vector3 point);
    void particleGrabStarted(rf::Vector3 point); // cloth / soft body / liquid particle under the cursor
    void grabMoved(rf::Vector3 target);
    void grabReleased();
    // The GPU / driver cannot do OpenGL 3.0: the window offers a restart in software mode.
    void openGLUnsupported(QString renderer);

protected:
    void initializeGL() override;
    void resizeGL(int w, int h) override;
    void paintGL() override;
    void mousePressEvent(QMouseEvent* e) override;
    void mouseMoveEvent(QMouseEvent* e) override;
    void mouseReleaseEvent(QMouseEvent* e) override;
    void wheelEvent(QWheelEvent* e) override;
    void mouseDoubleClickEvent(QMouseEvent* e) override;
    void keyPressEvent(QKeyEvent* e) override;

private:
    struct Buffer {
        QOpenGLVertexArrayObject vao;
        QOpenGLBuffer vbo{QOpenGLBuffer::VertexBuffer};
        int count = 0;
    };

    void buildShaders();
    void uploadObstacle(const rf::TriMesh& mesh);
    void drawBackground();
    void drawLines(const std::vector<float>& data, GLenum mode, const QMatrix4x4& vp, const QVector4D& color, float width);
    void drawObstacle(const QMatrix4x4& view, const QMatrix4x4& proj);
    void drawVesselGlass(const QMatrix4x4& view, const QMatrix4x4& proj);
    void drawBodies(const QMatrix4x4& view, const QMatrix4x4& proj);
    void drawParticles(const QMatrix4x4& view, const QMatrix4x4& proj);
    void drawCloths(const QMatrix4x4& view, const QMatrix4x4& proj); // cloth sheets + soft body surfaces
    void drawSlice(const QMatrix4x4& vp);
    void drawStreamlines(const QMatrix4x4& vp);
    void drawFieldLines(const QMatrix4x4& vp); // magnetic field lines
    void drawVoxelGrid(const QMatrix4x4& vp);
    void drawVectors(const QMatrix4x4& vp, const QVector3D& eye);
    void drawHeatSource(const QMatrix4x4& vp);
    void drawProbe(const QMatrix4x4& vp);
    void drawGrab(const QMatrix4x4& vp);
    void drawJoints(const QMatrix4x4& vp);
    ViewportTool& toolFor(Qt::MouseButton button, const QPointF& pos);
    void drawVolume(const QMatrix4x4& vp, const QVector3D& eye);
    void drawOverlay();
    void drawLegend(class QPainter& p, const QRect& r, float lo, float hi, const QString& label, int map);
    QMatrix4x4 viewMatrix() const { return camera_.view(); }
    QMatrix4x4 projMatrix() const { return camera_.projection(float(width()) / std::max(1, height())); }
    QVector3D eyePosition() const { return camera_.eye(); }

    std::shared_ptr<const rf::RenderSnapshot> snap_;

    QOpenGLShaderProgram meshProg_, sphereProg_, lineProg_, sliceProg_, volumeProg_, bgProg_, planetProg_;
    Buffer bg_, obstacle_, cube_, particles_, lines_, slice_, volumeBox_;
    QOpenGLBuffer obstacleScalar_{QOpenGLBuffer::VertexBuffer};
    std::vector<uint32_t> obstacleVertexMap_; // expanded vertex -> original vertex
    struct GpuMesh {
        std::unique_ptr<QOpenGLVertexArrayObject> vao;
        std::unique_ptr<QOpenGLBuffer> vbo;
        int count = 0;
        // Hull cache: the mesh is held here too, so its address (the cache key) cannot be taken by a
        // new mesh while the entry exists; the entry is dropped once only the cache holds the mesh.
        std::shared_ptr<const rf::TriMesh> mesh;
    };
    std::map<const rf::TriMesh*, GpuMesh> hullCache_;
    void pruneHullCache(); // needs the GL context current (paintGL)
    GpuMesh clothMesh_; // rebuilt every frame
    FluidSurfaceRenderer fluidSurface_; // liquid as a water surface (screen space)
    GpuMesh& hullMesh(const std::shared_ptr<const rf::TriMesh>& m, bool smooth = false);
    uint64_t obstacleVersion_ = 0;
    GLuint sliceTex_ = 0, volumeTex_ = 0;
    // Copy of the scene's depth (after the opaque pass): the volume rays stop at solid surfaces.
    GLuint depthTex_ = 0, depthFbo_ = 0;
    int depthW_ = 0, depthH_ = 0;
    bool copySceneDepth(int w, int h);
    // Every snapshot gets a new serial; the sprite, slice and smoke buffers are rebuilt once per new
    // snapshot (the frame number is no key: a paused reset or preset switch stays at frame 0).
    uint64_t snapSerial_ = 0, particleSerial_ = ~0ull, sliceSerial_ = ~0ull, volumeSerial_ = ~0ull;

    QString renderer_;   // OpenGL renderer (GPU name, or llvmpipe on the CPU)
    bool softwareRenderer() const {
        return renderer_.contains("llvmpipe", Qt::CaseInsensitive) || renderer_.contains("softpipe", Qt::CaseInsensitive);
    }
    bool glOk_ = true;   // OpenGL 3.0 or newer available
    QString glslHeader_; // "#version 330 core" on 3.3 core contexts, "#version 130" on 3.0
    OrbitCamera camera_;
    OrbitTool orbitTool_;
    DisturbTool disturbTool_;
    GrabTool grabTool_;
    ViewportTool* grabbed_ = nullptr; // tool that received the press, keeps the drag
    Ray probeRay_;
    rf::Vector3 probeHit_;
    qint64 probeTime_ = 0;
    float brushRadius_ = 0.15f;
    bool framedOnce_ = false;
    std::string lastSceneName_; // the scene is framed again when another one is loaded

    int colormap_ = Turbo;
    bool showDomain_ = true, showFloor_ = true, showConvexParts_ = false;
    float particleScale_ = 1.0f, smokeDensity_ = 16.0f, sliceOpacity_ = 0.92f;

    // fps
    qint64 fpsLast_ = 0;
    int fpsFrames_ = 0;
    float fps_ = 0;
};
