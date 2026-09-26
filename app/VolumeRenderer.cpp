// The smoke, fire and plasma volume of the viewport: ray marching through the 3D texture of the
// grid's tracer (and temperature) against the depth of the scene drawn so far, so bodies and the
// obstacle occlude it correctly. Methods of Viewport.
#include "Viewport.h"

using rf::Vector3;

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
    if (volumeSerial_ != snapSerial_) {
        glBindTexture(GL_TEXTURE_3D, volumeTex_);
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        const bool fire = snap_->volumeChannels == 2; // smoke + temperature
        glTexImage3D(GL_TEXTURE_3D, 0, fire ? GL_RG8 : GL_R8, snap_->volX, snap_->volY, snap_->volZ, 0, fire ? GL_RG : GL_RED,
                     GL_UNSIGNED_BYTE, snap_->volume.data());
        glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
        volumeSerial_ = snapSerial_;
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
