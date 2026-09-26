#include "MainWindow.h"

#include <QApplication>
#include <QCommandLineParser>
#include <QPalette>
#include <QStyleFactory>
#include <QSurfaceFormat>

static void applyDarkTheme(QApplication& app) {
    app.setStyle(QStyleFactory::create("Fusion"));
    QPalette p;
    const QColor bg(37, 40, 46), base(28, 30, 35), text(220, 223, 228), accent(64, 150, 255);
    p.setColor(QPalette::Window, bg);
    p.setColor(QPalette::WindowText, text);
    p.setColor(QPalette::Base, base);
    p.setColor(QPalette::AlternateBase, bg);
    p.setColor(QPalette::ToolTipBase, base);
    p.setColor(QPalette::ToolTipText, text);
    p.setColor(QPalette::Text, text);
    p.setColor(QPalette::Button, QColor(48, 52, 60));
    p.setColor(QPalette::ButtonText, text);
    p.setColor(QPalette::Highlight, accent);
    p.setColor(QPalette::HighlightedText, Qt::white);
    p.setColor(QPalette::Link, accent);
    p.setColor(QPalette::PlaceholderText, QColor(130, 135, 145));
    p.setColor(QPalette::Disabled, QPalette::Text, QColor(120, 124, 130));
    p.setColor(QPalette::Disabled, QPalette::ButtonText, QColor(120, 124, 130));
    p.setColor(QPalette::Disabled, QPalette::WindowText, QColor(120, 124, 130));
    app.setPalette(p);
    app.setStyleSheet(
        "QGroupBox { border: 1px solid #3a3f48; border-radius: 6px; margin-top: 14px; padding-top: 6px; font-weight: 600; }"
        "QGroupBox::title { subcontrol-origin: margin; left: 10px; padding: 0 4px; color: #8fc1ff; }"
        "QDockWidget::title { background: #2c3038; padding: 5px; font-weight: 600; }"
        "QToolBar { spacing: 6px; padding: 4px; border: none; }"
        "QToolButton { padding: 4px 8px; border-radius: 4px; }"
        "QToolButton:hover { background: #3a404b; }"
        "QToolButton:checked { background: #2d5c99; }"
        "QTableWidget { gridline-color: #333842; }"
        "QHeaderView::section { background: #2c3038; padding: 4px; border: none; }");
}

// Software rendering for PCs without a usable GPU / OpenGL driver: Qt's Mesa llvmpipe
// (opengl32sw.dll, an OpenGL 3.0 context) runs the same shaders on the CPU. Chosen by --software-gl or the
// environment variable RF_SOFTWARE_GL=1; it must be decided before the application object exists.
static bool wantsSoftwareGL(int argc, char** argv) {
    for (int i = 1; i < argc; ++i)
        if (QString::fromLocal8Bit(argv[i]) == "--software-gl") return true;
    return qEnvironmentVariableIntValue("RF_SOFTWARE_GL") != 0;
}

int main(int argc, char** argv) {
    const bool software = wantsSoftwareGL(argc, argv);
    if (software) QCoreApplication::setAttribute(Qt::AA_UseSoftwareOpenGL);
    QSurfaceFormat fmt;
    fmt.setVersion(3, 3);
    fmt.setProfile(QSurfaceFormat::CoreProfile);
    fmt.setDepthBufferSize(24);
    fmt.setSamples(software ? 0 : 4); // multisampling is costly on the CPU
    fmt.setSwapInterval(1);
    QSurfaceFormat::setDefaultFormat(fmt);

    QApplication app(argc, argv);
    QApplication::setApplicationName("PhysRealFlow");
    QApplication::setOrganizationName("PhysRealFlow");
    applyDarkTheme(app);

    QCommandLineParser cli;
    cli.setApplicationDescription("PhysRealFlow");
    cli.addHelpOption();
    QCommandLineOption presetOpt("preset", "Scene preset index to load.", "index");
    QCommandLineOption framesOpt("frames", "Frames to simulate before the screenshot.", "n", "120");
    QCommandLineOption shotOpt("screenshot", "Save a PNG of the viewport and quit.", "file");
    QCommandLineOption csvOpt("csv", "Save the plotted time series as CSV and quit.", "file");
    QCommandLineOption sizeOpt("size", "Window size WxH.", "size", "1600x950");
    QCommandLineOption softOpt("software-gl", "Render on the CPU (Mesa llvmpipe) - for PCs without a usable GPU.");
    cli.addOptions({presetOpt, framesOpt, shotOpt, csvOpt, sizeOpt, softOpt});
    cli.process(app);

    MainWindow w;
    QStringList wh = cli.value(sizeOpt).split('x');
    w.resize(wh.value(0).toInt() > 0 ? wh.value(0).toInt() : 1600, wh.value(1).toInt() > 0 ? wh.value(1).toInt() : 950);
    w.show();
    if (cli.isSet(shotOpt) || cli.isSet(csvOpt))
        w.runAutomation(cli.value(presetOpt).toInt(), cli.value(framesOpt).toInt(), cli.value(shotOpt), cli.value(csvOpt));
    else if (cli.isSet(presetOpt)) {
        int p = cli.value(presetOpt).toInt();
        w.runAutomation(p, 0, {}, {});
    }
    return app.exec();
}
