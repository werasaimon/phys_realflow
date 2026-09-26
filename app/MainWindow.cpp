#include "MainWindow.h"

#include "ParamForm.h"
#include "PlotPanel.h"
#include "Viewport.h"

#include <QActionGroup>
#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QDateTime>
#include <QDockWidget>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QMenuBar>
#include <QMessageBox>
#include <QPointer>
#include <QProcess>
#include <QPushButton>
#include <QScrollArea>
#include <QSplitter>
#include <QStatusBar>
#include <QStyle>
#include <QTableWidget>
#include <QTimer>
#include <QToolBar>
#include <QVBoxLayout>

#include <random>

using namespace rf;
using Snap = RenderSnapshot;

static const QStringList kShapes = {"Нет", "Сфера", "Куб", "Цилиндр", "Крыло NACA", "Обтекаемое тело",
                                    "Эллипсоид", "Конус", "Импортированная модель"};
static const QStringList kFields = {"Скорость |V|", "Давление p", "Коэф. давления Cp", "Завихренность |ω|",
                                    "Дым (плотность)", "Температура", "Скорость Vx", "Скорость Vy",
                                    "Магнитное поле |B|", "Плотность тока |J|"};
static const QStringList kBC = {"Стенка", "Вход потока", "Выход (p = 0)"};

MainWindow::MainWindow() {
    setWindowTitle("PhysRealFlow");
    ctrl_ = std::make_unique<SimController>();

    view_ = new Viewport(this);
    view_->setBrushRadius(brush_.radius);
    setCentralWidget(view_);
    connect(view_, &Viewport::disturbanceRequested, this, &MainWindow::onDisturbance);
    // Mouse joint: pick with the cursor ray, drag the target, release.
    connect(view_, &Viewport::openGLUnsupported, this, [this](QString renderer) {
        auto answer = QMessageBox::question(
            this, "OpenGL",
            QString("Видеокарта или драйвер не поддерживает OpenGL 3.0 (%1).\n\n"
                    "Перезапустить программу в программном режиме? Всё будет рисоваться процессором "
                    "(медленнее, но работает на любом ПК).").arg(renderer));
        if (answer != QMessageBox::Yes) return;
        QStringList args = QCoreApplication::arguments().mid(1);
        args << "--software-gl";
        QProcess::startDetached(QCoreApplication::applicationFilePath(), args);
        qApp->quit();
    });
    connect(view_, &Viewport::grabStarted, this, [this](int body, Vector3 p) {
        ctrl_->post([body, p](Simulation& s) { s.rigid.grab(body, p); });
    });
    connect(view_, &Viewport::particleGrabStarted, this, [this](Vector3 p) {
        ctrl_->post([p](Simulation& s) { s.particles.grab(p); });
    });
    connect(view_, &Viewport::grabMoved, this, [this](Vector3 t) {
        ctrl_->post([t](Simulation& s) {
            s.rigid.setGrabTarget(t);
            s.particles.setGrabTarget(t);
        });
    });
    connect(view_, &Viewport::grabReleased, this, [this] {
        ctrl_->post([](Simulation& s) {
            s.rigid.releaseGrab();
            s.particles.releaseGrab();
        });
    });

    buildActions();
    buildParameterDock();
    buildVisualDock();
    buildResultsDock();

    status_ = new QLabel;
    statusBar()->addPermanentWidget(status_, 1);

    timer_ = new QTimer(this);
    connect(timer_, &QTimer::timeout, this, &MainWindow::poll);
    timer_->start(15);
    poll();
}

MainWindow::~MainWindow() {
    timer_->stop();
    ctrl_.reset(); // stop the simulation thread before the widgets go away
}

// ---------------------------------------------------------------------------
// Toolbar & menu
// ---------------------------------------------------------------------------
void MainWindow::buildActions() {
    auto* tb = addToolBar("Управление");
    tb->setObjectName("mainToolbar");
    tb->setMovable(false);
    tb->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    tb->setIconSize(QSize(18, 18));

    playAct_ = new QAction(style()->standardIcon(QStyle::SP_MediaPlay), "Старт", this);
    playAct_->setShortcut(Qt::Key_Space);
    playAct_->setCheckable(true);
    connect(playAct_, &QAction::triggered, this, &MainWindow::togglePlay);

    auto* stepAct = new QAction(style()->standardIcon(QStyle::SP_MediaSkipForward), "Шаг", this);
    stepAct->setShortcut(Qt::Key_S);
    connect(stepAct, &QAction::triggered, this, [this] { ctrl_->requestStep(); });

    auto* resetAct = new QAction(style()->standardIcon(QStyle::SP_BrowserReload), "Сброс", this);
    resetAct->setShortcut(Qt::Key_R);
    connect(resetAct, &QAction::triggered, this, [this] { ctrl_->post([](Simulation& s) { s.reset(); }); });

    auto* importAct = new QAction(style()->standardIcon(QStyle::SP_DialogOpenButton), "Импорт модели…", this);
    importAct->setShortcut(QKeySequence::Open);
    connect(importAct, &QAction::triggered, this, &MainWindow::importMesh);
    auto* shotAct = new QAction("Скриншот…", this);
    shotAct->setShortcut(QKeySequence("Ctrl+P"));
    connect(shotAct, &QAction::triggered, this, &MainWindow::screenshot);
    auto* csvAct = new QAction("Экспорт графиков CSV…", this);
    csvAct->setShortcut(QKeySequence("Ctrl+E"));
    connect(csvAct, &QAction::triggered, this, &MainWindow::exportCsv);
    auto* loadsAct = new QAction("Экспорт нагрузок по полигонам CSV…", this);
    loadsAct->setToolTip("Cp, Cf и сила на каждом треугольнике тела в аэротрубе");
    connect(loadsAct, &QAction::triggered, this, &MainWindow::exportSurfaceLoads);

    presetCombo_ = new QComboBox;
    presetCombo_->setMinimumWidth(320);
    SimMode prevMode = presetMode(Preset(0));
    for (int p = 0; p < int(Preset::Count); ++p) {
        if (presetMode(Preset(p)) != prevMode) {
            presetCombo_->insertSeparator(presetCombo_->count());
            prevMode = presetMode(Preset(p));
        }
        presetCombo_->addItem(QString::fromStdString(presetName(Preset(p))), p);
    }
    connect(presetCombo_, &QComboBox::activated, this, [this](int idx) {
        int p = presetCombo_->itemData(idx).toInt();
        ctrl_->post([p](Simulation& s) { s.loadPreset(Preset(p)); });
    });

    tb->addWidget(new QLabel("  Сцена "));
    tb->addWidget(presetCombo_);
    tb->addSeparator();
    tb->addAction(playAct_);
    tb->addAction(stepAct);
    tb->addAction(resetAct);

    auto* file = menuBar()->addMenu("&Файл");
    file->addAction(importAct);
    file->addAction(csvAct);
    file->addAction(loadsAct);
    file->addAction(shotAct);
    file->addSeparator();
    file->addAction("Выход", QKeySequence::Quit, this, &QWidget::close);

    auto* sim = menuBar()->addMenu("&Симуляция");
    sim->addAction(playAct_);
    sim->addAction(stepAct);
    sim->addAction(resetAct);
    auto* rt = sim->addAction("Не быстрее реального времени");
    rt->setCheckable(true);
    rt->setChecked(true);
    connect(rt, &QAction::toggled, this, [this](bool on) { ctrl_->setRealtimeLimit(on); });

    auto* viewMenu = menuBar()->addMenu("&Вид");
    viewMenu->setObjectName("viewMenu");
    viewMenu->addAction("Показать всю сцену", QKeySequence(Qt::Key_F), view_, &Viewport::frameScene);
    viewMenu->addSeparator();

    auto* help = menuBar()->addMenu("&Справка");
    help->addAction("Управление и методы", QKeySequence::HelpContents, this, &MainWindow::showHelp);
    help->addAction("О программе", this, &MainWindow::showAbout);
}

void MainWindow::togglePlay() {
    bool on = !ctrl_->isRunning();
    ctrl_->setRunning(on);
    playAct_->setChecked(on);
    playAct_->setText(on ? "Пауза" : "Старт");
    playAct_->setIcon(style()->standardIcon(on ? QStyle::SP_MediaPause : QStyle::SP_MediaPlay));
}

// ---------------------------------------------------------------------------
// Parameter forms. Each form shows the few parameters people change most; everything else sits
// in its collapsed "Дополнительно" section. Rows marked * restart the scene.
// ---------------------------------------------------------------------------
QGroupBox* MainWindow::addGroup(QVBoxLayout* col, const QString& title, ParamForm* form) {
    auto* g = new QGroupBox(title);
    auto* l = new QVBoxLayout(g);
    l->setContentsMargins(4, 8, 4, 4);
    l->addWidget(form);
    col->addWidget(g);
    forms_.push_back(form);
    return g;
}

ParamForm* MainWindow::buildObjectForm() {
    auto* f = new ParamForm(ctrl_.get());
    f->addCombo("Форма", kShapes, [](const Snap& s) { return int(s.obstacleSettings.shape); },
                [](Simulation& s, int v) {
                    if (ObstacleShape(v) == ObstacleShape::Custom && !s.obstacle.customMesh) return;
                    s.obstacle.shape = ObstacleShape(v);
                    s.rebuildObstacle();
                });
    f->addDouble("Размер", 0.02, 5, 0.05, 3, [](const Snap& s) { return s.obstacleSettings.size; },
                 [](Simulation& s, double v) { s.obstacle.size = float(v); s.rebuildObstacle(); },
                 "Характерный размер: диаметр / хорда / ребро / длина", "м");
    f->addDouble("Угол атаки", -90, 90, 1, 1, [](const Snap& s) { return s.obstacleSettings.angleOfAttackDeg; },
                 [](Simulation& s, double v) { s.obstacle.angleOfAttackDeg = float(v); s.rebuildObstacle(); },
                 "Поворот вокруг оси Z, нос вверх при положительном угле", "°");
    auto* imp = new QPushButton("Импорт OBJ / STL…");
    connect(imp, &QPushButton::clicked, this, &MainWindow::importMesh);
    f->addRow(imp);

    f->beginAdvanced();
    f->addDouble("Размах / длина", 0.05, 5, 0.05, 3, [](const Snap& s) { return s.obstacleSettings.span; },
                 [](Simulation& s, double v) { s.obstacle.span = float(v); s.rebuildObstacle(); }, "Размах крыла, длина цилиндра", "м");
    f->addDouble("Относит. толщина", 0.05, 1, 0.02, 2, [](const Snap& s) { return s.obstacleSettings.thickness; },
                 [](Simulation& s, double v) { s.obstacle.thickness = float(v); s.rebuildObstacle(); });
    f->addText("Профиль NACA", [](const Snap& s) { return QString::fromStdString(s.obstacleSettings.nacaCode); },
               [](Simulation& s, QString v) {
                   if (v.size() == 4) { s.obstacle.nacaCode = v.toStdString(); s.rebuildObstacle(); }
               },
               "Четырёхзначный код: 0012, 2412, 4415…");
    f->addDouble("Рыскание", -180, 180, 5, 1, [](const Snap& s) { return s.obstacleSettings.yawDeg; },
                 [](Simulation& s, double v) { s.obstacle.yawDeg = float(v); s.rebuildObstacle(); }, {}, "°");
    f->addDouble("Крен", -180, 180, 5, 1, [](const Snap& s) { return s.obstacleSettings.rollDeg; },
                 [](Simulation& s, double v) { s.obstacle.rollDeg = float(v); s.rebuildObstacle(); }, {}, "°");
    for (int a = 0; a < 3; ++a) {
        const char* names[3] = {"Позиция X", "Позиция Y", "Позиция Z"};
        f->addDouble(names[a], -5, 5, 0.05, 3, [a](const Snap& s) { return s.obstacleSettings.position[a]; },
                     [a](Simulation& s, double v) { s.obstacle.position[a] = float(v); s.rebuildObstacle(); }, {}, "м");
    }
    f->addHint("Изменение геометрии перезапускает расчёт.");
    return f;
}

ParamForm* MainWindow::buildFluidForm() {
    auto* f = new ParamForm(ctrl_.get());
    f->addDouble("Радиус частицы *", 3, 60, 1, 1, [](const Snap& s) { return s.particleParams.particleRadius * 1000; },
                 [](Simulation& s, double v) { s.particles.params.particleRadius = float(v / 1000); s.reset(); },
                 "Меньше радиус — больше частиц, точнее и медленнее", "мм");
    f->addDouble("Вязкость", 0, 1, 0.01, 3, [](const Snap& s) { return s.particleParams.viscosity; },
                 [](Simulation& s, double v) { s.particles.params.viscosity = float(v); }, "XSPH: 0 — вода, 0.3+ — мёд");
    f->addBool("Струя из источника", [](const Snap& s) { return s.emitter.enabled; },
               [](Simulation& s, bool v) { s.particles.emitter.enabled = v; });
    auto* addBlock = new QPushButton("Добавить объём жидкости");
    connect(addBlock, &QPushButton::clicked, this, [this] {
        ctrl_->post([](Simulation& s) {
            AABB d = s.particles.domain();
            Vector3 c = d.center(), e = d.extent();
            s.particles.addBlock(AABB({c.x - 0.15f * e.x, d.hi.y - 0.35f * e.y, c.z - 0.25f * e.z},
                                {c.x + 0.15f * e.x, d.hi.y - 0.05f * e.y, c.z + 0.25f * e.z}));
        });
    });
    f->addRow(addBlock);

    f->beginAdvanced();
    // The particle mass is calibrated from ρ0 when the scene is built: a new density needs a rebuild.
    f->addDouble("Плотность ρ0 *", 1, 20000, 50, 0, [](const Snap& s) { return s.particleParams.restDensity; },
                 [](Simulation& s, double v) { s.particles.params.restDensity = float(v); s.reset(); }, {}, "кг/м³");
    f->addInt("Итерации несжимаемости", 1, 30, [](const Snap& s) { return s.particleParams.solverIterations; },
              [](Simulation& s, int v) { s.particles.params.solverIterations = v; });
    f->addInt("Подшагов на кадр", 1, 20, [](const Snap& s) { return s.particleParams.substeps; },
              [](Simulation& s, int v) { s.particles.params.substeps = v; });
    f->addDouble("Завихренность", 0, 50, 0.5, 2, [](const Snap& s) { return s.particleParams.vorticity; },
                 [](Simulation& s, double v) { s.particles.params.vorticity = float(v); });
    f->addDouble("Поверхностное натяжение", 0, 0.01, 0.0001, 4, [](const Snap& s) { return s.particleParams.tensileK; },
                 [](Simulation& s, double v) { s.particles.params.tensileK = float(v); });
    f->addDouble("Трение о стенки", 0, 1, 0.05, 2, [](const Snap& s) { return s.particleParams.wallFriction; },
                 [](Simulation& s, double v) { s.particles.params.wallFriction = float(v); });
    f->addDouble("Гравитация g", -30, 30, 0.5, 2, [](const Snap& s) { return -s.rigid.gravity.y; },
                 [](Simulation& s, double v) { s.setGravity({0.0f, float(-v), 0.0f}); s.touchParams(); }, "Одна для всей сцены",
                 "м/с²");
    f->addDouble("Скорость струи", 0.1, 20, 0.5, 2, [](const Snap& s) { return s.emitter.speed; },
                 [](Simulation& s, double v) { s.particles.emitter.speed = float(v); }, {}, "м/с");
    f->addDouble("Радиус струи", 0.01, 0.3, 0.01, 3, [](const Snap& s) { return s.emitter.radius; },
                 [](Simulation& s, double v) { s.particles.emitter.radius = float(v); }, {}, "м");
    return f;
}

ParamForm* MainWindow::buildTokamakForm() {
    auto* f = new ParamForm(ctrl_.get());
    f->addDouble("Запас устойчивости q(a) *", 0.2, 8, 0.1, 2, [](const Snap& s) { return s.tokamak.safetyFactorEdge; },
                 [](Simulation& s, double v) { s.tokamak.safetyFactorEdge = float(v); s.reset(); },
                 "Ток плазмы I_p = 2π a² B0 / (μ0 R0 q(a)). Между 2a²/(a²+b²) = 0.4 и 1 шнур скручивается в винт — "
                 "кинк-неустойчивость (Крускал–Шафранов); выше 1 держит натяжение линий, ниже 0.4 — стенка");
    f->addDouble("Тороидальное поле B0 *", 0.5, 20, 0.5, 1, [](const Snap& s) { return s.tokamak.toroidalField * 1000; },
                 [](Simulation& s, double v) { s.tokamak.toroidalField = float(v) / 1000; s.reset(); },
                 "Поле катушек на магнитной оси; B_φ = B0 R0 / R. Скорость Альфвена растёт с ним, шаг по времени падает", "мТл");
    f->addBool("Вертикальное поле Шафранова *", [](const Snap& s) { return s.tokamak.verticalField; },
               [](Simulation& s, bool v) { s.tokamak.verticalField = v; s.reset(); },
               "Держит кольцо от расширения по большому радиусу; без него его держат только токи изображения в стенке");
    f->addDouble("Затравка кинка *", 0, 0.2, 0.01, 2, [](const Snap& s) { return s.tokamak.seedDisplacement; },
                 [](Simulation& s, double v) { s.tokamak.seedDisplacement = float(v); s.reset(); },
                 "Винтовое смещение шнура в начале (m = 1, n = 1), доля малого радиуса a");
    return f;
}

ParamForm* MainWindow::buildGasForm() {
    auto* f = new ParamForm(ctrl_.get());
    f->addDouble("Скорость потока U", 0, 100, 1, 2, [](const Snap& s) { return s.ns.inflowSpeed; },
                 [](Simulation& s, double v) { s.grid.params.inflowSpeed = float(v); },
                 "Скорость на входной границе (и масштаб для Cp/Cd)", "м/с");
    f->addInt("Ячеек по X *", 8, 400, [](const Snap& s) { return s.ns.resolutionX; },
              [](Simulation& s, int v) { s.grid.params.resolutionX = v; s.reset(); }, "Разрешение сетки вдоль X; по Y/Z — пропорционально");
    f->addBool("Источник тепла и дыма (сфера)", [](const Snap& s) { return s.heat.enabled; },
               [](Simulation& s, bool v) { s.grid.source.enabled = v; });
    f->addDouble("Радиус источника", 0.02, 2, 0.02, 2, [](const Snap& s) { return s.heat.radius; },
                 [](Simulation& s, double v) { s.grid.source.radius = float(v); }, {}, "м");
    f->addDouble("Подъёмная сила тепла", -50, 50, 0.5, 2, [](const Snap& s) { return s.ns.heatBuoyancy; },
                 [](Simulation& s, double v) { s.grid.params.heatBuoyancy = float(v); }, "Буссинеск: ускорение на единицу T", "м/с²");
    f->addBool("Трение о поверхность", [](const Snap& s) { return s.ns.wallFriction; },
               [](Simulation& s, bool v) { s.grid.params.wallFriction = v; },
               "Касательное трение газа о тела (пристеночная функция: Cf пограничного слоя по Шлихтингу)");
    f->addBool("Газ действует на тела", [](const Snap& s) { return s.gasPushesBodies; },
               [](Simulation& s, bool v) { s.gasPushesBodies = v; },
               "Давление газа и архимедова сила на твёрдые тела. Тела всегда вытесняют газ; "
               "с плотностью воздуха 1.2 кг/м³ обратное влияние слабое — увеличьте плотность среды");

    f->beginAdvanced();
    for (int a = 0; a < 3; ++a) {
        const char* names[3] = {"Источник X", "Источник Y", "Источник Z"};
        f->addDouble(names[a], -20, 20, 0.05, 3, [a](const Snap& s) { return s.heat.center[a]; },
                     [a](Simulation& s, double v) { s.grid.source.center[a] = float(v); }, "Центр сферы-источника", "м");
    }
    // Smoke scenes: relative units (~1); fire: kelvin above the air (the burner is 400 K).
    f->addDouble("Температура источника", 0, 3000, 0.1, 2, [](const Snap& s) { return s.heat.temperature; },
                 [](Simulation& s, double v) { s.grid.source.temperature = float(v); },
                 "Дым: относительные единицы (~1); огонь: кельвины над температурой воздуха");
    f->addDouble("Плотность дыма источника", 0, 5, 0.1, 2, [](const Snap& s) { return s.heat.smoke; },
                 [](Simulation& s, double v) { s.grid.source.smoke = float(v); });
    f->addDouble("Вес дыма", -50, 50, 0.1, 2, [](const Snap& s) { return s.ns.smokeBuoyancy; },
                 [](Simulation& s, double v) { s.grid.params.smokeBuoyancy = float(v); }, {}, "м/с²");
    f->addDouble("Остывание газа", 0, 10, 0.05, 3, [](const Snap& s) { return s.ns.temperatureDissipation; },
                 [](Simulation& s, double v) { s.grid.params.temperatureDissipation = float(v); }, {}, "1/с");
    f->addDouble("Рассеяние дыма", 0, 5, 0.01, 3, [](const Snap& s) { return s.ns.smokeDissipation; },
                 [](Simulation& s, double v) { s.grid.params.smokeDissipation = float(v); }, {}, "1/с");
    f->addDouble("Сохранение вихрей", 0, 20, 0.1, 2, [](const Snap& s) { return s.ns.vorticityConfinement; },
                 [](Simulation& s, double v) { s.grid.params.vorticityConfinement = float(v); },
                 "Vorticity confinement: компенсирует численную диссипацию вихрей");
    f->addDouble("Плотность среды ρ", 0.01, 5000, 0.1, 3, [](const Snap& s) { return s.ns.fluidDensity; },
                 [](Simulation& s, double v) { s.grid.params.fluidDensity = float(v); }, "Воздух 1.225, вода 998", "кг/м³");
    f->addDouble("Кинем. вязкость ν", 0, 1, 1e-5, 6, [](const Snap& s) { return s.ns.kinematicViscosity; },
                 [](Simulation& s, double v) { s.grid.params.kinematicViscosity = float(v); }, "Воздух 1.5e-5 м²/с", "м²/с");
    for (int a = 0; a < 3; ++a) {
        const char* names[3] = {"Размер области X *", "Размер области Y *", "Размер области Z *"};
        f->addDouble(names[a], 0.2, 50, 0.1, 2, [a](const Snap& s) { return s.ns.domainSize[a]; },
                     [a](Simulation& s, double v) { s.grid.params.domainSize[a] = float(v); s.reset(); }, {}, "м");
    }
    const char* faces[6] = {"Граница −X *", "Граница +X *", "Граница −Y *", "Граница +Y *", "Граница −Z *", "Граница +Z *"};
    for (int i = 0; i < 6; ++i)
        f->addCombo(faces[i], kBC, [i](const Snap& s) { return int(s.ns.bc[i]); },
                    [i](Simulation& s, int v) { s.grid.params.bc[i] = BoundaryType(v); s.reset(); });
    f->addDouble("Число Куранта (CFL)", 0.2, 8, 0.1, 2, [](const Snap& s) { return s.ns.cfl; },
                 [](Simulation& s, double v) { s.grid.params.cfl = float(v); });
    f->addInt("Макс. итераций давления", 10, 5000, [](const Snap& s) { return s.ns.maxPressureIterations; },
              [](Simulation& s, int v) { s.grid.params.maxPressureIterations = v; });
    f->addDouble("Точность давления", 1e-7, 1e-1, 1e-5, 7, [](const Snap& s) { return s.ns.pressureTolerance; },
                 [](Simulation& s, double v) { s.grid.params.pressureTolerance = float(v); }, "Относительная невязка PCG");
    f->addBool("Адвекция MacCormack (2-й порядок)", [](const Snap& s) { return s.ns.maccormack; },
               [](Simulation& s, bool v) { s.grid.params.maccormack = v; });
    f->addBool("Дымовые струйки на входе", [](const Snap& s) { return s.ns.smokeRake; },
               [](Simulation& s, bool v) { s.grid.params.smokeRake = v; });
    f->addBool("Cd/Cl по площади в плане", [](const Snap& s) { return s.ns.usePlanformArea; },
               [](Simulation& s, bool v) { s.grid.params.usePlanformArea = v; }, "Иначе — по площади миделя");
    return f;
}

ParamForm* MainWindow::buildBrushForm() {
    auto* f = new ParamForm(ctrl_.get());
    auto* radius = new QDoubleSpinBox;
    radius->setRange(0.02, 2.0);
    radius->setSingleStep(0.01);
    radius->setValue(brush_.radius);
    radius->setSuffix(" м");
    connect(radius, &QDoubleSpinBox::valueChanged, this, [this](double v) {
        brush_.radius = float(v);
        view_->setBrushRadius(float(v));
    });
    f->addRow("Радиус", radius);
    auto* strength = new QDoubleSpinBox;
    strength->setRange(0, 10);
    strength->setSingleStep(0.1);
    strength->setValue(brush_.strength);
    strength->setToolTip("Скорость газа = скорость движения точки под курсором × сила");
    connect(strength, &QDoubleSpinBox::valueChanged, this, [this](double v) { brush_.strength = float(v); });
    f->addRow("Сила", strength);
    auto* smoke = new QDoubleSpinBox;
    smoke->setRange(0, 1);
    smoke->setSingleStep(0.05);
    smoke->setValue(brush_.smoke);
    connect(smoke, &QDoubleSpinBox::valueChanged, this, [this](double v) { brush_.smoke = float(v); });
    f->addRow("Дым", smoke);
    f->beginAdvanced();
    auto* heat = new QDoubleSpinBox;
    heat->setRange(0, 5);
    heat->setSingleStep(0.1);
    heat->setValue(brush_.heat);
    connect(heat, &QDoubleSpinBox::valueChanged, this, [this](double v) { brush_.heat = float(v); });
    f->addRow("Тепло", heat);
    f->addHint("ЛКМ по газу: луч из курсора пересекается с плоскостью сечения; движение точки "
               "пересечения задаёт скорость газа. Камера — зажатое колесо.");
    return f;
}

// Region where the body buttons spawn / throw bodies in the current mode.
static AABB bodyArea(const Simulation& s) {
    if (s.mode() == SimMode::Fluid) return s.particles.domain();
    if (s.mode() == SimMode::WindTunnel) return s.grid.domain();
    return AABB({-2, 0, -2}, {2, 5, 2});
}

ParamForm* MainWindow::buildRigidForm() {
    auto* f = new ParamForm(ctrl_.get());
    auto* row = new QWidget;
    auto* h = new QHBoxLayout(row);
    h->setContentsMargins(0, 0, 0, 0);
    auto* bSphere = new QPushButton("+ Шар");
    auto* bBox = new QPushButton("+ Куб");
    auto* bThrow = new QPushButton("Бросить");
    auto* bBullet = new QPushButton("Пуля");
    bBullet->setToolTip("Маленький шар на 150 м/с — проверка CCD");
    auto* bHull = new QPushButton("+ Многогранник");
    bHull->setToolTip("Выпуклый многогранник (цилиндр/конус/гранёный шар по очереди); столкновения через GJK-EPA");
    auto* bTeapot = new QPushButton("+ Чайник");
    auto* bBunny = new QPushButton("+ Кролик");
    bTeapot->setToolTip("Невыпуклое тело: чайник, разбитый на выпуклые части");
    bBunny->setToolTip("Невыпуклое тело: кролик, разбитый на выпуклые части");
    h->addWidget(bSphere);
    h->addWidget(bBox);
    h->addWidget(bHull);
    h->addWidget(bBullet);
    connect(bBullet, &QPushButton::clicked, this, [this] {
        ctrl_->post([](Simulation& s) {
            AABB d = bodyArea(s);
            Vector3 e = d.extent();
            int i = s.rigid.addSphere({d.lo.x + 0.05f * e.x, d.lo.y + 0.3f * e.y, d.center().z}, 0.03f, 8000.0f, {1.0f, 0.95f, 0.4f});
            s.rigid.bodies()[i].vel = {150.0f, 0.0f, 0.0f};
        });
    });
    h->addWidget(bThrow);
    auto spawnPos = [](Simulation& s) {
        static std::mt19937 rng(42);
        std::uniform_real_distribution<float> U(-0.3f, 0.3f);
        AABB d = bodyArea(s);
        Vector3 c = d.center(), e = d.extent();
        return Vector3(c.x + U(rng) * e.x, d.hi.y - 0.15f * e.y, c.z + U(rng) * e.z);
    };
    connect(bSphere, &QPushButton::clicked, this, [this, spawnPos] {
        ctrl_->post([spawnPos](Simulation& s) {
            s.rigid.addSphere(spawnPos(s), s.mode() == SimMode::Rigid ? 0.2f : s.mode() == SimMode::Fluid ? 0.08f : 0.12f, 500.0f, {0.9f, 0.35f, 0.3f});
        });
    });
    connect(bBox, &QPushButton::clicked, this, [this, spawnPos] {
        ctrl_->post([spawnPos](Simulation& s) {
            float hh = s.mode() == SimMode::Rigid ? 0.18f : s.mode() == SimMode::Fluid ? 0.07f : 0.11f;
            s.rigid.addBox(spawnPos(s), Vector3(hh), Quaternion::fromAxisAngle({1, 1, 0.3f}, 0.8f), 500.0f, {0.3f, 0.7f, 0.9f});
        });
    });
    connect(bHull, &QPushButton::clicked, this, [this, spawnPos] {
        ctrl_->post([spawnPos](Simulation& s) {
            static int k = 0;
            const int kind = k++ % 3; // cylinder, cone, faceted ball in turn
            float sc = s.mode() == SimMode::Rigid ? 1.0f : s.mode() == SimMode::Fluid ? 0.5f : 0.7f;
            TriMesh m = kind == 0   ? primitives::cylinder(0.14f * sc, 0.3f * sc, 12)
                        : kind == 1 ? primitives::cone(0.16f * sc, 0.35f * sc, 12)
                                    : primitives::sphere(0.16f * sc, 8, 5);
            s.rigid.addConvex(m, spawnPos(s), Quaternion::fromAxisAngle({1, 0.3f, 0.2f}, 0.9f), 600.0f, {0.8f, 0.8f, 0.3f});
        });
    });
    connect(bTeapot, &QPushButton::clicked, this, [this, spawnPos] {
        ctrl_->post([spawnPos](Simulation& s) {
            s.rigid.addCompound(teapotShape(), spawnPos(s), Quaternion::fromAxisAngle({0.3f, 1, 0.2f}, 0.7f), 500.0f,
                                {0.8f, 0.55f, 0.85f});
        });
    });
    connect(bBunny, &QPushButton::clicked, this, [this, spawnPos] {
        ctrl_->post([spawnPos](Simulation& s) {
            s.rigid.addCompound(bunnyShape(), spawnPos(s), Quaternion::fromAxisAngle({0, 1, 0}, 0.5f), 500.0f,
                                {0.92f, 0.9f, 0.86f});
        });
    });
    connect(bThrow, &QPushButton::clicked, this, [this] {
        ctrl_->post([](Simulation& s) {
            AABB d = bodyArea(s);
            Vector3 e = d.extent();
            int i = s.rigid.addSphere({d.lo.x + 0.1f * e.x, d.lo.y + 0.7f * e.y, d.center().z},
                                      s.mode() == SimMode::Rigid ? 0.22f : s.mode() == SimMode::Fluid ? 0.07f : 0.12f, 3000.0f, {0.95f, 0.8f, 0.2f});
            s.rigid.bodies()[i].vel = Vector3(1.8f * e.x, 0.3f * e.y, 0.0f);
        });
    });
    f->addRow(row);
    auto* row2 = new QWidget;
    auto* h2 = new QHBoxLayout(row2);
    h2->setContentsMargins(0, 0, 0, 0);
    h2->addWidget(bTeapot);
    h2->addWidget(bBunny);
    f->addRow(row2);
    // Soft bodies and cloth: the particle system runs in every mode (liquid, gas, rigid).
    auto* softRow = new QWidget;
    auto* sh = new QHBoxLayout(softRow);
    sh->setContentsMargins(0, 0, 0, 0);
    auto* bSoftCube = new QPushButton("+ Мягкий куб");
    auto* bSoftBall = new QPushButton("+ Желе-шар");
    auto* bCloth = new QPushButton("+ Ткань");
    bSoftCube->setToolTip("Поролон 150 кг/м³: частицы, форма держится сопоставлением формы (Müller 2005)");
    bSoftBall->setToolTip("Мягкий шар: низкая жёсткость формы");
    bCloth->setToolTip("Свободный лоскут хлопка 0.3 кг/м² (XPBD): падает, ложится на всё, рвётся по нитям");
    auto* bWater = new QPushButton("+ Вода");
    bWater->setToolTip("Объём воды (частицы PBF) сверху: в режиме газа она связана с воздухом");
    sh->addWidget(bWater);
    sh->addWidget(bSoftCube);
    sh->addWidget(bSoftBall);
    sh->addWidget(bCloth);
    connect(bWater, &QPushButton::clicked, this, [this] {
        ctrl_->post([](Simulation& s) {
            AABB d = s.particles.domain();
            Vector3 c = d.center(), e = d.extent();
            s.particles.addBlock(AABB({c.x - 0.12f * e.x, d.hi.y - 0.4f * e.y, c.z - 0.25f * e.z},
                                      {c.x + 0.12f * e.x, d.hi.y - 0.1f * e.y, c.z + 0.25f * e.z}));
        });
    });
    auto top = [](Simulation& s) {
        static std::mt19937 rng(7);
        std::uniform_real_distribution<float> U(-0.3f, 0.3f);
        AABB d = s.particles.domain();
        Vector3 c = d.center(), e = d.extent();
        return Vector3(c.x + U(rng) * e.x, d.hi.y - 0.2f * e.y, c.z + 0.3f * U(rng) * e.z);
    };
    connect(bSoftCube, &QPushButton::clicked, this, [this, top] {
        ctrl_->post([top](Simulation& s) {
            TriMesh m = primitives::box(Vector3(0.08f));
            m.translate(top(s));
            s.particles.addSoftBody(m, 150.0f, 0.4f, {0.3f, 0.75f, 0.95f});
        });
    });
    connect(bSoftBall, &QPushButton::clicked, this, [this, top] {
        ctrl_->post([top](Simulation& s) {
            TriMesh m = primitives::sphere(0.08f, 16, 8);
            m.translate(top(s));
            s.particles.addSoftBody(m, 150.0f, 0.15f, {0.55f, 0.9f, 0.35f});
        });
    });
    connect(bCloth, &QPushButton::clicked, this, [this, top] {
        ctrl_->post([top](Simulation& s) {
            Vector3 p = top(s);
            ClothMaterial cotton; // defaults: cotton, tears at ~4 kN/m
            s.particles.addCloth(p - Vector3(0.2f, 0, 0.2f), {0.4f, 0, 0}, {0, 0, 0.4f}, cotton, 0, {0.9f, 0.85f, 0.75f});
        });
    });
    f->addRow(softRow);

    f->beginAdvanced();
    f->addCombo("Решатель", {"XPBD (Müller 2020, эксперим.)", "Импульсы + ударная волна"},
                [](const Snap& s) { return int(s.rigid.solver); },
                [](Simulation& s, int v) {
                    s.rigid.params.solver = RigidSolver(v);
                    s.rigid.params.substeps = v == 0 ? 30 : 10;
                    s.touchParams();
                },
                "Импульсный: тёплый старт, split impulse, ударная волна, трение манифолда и замок SO(3)");
    f->addBool("Ударная волна (стопки)", [](const Snap& s) { return s.rigid.shockPropagation; },
               [](Simulation& s, bool v) { s.rigid.params.shockPropagation = v; });
    f->addBool("Замок вращения SO(3)", [](const Snap& s) { return s.rigid.rotationalLock; },
               [](Simulation& s, bool v) { s.rigid.params.rotationalLock = v; });
    f->addBool("CCD — непрерывные столкновения", [](const Snap& s) { return s.rigid.ccd; },
               [](Simulation& s, bool v) { s.rigid.params.ccd = v; },
               "Консервативное продвижение на GJK: быстрые тела не пролетают сквозь тонкие стенки");
    f->addBool("Блочный решатель манифолда (LCP)", [](const Snap& s) { return s.rigid.blockSolver; },
               [](Simulation& s, bool v) { s.rigid.params.blockSolver = v; });
    f->addBool("Засыпание островов тел", [](const Snap& s) { return s.rigid.sleeping; },
               [](Simulation& s, bool v) {
                   s.rigid.params.sleeping = v;
                   for (int i = 0; i < int(s.rigid.bodies().size()); ++i) s.rigid.wake(i);
               },
               "Покоящиеся острова (граф контактов) перестают считаться до касания; рисуются темнее");
    f->addInt("Итерации контактов", 1, 100, [](const Snap& s) { return s.rigid.iterations; },
              [](Simulation& s, int v) { s.rigid.params.iterations = v; });
    f->addInt("Подшагов", 1, 60, [](const Snap& s) { return s.rigid.substeps; }, // XPBD uses 30
              [](Simulation& s, int v) { s.rigid.params.substeps = v; });
    f->addDouble("Лин. демпфирование", 0, 5, 0.01, 3, [](const Snap& s) { return s.rigid.linearDamping; },
                 [](Simulation& s, double v) { s.rigid.params.linearDamping = float(v); }, {}, "1/с");
    f->addDouble("Угл. демпфирование", 0, 5, 0.01, 3, [](const Snap& s) { return s.rigid.angularDamping; },
                 [](Simulation& s, double v) { s.rigid.params.angularDamping = float(v); }, {}, "1/с");
    // One gravity for bodies, particles and the flame (Simulation::setGravity); the same row as in
    // the liquid form - touchParams refreshes both.
    f->addDouble("Гравитация g", -30, 30, 0.5, 2, [](const Snap& s) { return -s.rigid.gravity.y; },
                 [](Simulation& s, double v) { s.setGravity({0.0f, float(-v), 0.0f}); s.touchParams(); }, "Одна для всей сцены",
                 "м/с²");
    return f;
}

void MainWindow::buildParameterDock() {
    auto* dock = new QDockWidget("Параметры", this);
    dock->setObjectName("paramsDock");
    auto* scroll = new QScrollArea;
    scroll->setWidgetResizable(true);
    scroll->setFrameShape(QFrame::NoFrame);
    auto* host = new QWidget;
    auto* col = new QVBoxLayout(host);
    col->setContentsMargins(6, 6, 6, 6);
    col->setSpacing(8);

    gasBox_ = addGroup(col, "Газ / поток — сетка Навье–Стокса", buildGasForm());
    tokamakBox_ = addGroup(col, "Токамак", buildTokamakForm());
    brushBox_ = addGroup(col, "Кисть возмущения", buildBrushForm());
    fluidBox_ = addGroup(col, "Жидкость — частицы SPH", buildFluidForm());
    objectBox_ = addGroup(col, "Объект / препятствие", buildObjectForm());
    rigidBox_ = addGroup(col, "Твёрдые тела", buildRigidForm());
    col->addStretch(1);

    scroll->setWidget(host);
    scroll->setMinimumWidth(340);
    dock->setWidget(scroll);
    addDockWidget(Qt::LeftDockWidgetArea, dock);
    menuBar()->findChild<QMenu*>("viewMenu")->addAction(dock->toggleViewAction());
}

// ---------------------------------------------------------------------------
// Visualisation dock
// ---------------------------------------------------------------------------
void MainWindow::buildVisualDock() {
    auto* dock = new QDockWidget("Визуализация", this);
    dock->setObjectName("visDock");
    auto* scroll = new QScrollArea;
    scroll->setWidgetResizable(true);
    scroll->setFrameShape(QFrame::NoFrame);
    auto* host = new QWidget;
    auto* col = new QVBoxLayout(host);
    col->setContentsMargins(6, 6, 6, 6);

    {
        auto* f = new ParamForm(ctrl_.get());
        f->addCombo("Поле", kFields, [](const Snap& s) { return int(s.vis.sliceField); },
                    [](Simulation& s, int v) { s.vis.sliceField = GridField(v); });
        f->addCombo("Сечение ⟂ оси", {"X", "Y", "Z"}, [](const Snap& s) { return s.vis.sliceAxis; },
                    [](Simulation& s, int v) { s.vis.sliceAxis = v; });
        f->addSlider("Положение", 0, 1, 200, [](const Snap& s) { return s.vis.slicePosition; },
                     [](Simulation& s, double v) { s.vis.slicePosition = float(v); });
        f->addBool("Цветное сечение поля", [](const Snap& s) { return s.vis.showSlice; },
                   [](Simulation& s, bool v) { s.vis.showSlice = v; });
        f->addBool("Дым (объёмный рендер)", [](const Snap& s) { return s.vis.showSmoke; },
                   [](Simulation& s, bool v) { s.vis.showSmoke = v; });
        f->addBool("Силовые линии магнитного поля", [](const Snap& s) { return s.vis.showFieldLines; },
                   [](Simulation& s, bool v) { s.vis.showFieldLines = v; });
        f->addCombo("Сетка вокселей", {"Скрыть", "В плоскости сечения", "Вся 3D-сетка"},
                    [](const Snap& s) { return s.vis.gridDisplay; }, [](Simulation& s, int v) { s.vis.gridDisplay = v; },
                    "Границы ячеек расчётной сетки");
        f->addCombo("Векторы скорости", {"Скрыть", "В плоскости сечения", "Во всём объёме"},
                    [](const Snap& s) { return s.vis.vectorDisplay; }, [](Simulation& s, int v) { s.vis.vectorDisplay = v; },
                    "Стрелка в центре ячейки: скорость, усреднённая по граням MAC-ячейки");

        f->beginAdvanced();
        f->addInt("Прореживание векторов", 1, 16, [](const Snap& s) { return s.vis.vectorStride; },
                  [](Simulation& s, int v) { s.vis.vectorStride = v; }, "Каждая n-я ячейка", "яч.");
        f->addDouble("Длина векторов", 0.1, 10, 0.1, 2, [](const Snap& s) { return s.vis.vectorScale; },
                     [](Simulation& s, double v) { s.vis.vectorScale = float(v); });
        f->addBool("Векторы только в дыму", [](const Snap& s) { return s.vis.vectorsWhereSmoke; },
                   [](Simulation& s, bool v) { s.vis.vectorsWhereSmoke = v; });
        f->addBool("Линии тока", [](const Snap& s) { return s.vis.showStreamlines; },
                   [](Simulation& s, bool v) { s.vis.showStreamlines = v; });
        f->addInt("Густота линий тока", 2, 40, [](const Snap& s) { return s.vis.streamlineSeeds; },
                  [](Simulation& s, int v) { s.vis.streamlineSeeds = v; });
        f->addBool("Cp на поверхности тела", [](const Snap& s) { return s.vis.surfacePressure; },
                   [](Simulation& s, bool v) { s.vis.surfacePressure = v; });
        f->addBool("Автодиапазон шкалы", [](const Snap& s) { return s.vis.autoRange; },
                   [](Simulation& s, bool v) { s.vis.autoRange = v; });
        f->addDouble("Минимум шкалы", -1e6, 1e6, 0.1, 3, [](const Snap& s) { return s.vis.rangeMin; },
                     [](Simulation& s, double v) { s.vis.rangeMin = float(v); s.vis.autoRange = false; s.touchParams(); });
        f->addDouble("Максимум шкалы", -1e6, 1e6, 0.1, 3, [](const Snap& s) { return s.vis.rangeMax; },
                     [](Simulation& s, double v) { s.vis.rangeMax = float(v); s.vis.autoRange = false; s.touchParams(); });
        auto* smoke = new QDoubleSpinBox;
        smoke->setRange(0.5, 100);
        smoke->setValue(16);
        connect(smoke, &QDoubleSpinBox::valueChanged, view_, [this](double v) { view_->setSmokeDensity(float(v)); });
        f->addRow("Плотность дыма (рендер)", smoke);
        auto* alpha = new QDoubleSpinBox;
        alpha->setRange(0.1, 1);
        alpha->setSingleStep(0.05);
        alpha->setValue(0.92);
        connect(alpha, &QDoubleSpinBox::valueChanged, view_, [this](double v) { view_->setSliceOpacity(float(v)); });
        f->addRow("Непрозрачность сечения", alpha);
        fieldBox_ = addGroup(col, "Поле течения", f);
    }
    {
        auto* f = new ParamForm(ctrl_.get());
        f->addCombo("Цвет частиц", {"Скорость", "Плотность", "Однотонный"},
                    [](const Snap& s) { return int(s.vis.particleColoring); },
                    [](Simulation& s, int v) { s.vis.particleColoring = ParticleColoring(v); });
        f->beginAdvanced();
        auto* pscale = new QDoubleSpinBox;
        pscale->setRange(0.2, 3.0);
        pscale->setSingleStep(0.1);
        pscale->setValue(1.0);
        connect(pscale, &QDoubleSpinBox::valueChanged, view_, [this](double v) { view_->setParticleScale(float(v)); });
        f->addRow("Размер частиц", pscale);
        particleBox_ = addGroup(col, "Частицы", f);
    }
    {
        auto* f = new ParamForm(ctrl_.get());
        auto* cmap = new QComboBox;
        cmap->addItems({"Turbo", "Viridis", "Холодный–тёплый", "Оттенки серого"});
        connect(cmap, &QComboBox::currentIndexChanged, view_, &Viewport::setColormap);
        f->addRow("Палитра", cmap);
        auto* parts = new QCheckBox("Выпуклые части тел");
        parts->setToolTip("Показывать невыпуклые тела как набор выпуклых оболочек, с которыми работает коллизия");
        connect(parts, &QCheckBox::toggled, view_, &Viewport::setShowConvexParts);
        f->addRow(parts);
        // Liquid exists in every mode (pool scenes, water in the gas): the switch is always shown.
        f->addBool("Поверхность воды (шейдер)", [](const Snap& s) { return s.vis.liquidSurface; },
                   [](Simulation& s, bool v) { s.vis.liquidSurface = v; });
        f->beginAdvanced();
        auto* dom = new QCheckBox("Границы области");
        dom->setChecked(true);
        connect(dom, &QCheckBox::toggled, view_, &Viewport::setShowDomain);
        f->addRow(dom);
        auto* floor = new QCheckBox("Сетка пола");
        floor->setChecked(true);
        connect(floor, &QCheckBox::toggled, view_, &Viewport::setShowFloor);
        f->addRow(floor);
        addGroup(col, "Отображение", f);
    }
    col->addStretch(1);
    scroll->setWidget(host);
    scroll->setMinimumWidth(290);
    dock->setWidget(scroll);
    addDockWidget(Qt::RightDockWidgetArea, dock);
    menuBar()->findChild<QMenu*>("viewMenu")->addAction(dock->toggleViewAction());
}

// ---------------------------------------------------------------------------
// Results dock: numbers on the left, one chart per quantity on the right
// ---------------------------------------------------------------------------
void MainWindow::buildResultsDock() {
    auto* dock = new QDockWidget("Результаты", this);
    dock->setObjectName("resultsDock");

    info_ = new QTableWidget(0, 2);
    info_->setHorizontalHeaderLabels({"Величина", "Значение"});
    info_->horizontalHeader()->setSectionResizeMode(0, QHeaderView::ResizeToContents);
    info_->horizontalHeader()->setStretchLastSection(true);
    info_->verticalHeader()->setVisible(false);
    info_->verticalHeader()->setDefaultSectionSize(22);
    info_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    info_->setSelectionMode(QAbstractItemView::NoSelection);
    info_->setShowGrid(false);
    info_->setAlternatingRowColors(true);

    plots_ = new PlotPanel;
    auto* plotHost = new QWidget;
    auto* pl = new QVBoxLayout(plotHost);
    pl->setContentsMargins(0, 0, 0, 0);
    pl->addWidget(plots_, 1);
    auto* bar = new QHBoxLayout;
    bar->addStretch(1);
    auto* clr = new QPushButton("Очистить");
    connect(clr, &QPushButton::clicked, plots_, &PlotPanel::clear);
    auto* csv = new QPushButton("Экспорт CSV…");
    connect(csv, &QPushButton::clicked, this, &MainWindow::exportCsv);
    bar->addWidget(clr);
    bar->addWidget(csv);
    pl->addLayout(bar);

    auto* split = new QSplitter(Qt::Horizontal);
    split->addWidget(info_);
    split->addWidget(plotHost);
    split->setStretchFactor(0, 0);
    split->setStretchFactor(1, 1);
    split->setSizes({340, 1000});
    dock->setWidget(split);
    addDockWidget(Qt::BottomDockWidgetArea, dock);
    menuBar()->findChild<QMenu*>("viewMenu")->addAction(dock->toggleViewAction());
}

// ---------------------------------------------------------------------------
// Runtime
// ---------------------------------------------------------------------------
void MainWindow::onDisturbance(Vector3 pos, Vector3 vel) {
    Disturbance d;
    d.center = pos;
    d.radius = brush_.radius;
    d.velocity = vel * brush_.strength;
    d.velocityBlend = rf::length(vel) > 1e-6f ? 1.0f : 0.0f; // a click adds smoke without stopping the flow
    d.smoke = brush_.smoke;
    d.heat = brush_.heat;
    ctrl_->post([d](Simulation& s) {
        if (s.mode() == SimMode::WindTunnel) s.grid.applyDisturbance(d);
    });
    if (!ctrl_->isRunning()) ctrl_->requestStep(); // make the effect visible even when paused
}

void MainWindow::poll() {
    uint64_t serial = ctrl_->serial();
    if (serial == lastSerial_) return;
    lastSerial_ = serial;
    if (auto s = ctrl_->snapshot()) onSnapshot(s);
}

void MainWindow::updateModeVisibility(SimMode m) {
    gasBox_->setVisible(m == SimMode::WindTunnel);
    brushBox_->setVisible(m == SimMode::WindTunnel);
    fluidBox_->setVisible(m == SimMode::Fluid);
    rigidBox_->setVisible(true); // bodies can be added in every mode (fluid, gas, rigid)
    fieldBox_->setVisible(m == SimMode::WindTunnel);
    particleBox_->setVisible(m == SimMode::Fluid);
}

void MainWindow::updateInfo(const Snap& s) {
    info_->setRowCount(int(s.info.size()));
    for (int r = 0; r < int(s.info.size()); ++r)
        for (int c = 0; c < 2; ++c) {
            QTableWidgetItem* it = info_->item(r, c);
            if (!it) info_->setItem(r, c, it = new QTableWidgetItem);
            it->setText(QString::fromStdString(c == 0 ? s.info[r].first : s.info[r].second));
        }
    status_->setText(QString("%1  ·  кадр %2  ·  t = %3 с  ·  %4")
                         .arg(QString::fromStdString(presetName(s.preset)))
                         .arg(s.frame)
                         .arg(s.time, 0, 'f', 3)
                         .arg(ctrl_->isRunning() ? "идёт расчёт" : "пауза"));
}

void MainWindow::onSnapshot(const std::shared_ptr<const Snap>& s) {
    view_->setSnapshot(s);

    if (s->paramsVersion != lastParamsVersion_) {
        lastParamsVersion_ = s->paramsVersion;
        for (ParamForm* f : forms_) f->refresh(*s);
        updateModeVisibility(s->mode);
        tokamakBox_->setVisible(s->preset == Preset::Tokamak);
        int idx = presetCombo_->findData(int(s->preset));
        if (idx >= 0) presetCombo_->setCurrentIndex(idx);
    }
    if (s->preset != lastPreset_ || s->frame < lastFrame_) {
        plots_->clear();
        lastPreset_ = s->preset;
    }
    if (s->frame != lastFrame_ && s->frame > 0) plots_->append(s->time, s->plots);
    lastFrame_ = s->frame;

    qint64 now = QDateTime::currentMSecsSinceEpoch();
    if (now - lastInfoUpdate_ > 120) {
        lastInfoUpdate_ = now;
        updateInfo(*s);
    }

    if (auto_.active && int(s->preset) == auto_.preset && int(s->frame) >= auto_.frames && ++auto_.settle == 1) {
        ctrl_->setRunning(false);
        QTimer::singleShot(300, this, [this] {
            if (!auto_.shot.isEmpty()) view_->grabFramebuffer().save(auto_.shot);
            if (!auto_.csv.isEmpty()) plots_->writeCsv(auto_.csv);
            qApp->quit();
        });
    }
}

void MainWindow::runAutomation(int preset, int frames, const QString& shot, const QString& csv) {
    preset = std::clamp(preset, 0, int(Preset::Count) - 1);
    ctrl_->post([preset](Simulation& s) { s.loadPreset(Preset(preset)); });
    auto_.active = !shot.isEmpty() || !csv.isEmpty();
    if (!auto_.active) {
        if (!ctrl_->isRunning()) togglePlay(); // opened from the command line: start simulating
        return;
    }
    auto_.preset = preset;
    auto_.frames = std::max(1, frames);
    auto_.shot = shot;
    auto_.csv = csv;
    ctrl_->setRealtimeLimit(false);
    ctrl_->setRunning(true);
}

// ---------------------------------------------------------------------------
// File operations
// ---------------------------------------------------------------------------
void MainWindow::importMesh() {
    QString path = QFileDialog::getOpenFileName(this, "Импорт 3D-модели", QString(), "3D-модели (*.obj *.stl);;Все файлы (*)");
    if (path.isEmpty()) return;
    std::string p = path.toStdString();
    QPointer<MainWindow> self(this);
    ctrl_->post([p, self](Simulation& s) {
        std::string err;
        if (s.loadCustomMesh(p, err)) return;
        QString msg = QString::fromStdString(err);
        QMetaObject::invokeMethod(qApp, [self, msg] {
            if (self) QMessageBox::warning(self, "Импорт модели", "Не удалось загрузить модель:\n" + msg);
        }, Qt::QueuedConnection);
    });
}

void MainWindow::exportSurfaceLoads() {
    QString path = QFileDialog::getSaveFileName(this, "Нагрузки по полигонам", "surface_loads.csv", "CSV (*.csv)");
    if (path.isEmpty()) return;
    std::string p = path.toStdString();
    QPointer<MainWindow> self(this);
    ctrl_->post([p, path, self](Simulation& s) {
        bool empty = s.surfaceLoads().triangles.empty();
        bool ok = !empty && saveSurfaceLoadsCsv(p, s.surfaceLoads());
        QMetaObject::invokeMethod(qApp, [self, ok, empty, path] {
            if (!self) return;
            if (empty) QMessageBox::information(self, "Экспорт", "Нет тела в аэротрубе (или расчёт не запущен).");
            else if (!ok) QMessageBox::warning(self, "Экспорт", "Не удалось записать файл.");
            else self->statusBar()->showMessage("Сохранено: " + path, 4000);
        }, Qt::QueuedConnection);
    });
}

void MainWindow::exportCsv() {
    if (plots_->empty()) {
        QMessageBox::information(this, "Экспорт", "Нет данных: запустите расчёт.");
        return;
    }
    QString path = QFileDialog::getSaveFileName(this, "Экспорт графиков", "results.csv", "CSV (*.csv)");
    if (path.isEmpty()) return;
    if (!plots_->writeCsv(path)) QMessageBox::warning(this, "Экспорт", "Не удалось записать файл.");
    else statusBar()->showMessage("Сохранено: " + path, 4000);
}

void MainWindow::screenshot() {
    QImage img = view_->grabFramebuffer();
    QString def = QString("realflow_%1.png").arg(QDateTime::currentDateTime().toString("yyyyMMdd_hhmmss"));
    QString path = QFileDialog::getSaveFileName(this, "Сохранить скриншот", def, "PNG (*.png)");
    if (!path.isEmpty() && img.save(path)) statusBar()->showMessage("Сохранено: " + path, 4000);
}

void MainWindow::showHelp() {
    QMessageBox box(this);
    box.setWindowTitle("Управление и методы");
    box.setTextFormat(Qt::RichText);
    box.setText(
        "<h3>Управление</h3>"
        "<b>Пробел</b> — старт/пауза, <b>S</b> — шаг, <b>R</b> — сброс, <b>F</b> — показать всю сцену.<br>"
        "Камера: <b>зажатое колесо</b> — вращение, <b>Shift+колесо</b> или ПКМ — сдвиг, прокрутка — масштаб.<br>"
        "<b>ЛКМ</b> (в сценах с газом) — возмущение: из точки экрана строится луч "
        "(обратная матрица проекция×вид), он пересекается с плоскостью сечения; "
        "скорость газа = скорость движения этой точки × «Сила».<br><br>"
        "<h3>Сетка Навье–Стокса (газ, дым)</h3>"
        "Разнесённая MAC-сетка: u, v, w на гранях, p, дым, T — в центрах ячеек. Шаг: "
        "источники → адвекция (полу-лагранжева RK2 / MacCormack) → силы (плавучесть Буссинеска, "
        "сохранение вихрей) → вязкость → проекция давления (PCG) → div u ≈ 0.<br>"
        "Качество проекции видно в таблице: «Макс. |div u| после проекции».<br><br>"
        "<h3>Частицы SPH и твёрдые тела</h3>"
        "Position Based Fluids, ядра poly6/spiky, соседи через хеш-сетку, столкновения с моделями через BVH; "
        "твёрдые тела — импульсный решатель контактов, широкая фаза — BVH.");
    box.exec();
}

void MainWindow::showAbout() {
    QMessageBox::about(this, "О программе",
                       "<b>PhysRealFlow 0.2</b><br>Сетка Навье–Стокса, частицы SPH, BVH, твёрдые тела.<br>"
                       "C++17 · Qt 6 · OpenGL 3.3");
}
