#pragma once

#include "SimController.h"

#include <QMainWindow>

#include <memory>
#include <vector>

class ParamForm;
class PlotPanel;
class Viewport;
class QAction;
class QComboBox;
class QGroupBox;
class QLabel;
class QTableWidget;
class QTimer;

class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    MainWindow();
    ~MainWindow() override;

    // Automation: load a preset, simulate `frames` frames, save a screenshot and/or CSV, quit.
    void runAutomation(int preset, int frames, const QString& screenshotPath, const QString& csvPath);

private:
    // UI construction
    void buildActions();
    void buildParameterDock();
    void buildVisualDock();
    void buildResultsDock();
    ParamForm* buildObjectForm();
    ParamForm* buildFluidForm();
    ParamForm* buildGasForm();
    ParamForm* buildBrushForm();
    ParamForm* buildRigidForm();
    QGroupBox* addGroup(class QVBoxLayout* col, const QString& title, ParamForm* form);

    // Runtime
    void poll();
    void onSnapshot(const std::shared_ptr<const rf::RenderSnapshot>& s);
    void updateModeVisibility(rf::SimMode m);
    void updateInfo(const rf::RenderSnapshot& s);
    void togglePlay();
    void onDisturbance(rf::Vector3 pos, rf::Vector3 vel);

    // File operations
    void importMesh();
    void exportCsv();
    void exportSurfaceLoads();
    void screenshot();
    void showHelp();
    void showAbout();

    std::unique_ptr<SimController> ctrl_;
    Viewport* view_ = nullptr;
    QTimer* timer_ = nullptr;
    uint64_t lastSerial_ = 0, lastParamsVersion_ = 0, lastFrame_ = 0;
    rf::Preset lastPreset_ = rf::Preset::Count;
    qint64 lastInfoUpdate_ = 0;

    QComboBox* presetCombo_ = nullptr;
    QAction* playAct_ = nullptr;
    std::vector<ParamForm*> forms_;
    QGroupBox *objectBox_ = nullptr, *fluidBox_ = nullptr, *gasBox_ = nullptr, *brushBox_ = nullptr,
              *rigidBox_ = nullptr, *fieldBox_ = nullptr, *particleBox_ = nullptr;
    QTableWidget* info_ = nullptr;
    PlotPanel* plots_ = nullptr;
    QLabel* status_ = nullptr;

    // Mouse brush (GUI-side settings, applied through Disturbance commands)
    struct Brush {
        float radius = 0.15f;   // m
        float strength = 1.0f;  // multiplier of the cursor speed
        float smoke = 0.6f;     // added per stroke sample
        float heat = 0.0f;
    } brush_;

    struct Automation {
        bool active = false;
        int preset = 0;
        int frames = 0;
        int settle = 0;
        QString shot, csv;
    } auto_;
};
