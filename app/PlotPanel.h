#pragma once
// Time-series plots as small multiples: one chart per quantity (never two y-scales on one chart),
// titled with the quantity, its unit and the current value; hover shows the exact point.
// The quantities come from the engine's Probe (every channel it reports, by name) and from the
// scene's own plots; the "Каналы" menu chooses which of them are drawn, and the CSV export
// writes all of them.

#include "core/Probe.h"

#include <QColor>
#include <QWidget>

#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

class QChart;
class QChartView;
class QHBoxLayout;
class QLabel;
class QLineSeries;
class QMenu;
class QToolButton;
class QValueAxis;

class TimeSeriesChart : public QWidget {
    Q_OBJECT
public:
    TimeSeriesChart(const QString& name, const QColor& color, QWidget* parent = nullptr);
    void append(double t, double v);
    void clear();

private:
    void updateAxes();

    QString name_;
    QLabel* title_;
    QChart* chart_;
    QChartView* view_;
    QLineSeries* series_;
    QValueAxis *ax_, *ay_;
    QList<QPointF> points_;
    int sinceAxisUpdate_ = 0;
};

class PlotPanel : public QWidget {
    Q_OBJECT
public:
    explicit PlotPanel(QWidget* parent = nullptr);
    // A new scene: series are dropped; its remembered channel selection (or the default: the
    // scene's own plots, else the frame time) is restored.
    void setScene(const std::string& name);
    // One frame: the scene's plots and every probe channel are recorded (for the CSV); the
    // selected ones are drawn.
    void append(double t, const std::vector<std::pair<std::string, float>>& plots, const rf::Probe::Snapshot& probe);
    void clear();
    bool empty() const { return times_.empty(); }
    // All recorded quantities, one column each, one row per frame.
    bool writeCsv(const QString& path) const;
    QWidget* channelButton() const;

private:
    void rebuildMenu(const std::vector<std::string>& names);
    void applySelection();
    void fillFromHistory(TimeSeriesChart* chart, const std::string& name);
    QColor nextColor();

    QHBoxLayout* row_;
    QLabel* placeholder_;
    QToolButton* channelBtn_;
    QMenu* channelMenu_;
    std::vector<std::string> menuNames_; // the names the menu was built from
    std::vector<std::pair<std::string, TimeSeriesChart*>> charts_; // in selection order
    std::vector<double> times_;
    std::map<std::string, std::vector<double>> history_; // every recorded quantity
    std::vector<std::string> plotNames_;                  // the scene's own plots (first in the CSV)
    std::set<std::string> selected_;
    std::map<std::string, std::set<std::string>> selectionByScene_;
    std::string scene_;
    int colorsUsed_ = 0;
};
