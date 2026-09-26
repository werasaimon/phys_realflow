#pragma once
// Time-series plots as small multiples: one chart per quantity (never two y-scales on one chart),
// titled with the quantity, its unit and the current value; hover shows the exact point.

#include <QColor>
#include <QWidget>

#include <map>
#include <string>
#include <utility>
#include <vector>

class QChart;
class QChartView;
class QHBoxLayout;
class QLabel;
class QLineSeries;
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
    void append(double t, const std::vector<std::pair<std::string, float>>& values);
    void clear();
    bool empty() const { return times_.empty(); }
    bool writeCsv(const QString& path) const;

private:
    QHBoxLayout* row_;
    QLabel* placeholder_;
    std::vector<std::pair<std::string, TimeSeriesChart*>> charts_; // creation order = colour order
    std::vector<double> times_;
    std::map<std::string, std::vector<double>> history_;
};
