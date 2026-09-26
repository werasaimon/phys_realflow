#include "PlotPanel.h"

#include <QCursor>
#include <QFile>
#include <QHBoxLayout>
#include <QLabel>
#include <QTextStream>
#include <QToolTip>
#include <QVBoxLayout>
#include <QtCharts/QChart>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QtCharts/QValueAxis>

#include <cmath>

namespace {
// Categorical slots (dark-surface steps), assigned in fixed order, never cycled by rank.
const QColor kSeries[] = {QColor("#3987e5"), QColor("#d95926"), QColor("#199e70"), QColor("#c98500")};
const QColor kSurface(28, 30, 35);
const QColor kGrid(52, 56, 64);
const QColor kAxisText(150, 157, 170);
const QColor kInk(225, 228, 235);
constexpr int kMaxPoints = 6000;

QString formatValue(double v) {
    double a = std::fabs(v);
    if (a != 0 && (a < 1e-3 || a >= 1e5)) return QString::number(v, 'e', 3);
    return QString::number(v, 'g', 4);
}
} // namespace

// ---------------------------------------------------------------------------
TimeSeriesChart::TimeSeriesChart(const QString& name, const QColor& color, QWidget* parent)
    : QWidget(parent), name_(name) {
    auto* col = new QVBoxLayout(this);
    col->setContentsMargins(0, 0, 0, 0);
    col->setSpacing(0);

    title_ = new QLabel;
    title_->setStyleSheet(QString("color: %1; font-weight: 600; padding: 4px 8px;").arg(kInk.name()));
    col->addWidget(title_);

    chart_ = new QChart;
    chart_->setBackgroundBrush(kSurface);
    chart_->setBackgroundRoundness(0);
    chart_->setMargins(QMargins(2, 2, 6, 2));
    chart_->legend()->hide(); // single series: the title names it
    series_ = new QLineSeries;
    QPen pen(color, 2);
    pen.setCapStyle(Qt::RoundCap);
    series_->setPen(pen);
    chart_->addSeries(series_);

    auto styleAxis = [](QValueAxis* a) {
        a->setLabelsColor(kAxisText);
        a->setGridLineColor(kGrid);
        a->setLinePenColor(kGrid);
        a->setMinorGridLineVisible(false);
        QFont f = a->labelsFont();
        f.setPointSizeF(8);
        a->setLabelsFont(f);
    };
    ax_ = new QValueAxis;
    ax_->setTitleText("t, с");
    ax_->setTitleBrush(kAxisText);
    ax_->setLabelFormat("%.2g");
    ax_->setTickCount(5);
    styleAxis(ax_);
    ay_ = new QValueAxis;
    ay_->setLabelFormat("%.3g");
    ay_->setTickCount(5);
    styleAxis(ay_);
    chart_->addAxis(ax_, Qt::AlignBottom);
    chart_->addAxis(ay_, Qt::AlignLeft);
    series_->attachAxis(ax_);
    series_->attachAxis(ay_);

    view_ = new QChartView(chart_);
    view_->setRenderHint(QPainter::Antialiasing);
    view_->setMinimumSize(260, 150);
    col->addWidget(view_, 1);

    // Hover: exact value at the nearest sample.
    connect(series_, &QLineSeries::hovered, this, [this](const QPointF& p, bool state) {
        if (!state) { QToolTip::hideText(); return; }
        QToolTip::showText(QCursor::pos(), QString("%1\nt = %2 с\n%3").arg(name_, formatValue(p.x()), formatValue(p.y())));
    });
    title_->setText(name_);
}

void TimeSeriesChart::append(double t, double v) {
    if (!std::isfinite(v)) return;
    points_.append(QPointF(t, v));
    if (points_.size() > kMaxPoints) points_.remove(0, points_.size() - kMaxPoints);
    series_->append(t, v);
    if (series_->count() > kMaxPoints) series_->removePoints(0, series_->count() - kMaxPoints);
    title_->setText(QString("%1  <span style='color:#9aa3b2;font-weight:400'>·  %2</span>").arg(name_, formatValue(v)));
    if (++sinceAxisUpdate_ >= 5 || points_.size() < 10) updateAxes();
}

void TimeSeriesChart::updateAxes() {
    sinceAxisUpdate_ = 0;
    if (points_.isEmpty()) return;
    double t0 = points_.front().x(), t1 = points_.back().x();
    double y0 = points_.front().y(), y1 = y0;
    for (const QPointF& p : points_) {
        y0 = std::min(y0, p.y());
        y1 = std::max(y1, p.y());
    }
    if (t1 <= t0) t1 = t0 + 1e-3;
    double pad = std::max((y1 - y0) * 0.1, std::max(std::fabs(y1), 1e-9) * 0.02);
    ax_->setRange(t0, t1);
    ay_->setRange(y0 - pad, y1 + pad);
}

void TimeSeriesChart::clear() {
    points_.clear();
    series_->clear();
    title_->setText(name_);
}

// ---------------------------------------------------------------------------
PlotPanel::PlotPanel(QWidget* parent) : QWidget(parent) {
    row_ = new QHBoxLayout(this);
    row_->setContentsMargins(0, 0, 0, 0);
    row_->setSpacing(6);
    placeholder_ = new QLabel("Графики появятся после запуска расчёта (Пробел).");
    placeholder_->setAlignment(Qt::AlignCenter);
    placeholder_->setStyleSheet("color: #8a93a3;");
    row_->addWidget(placeholder_);
}

void PlotPanel::append(double t, const std::vector<std::pair<std::string, float>>& values) {
    if (values.empty()) return;
    placeholder_->hide();
    times_.push_back(t);
    for (const auto& [name, v] : values) {
        auto it = std::find_if(charts_.begin(), charts_.end(), [&](auto& c) { return c.first == name; });
        if (it == charts_.end()) {
            QColor color = kSeries[charts_.size() % (sizeof(kSeries) / sizeof(kSeries[0]))];
            auto* c = new TimeSeriesChart(QString::fromStdString(name), color);
            row_->addWidget(c, 1);
            charts_.emplace_back(name, c);
            it = charts_.end() - 1;
        }
        it->second->append(t, v);
        auto& h = history_[name];
        h.resize(times_.size() - 1, std::nan(""));
        h.push_back(v);
    }
}

void PlotPanel::clear() {
    for (auto& [name, c] : charts_) {
        row_->removeWidget(c);
        c->deleteLater();
    }
    charts_.clear();
    times_.clear();
    history_.clear();
    placeholder_->show();
}

bool PlotPanel::writeCsv(const QString& path) const {
    QFile f(path);
    if (!f.open(QIODevice::WriteOnly | QIODevice::Text)) return false;
    QTextStream out(&f);
    out << "t";
    for (auto& [name, c] : charts_) out << ",\"" << QString::fromStdString(name) << "\"";
    out << "\n";
    for (size_t i = 0; i < times_.size(); ++i) {
        out << QString::number(times_[i], 'g', 8);
        for (auto& [name, c] : charts_) {
            const auto& h = history_.at(name);
            out << ",";
            if (i < h.size() && std::isfinite(h[i])) out << QString::number(h[i], 'g', 8);
        }
        out << "\n";
    }
    return true;
}
