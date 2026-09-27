#include "PlotPanel.h"

#include <QAction>
#include <QCursor>
#include <QFile>
#include <QHBoxLayout>
#include <QLabel>
#include <QMenu>
#include <QTextStream>
#include <QToolButton>
#include <QToolTip>
#include <QVBoxLayout>
#include <QtCharts/QChart>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QtCharts/QValueAxis>

#include <algorithm>
#include <cmath>

namespace {
// Categorical slots (dark-surface steps), assigned in fixed order, never cycled by rank.
const QColor kSeries[] = {QColor("#3987e5"), QColor("#d95926"), QColor("#199e70"), QColor("#c98500"),
                          QColor("#8e6ad6"), QColor("#2aa7b8")};
const QColor kSurface(28, 30, 35);
const QColor kGrid(52, 56, 64);
const QColor kAxisText(150, 157, 170);
const QColor kInk(225, 228, 235);
constexpr int kMaxPoints = 6000;
const char* kDefaultChannel = "frame/step ms";

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
    placeholder_ = new QLabel("Графики появятся после запуска расчёта (Пробел). Величины — в меню «Каналы».");
    placeholder_->setAlignment(Qt::AlignCenter);
    placeholder_->setStyleSheet("color: #8a93a3;");
    row_->addWidget(placeholder_);

    // The chooser lives in the results dock's button bar (channelButton()); its menu lists every
    // quantity the engine has reported, checked ones are drawn.
    channelBtn_ = new QToolButton;
    channelBtn_->setText("Каналы…");
    channelBtn_->setToolTip("Какие величины рисовать: всё, что решатели сообщили отладчику (Probe), по имени");
    channelBtn_->setPopupMode(QToolButton::InstantPopup);
    channelMenu_ = new QMenu(channelBtn_);
    channelBtn_->setMenu(channelMenu_);
}

QWidget* PlotPanel::channelButton() const { return channelBtn_; }

void PlotPanel::setScene(const std::string& name) {
    if (!scene_.empty()) selectionByScene_[scene_] = selected_;
    scene_ = name;
    clear();
    auto it = selectionByScene_.find(name);
    selected_ = it != selectionByScene_.end() ? it->second : std::set<std::string>();
    menuNames_.clear(); // the menu is rebuilt from the first frame of the new scene
}

void PlotPanel::append(double t, const std::vector<std::pair<std::string, float>>& plots, const rf::Probe::Snapshot& probe) {
    times_.push_back(t);
    // Record everything: the scene's plots and every probe channel.
    plotNames_.clear();
    auto record = [&](const std::string& name, double v) {
        auto& h = history_[name];
        h.resize(times_.size() - 1, std::nan(""));
        h.push_back(v);
    };
    for (const auto& [name, v] : plots) {
        plotNames_.push_back(name);
        record(name, v);
    }
    for (const rf::Probe::Channel& c : probe.channels) record(c.name, c.value);

    // The menu: rebuilt when a new quantity appears. A first frame with nothing selected picks
    // the scene's own plots, else the frame time.
    std::vector<std::string> names;
    names.reserve(history_.size());
    for (const auto& [name, h] : history_) names.push_back(name);
    if (names != menuNames_) {
        if (selected_.empty()) {
            for (const std::string& n : plotNames_) selected_.insert(n);
            if (selected_.empty() && history_.count(kDefaultChannel)) selected_.insert(kDefaultChannel);
        }
        rebuildMenu(names);
        applySelection();
    }
    if (history_.empty()) return;
    placeholder_->hide();
    for (auto& [name, chart] : charts_) {
        const auto& h = history_[name];
        if (!h.empty() && h.size() == times_.size()) chart->append(t, h.back());
    }
}

void PlotPanel::rebuildMenu(const std::vector<std::string>& names) {
    menuNames_ = names;
    channelMenu_->clear();
    std::string group;
    for (const std::string& name : names) {
        // Names are "part/quantity": a separator between the parts keeps the long list readable.
        const std::string part = name.substr(0, name.find('/'));
        if (part != group && !group.empty()) channelMenu_->addSeparator();
        group = part;
        auto* act = channelMenu_->addAction(QString::fromStdString(name));
        act->setCheckable(true);
        act->setChecked(selected_.count(name) > 0);
        connect(act, &QAction::toggled, this, [this, name](bool on) {
            if (on) selected_.insert(name);
            else selected_.erase(name);
            applySelection();
        });
    }
    channelBtn_->setText(QString("Каналы (%1 из %2)").arg(selected_.size()).arg(names.size()));
}

QColor PlotPanel::nextColor() {
    return kSeries[colorsUsed_++ % (sizeof(kSeries) / sizeof(kSeries[0]))];
}

void PlotPanel::fillFromHistory(TimeSeriesChart* chart, const std::string& name) {
    const auto& h = history_[name];
    for (size_t i = 0; i < h.size() && i < times_.size(); ++i)
        if (std::isfinite(h[i])) chart->append(times_[i], h[i]);
}

void PlotPanel::applySelection() {
    // Drop the charts that are no longer selected, add the newly selected ones (filled from the
    // recorded history so the curve starts at the beginning, not at the click).
    for (auto it = charts_.begin(); it != charts_.end();) {
        if (selected_.count(it->first)) { ++it; continue; }
        row_->removeWidget(it->second);
        it->second->deleteLater();
        it = charts_.erase(it);
    }
    for (const std::string& name : selected_) {
        if (std::any_of(charts_.begin(), charts_.end(), [&](auto& c) { return c.first == name; })) continue;
        auto* c = new TimeSeriesChart(QString::fromStdString(name), nextColor());
        row_->addWidget(c, 1);
        charts_.emplace_back(name, c);
        fillFromHistory(c, name);
    }
    placeholder_->setVisible(charts_.empty());
    channelBtn_->setText(QString("Каналы (%1 из %2)").arg(selected_.size()).arg(menuNames_.size()));
}

void PlotPanel::clear() {
    for (auto& [name, c] : charts_) {
        row_->removeWidget(c);
        c->deleteLater();
    }
    charts_.clear();
    times_.clear();
    history_.clear();
    plotNames_.clear();
    colorsUsed_ = 0;
    placeholder_->show();
}

bool PlotPanel::writeCsv(const QString& path) const {
    QFile f(path);
    if (!f.open(QIODevice::WriteOnly | QIODevice::Text)) return false;
    // Columns: the scene's plots first, then every probe channel by name.
    std::vector<std::string> columns = plotNames_;
    for (const auto& [name, h] : history_)
        if (std::find(columns.begin(), columns.end(), name) == columns.end()) columns.push_back(name);
    QTextStream out(&f);
    out << "t";
    for (const std::string& name : columns) out << ",\"" << QString::fromStdString(name) << "\"";
    out << "\n";
    for (size_t i = 0; i < times_.size(); ++i) {
        out << QString::number(times_[i], 'g', 8);
        for (const std::string& name : columns) {
            const auto& h = history_.at(name);
            out << ",";
            if (i < h.size() && std::isfinite(h[i])) out << QString::number(h[i], 'g', 8);
        }
        out << "\n";
    }
    return true;
}
