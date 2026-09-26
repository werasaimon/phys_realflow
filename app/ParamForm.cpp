#include "ParamForm.h"

#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QFormLayout>
#include <QLabel>
#include <QLineEdit>
#include <QSignalBlocker>
#include <QSlider>
#include <QSpinBox>
#include <QToolButton>
#include <QVBoxLayout>

#include <cmath>

ParamForm::ParamForm(SimController* ctrl, QWidget* parent) : QWidget(parent), ctrl_(ctrl) {
    root_ = new QVBoxLayout(this);
    root_->setContentsMargins(6, 6, 6, 6);
    root_->setSpacing(4);
    form_ = makeForm();
    root_->addLayout(form_);
}

QFormLayout* ParamForm::makeForm() {
    auto* f = new QFormLayout;
    f->setContentsMargins(0, 0, 0, 0);
    f->setHorizontalSpacing(10);
    f->setVerticalSpacing(5);
    f->setFieldGrowthPolicy(QFormLayout::AllNonFixedFieldsGrow);
    return f;
}

void ParamForm::beginAdvanced(const QString& title) {
    auto* toggle = new QToolButton;
    toggle->setText(title);
    toggle->setCheckable(true);
    toggle->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    toggle->setArrowType(Qt::RightArrow);
    toggle->setAutoRaise(true);
    toggle->setStyleSheet("QToolButton { color: #8a93a3; border: none; padding: 2px 0; }");
    auto* box = new QWidget;
    box->setVisible(false);
    form_ = makeForm();
    box->setLayout(form_);
    connect(toggle, &QToolButton::toggled, box, [toggle, box](bool on) {
        box->setVisible(on);
        toggle->setArrowType(on ? Qt::DownArrow : Qt::RightArrow);
    });
    root_->addWidget(toggle);
    root_->addWidget(box);
}

QDoubleSpinBox* ParamForm::addDouble(const QString& label, double min, double max, double step, int decimals,
                                     Get<double> get, Set<double> set, const QString& tip, const QString& suffix) {
    auto* w = new QDoubleSpinBox;
    w->setRange(min, max);
    w->setSingleStep(step);
    w->setDecimals(decimals);
    w->setKeyboardTracking(false);
    w->setAccelerated(true);
    if (!suffix.isEmpty()) w->setSuffix(" " + suffix);
    w->setToolTip(tip);
    connect(w, &QDoubleSpinBox::valueChanged, this, [this, set](double v) { ctrl_->post([set, v](rf::Simulation& s) { set(s, v); }); });
    refreshers_.push_back([w, get](const Snap& s) { w->setValue(get(s)); });
    auto* l = new QLabel(label);
    l->setToolTip(tip);
    form_->addRow(l, w);
    return w;
}

QSpinBox* ParamForm::addInt(const QString& label, int min, int max, Get<int> get, Set<int> set, const QString& tip,
                            const QString& suffix) {
    auto* w = new QSpinBox;
    w->setRange(min, max);
    w->setKeyboardTracking(false);
    if (!suffix.isEmpty()) w->setSuffix(" " + suffix);
    w->setToolTip(tip);
    connect(w, &QSpinBox::valueChanged, this, [this, set](int v) { ctrl_->post([set, v](rf::Simulation& s) { set(s, v); }); });
    refreshers_.push_back([w, get](const Snap& s) { w->setValue(get(s)); });
    auto* l = new QLabel(label);
    l->setToolTip(tip);
    form_->addRow(l, w);
    return w;
}

QCheckBox* ParamForm::addBool(const QString& label, Get<bool> get, Set<bool> set, const QString& tip) {
    auto* w = new QCheckBox(label);
    w->setToolTip(tip);
    connect(w, &QCheckBox::toggled, this, [this, set](bool v) { ctrl_->post([set, v](rf::Simulation& s) { set(s, v); }); });
    refreshers_.push_back([w, get](const Snap& s) { w->setChecked(get(s)); });
    form_->addRow(w);
    return w;
}

QComboBox* ParamForm::addCombo(const QString& label, const QStringList& items, Get<int> get, Set<int> set,
                               const QString& tip) {
    auto* w = new QComboBox;
    w->addItems(items);
    w->setToolTip(tip);
    connect(w, &QComboBox::currentIndexChanged, this,
            [this, set](int v) { if (v >= 0) ctrl_->post([set, v](rf::Simulation& s) { set(s, v); }); });
    refreshers_.push_back([w, get](const Snap& s) { w->setCurrentIndex(get(s)); });
    auto* l = new QLabel(label);
    l->setToolTip(tip);
    form_->addRow(l, w);
    return w;
}

QSlider* ParamForm::addSlider(const QString& label, double min, double max, int steps, Get<double> get, Set<double> set,
                              const QString& tip) {
    auto* w = new QSlider(Qt::Horizontal);
    w->setRange(0, steps);
    w->setToolTip(tip);
    connect(w, &QSlider::valueChanged, this, [this, set, min, max, steps](int i) {
        double v = min + (max - min) * i / steps;
        ctrl_->post([set, v](rf::Simulation& s) { set(s, v); });
    });
    refreshers_.push_back([w, get, min, max, steps](const Snap& s) {
        w->setValue(int(std::lround((get(s) - min) / (max - min) * steps)));
    });
    auto* l = new QLabel(label);
    l->setToolTip(tip);
    form_->addRow(l, w);
    return w;
}

QLineEdit* ParamForm::addText(const QString& label, Get<QString> get, Set<QString> set, const QString& tip) {
    auto* w = new QLineEdit;
    w->setToolTip(tip);
    connect(w, &QLineEdit::editingFinished, this, [this, w, set]() {
        QString v = w->text();
        ctrl_->post([set, v](rf::Simulation& s) { set(s, v); });
    });
    refreshers_.push_back([w, get](const Snap& s) { if (!w->hasFocus()) w->setText(get(s)); });
    auto* l = new QLabel(label);
    l->setToolTip(tip);
    form_->addRow(l, w);
    return w;
}

void ParamForm::addRow(QWidget* w) { form_->addRow(w); }
void ParamForm::addRow(const QString& label, QWidget* w) { form_->addRow(label, w); }

void ParamForm::addHint(const QString& text) {
    auto* l = new QLabel(text);
    l->setWordWrap(true);
    l->setStyleSheet("color: #8a93a3; font-size: 8.5pt;");
    form_->addRow(l);
}

void ParamForm::refresh(const Snap& s) {
    // Block signals of all children so programmatic updates are not echoed back as commands.
    QList<QWidget*> children = findChildren<QWidget*>();
    std::vector<std::unique_ptr<QSignalBlocker>> blockers;
    blockers.reserve(children.size());
    for (QWidget* c : children) blockers.push_back(std::make_unique<QSignalBlocker>(c));
    for (auto& r : refreshers_) r(s);
}
