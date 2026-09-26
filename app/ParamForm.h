#pragma once
// Declarative parameter form: each row binds a widget to a getter (reads the published snapshot)
// and a setter (posted to the simulation thread).

#include "SimController.h"

#include <QWidget>

#include <functional>
#include <vector>

class QCheckBox;
class QComboBox;
class QDoubleSpinBox;
class QFormLayout;
class QLineEdit;
class QSpinBox;

class ParamForm : public QWidget {
    Q_OBJECT
public:
    using Snap = rf::RenderSnapshot;
    template <class T> using Get = std::function<T(const Snap&)>;
    template <class T> using Set = std::function<void(rf::Simulation&, T)>;

    explicit ParamForm(SimController* ctrl, QWidget* parent = nullptr);

    QDoubleSpinBox* addDouble(const QString& label, double min, double max, double step, int decimals, Get<double> get,
                              Set<double> set, const QString& tip = {}, const QString& suffix = {});
    QSpinBox* addInt(const QString& label, int min, int max, Get<int> get, Set<int> set, const QString& tip = {},
                     const QString& suffix = {});
    QCheckBox* addBool(const QString& label, Get<bool> get, Set<bool> set, const QString& tip = {});
    QComboBox* addCombo(const QString& label, const QStringList& items, Get<int> get, Set<int> set,
                        const QString& tip = {});
    // Slider over [min, max] with `steps` increments.
    class QSlider* addSlider(const QString& label, double min, double max, int steps, Get<double> get, Set<double> set,
                             const QString& tip = {});
    QLineEdit* addText(const QString& label, Get<QString> get, Set<QString> set, const QString& tip = {});
    void addRow(QWidget* w);
    void addRow(const QString& label, QWidget* w);
    void addHint(const QString& text);
    // Rows added after this call go into a collapsed "advanced" section.
    void beginAdvanced(const QString& title = QStringLiteral("Дополнительно"));

    void refresh(const Snap& s);

private:
    SimController* ctrl_;
    QFormLayout* makeForm();

    class QVBoxLayout* root_;
    QFormLayout* form_; // current target (basic or advanced)
    std::vector<std::function<void(const Snap&)>> refreshers_;
};
