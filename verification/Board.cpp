// The generated block of the board (docs/09-benchmarks.md): one table per category, the
// significance of the validation set (chi^2, look elsewhere), the state of the blind analysis,
// and the numbers behind every row. Written only between the two markers; the rest of the file
// is read from disk at the moment of writing and kept as it is.
#include "Report.h"

#include "Significance.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>

namespace rf::verify {

static const char* kBegin = "<!-- rf_verify: начало (генерирует rf_verify; правки руками здесь пропадут) -->";
static const char* kEnd = "<!-- rf_verify: конец -->";

static std::string printfString(const char* f, double a, double b = 0) {
    char buf[96];
    std::snprintf(buf, sizeof(buf), f, a, b);
    return buf;
}

std::string regressionMark(const std::vector<HistoryPoint>& earlier, double value, double sigma) {
    if (earlier.empty()) return "первый прогон";
    const HistoryPoint& prev = earlier.back();
    const double delta = value - prev.value, s = std::sqrt(prev.sigma * prev.sigma + sigma * sigma);
    if (!std::isfinite(delta)) return "—";
    if (s > 0) {
        const double n = std::fabs(delta) / s;
        if (n > 3) return (delta > 0 ? "↑ " : "↓ ") + printfString("%.1fσ", n);
        return printfString("≈ (%.1fσ)", n);
    }
    // A deterministic case: equal up to the 10 significant digits the result files keep.
    if (std::fabs(delta) <= 1e-9 * std::max(std::fabs(value), std::fabs(prev.value))) return "= (детерм.)";
    return std::string(delta > 0 ? "↑" : "↓") + " (детерм., Δ = " + num(delta) + ")";
}

static std::string trend(const History& history, const Outcome& o) {
    const auto it = history.find(o.c->id);
    const std::vector<HistoryPoint> none;
    const std::vector<HistoryPoint>& v = it == history.end() ? none : it->second;
    std::string s;
    for (size_t i = v.size() > 2 ? v.size() - 2 : 0; i < v.size(); ++i) s += num(v[i].value) + " → ";
    const double sigma = o.r.runs > 1 ? o.r.spread / std::sqrt(double(o.r.runs)) : 0;
    return s + "**" + num(o.r.value) + "** " + regressionMark(v, o.r.value, sigma);
}

// Text for a table cell: a bare | would end the cell.
static std::string cell(const std::string& s) {
    std::string out;
    for (char ch : s) out += ch == '|' ? std::string("\\|") : std::string(1, ch);
    return out;
}

static std::string sourceCell(const Reference& ref) {
    if (ref.doi_or_url.empty()) return cell(ref.source);
    return "[" + cell(ref.source) + "](" + ref.doi_or_url + ")";
}

static std::string oursCell(const Result& r) {
    std::string s = num(r.value);
    if (r.numericalUncertainty > 0) s += " ± " + num(r.numericalUncertainty);
    if (r.inputUncertainty > 0) s += " ± " + num(r.inputUncertainty) + " (UQ)";
    if (r.runs > 1) s += " (" + std::to_string(r.runs) + " зёрен, σ " + num(r.spread) + ")";
    return s + (r.unit.empty() ? "" : " " + r.unit);
}

static std::string referenceCell(const Outcome& o) {
    if (o.status == BlindStatus::Blinded) return "🔒 запечатан, SHA-256 `" + o.c->commitment.substr(0, 12) + "…`";
    return num(o.ref.value) + (o.ref.uncertainty > 0 ? " ± " + num(o.ref.uncertainty) : "");
}

// A pass whose pull is 2 sigma or more: the bias is resolved - the uncertainty is smaller than the
// bias itself, so the bias is measured, not noise - and it is still within the case's tolerance.
// The verdict stays a pass; the board says both, so a green tick never stands silently next to
// "7.6 sigma".
static bool resolvedBias(const Outcome& o) {
    return o.a.verdict == Verdict::Pass && std::isfinite(o.a.z) && std::fabs(o.a.z) >= 2;
}

// Footnote marks as superscript digits: 1 -> ¹, 12 -> ¹².
static std::string superscript(int n) {
    static const char* digits[10] = {"⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹"};
    std::string s;
    for (char c : std::to_string(n)) s += digits[c - '0'];
    return s;
}

struct Footnotes {
    int count = 0;    // marks given so far (numbered through the whole block)
    std::string text; // the notes of the current section, printed under its table
};

// The note under the table: E in % of the reference, the tolerance, and why z is large anyway.
static std::string biasNote(const Outcome& o, int mark) {
    const Acceptance& acc = o.c->acceptance;
    const std::string unit = o.r.unit.empty() ? "" : " " + o.r.unit;
    const std::string tolerance = acc.relative ? num(100 * acc.tolerance) + " %" : num(acc.tolerance) + unit;
    const bool numericalOnly = o.r.inputUncertainty == 0 && !(o.ref.uncertainty > 0);
    std::string s = "- " + superscript(mark) + " `" + o.c->id + "`: E = " + num(o.a.E) + unit;
    if (std::isfinite(o.a.relativeE)) s += " (" + printfString("%+.2f", 100 * o.a.relativeE) + " % от эталона)";
    s += ", в пределах допуска ±" + tolerance + ". z = " + printfString("%+.1f", o.a.z) + " σ, потому что ";
    s += numericalOnly ? "численная неопределённость u_num = " + num(o.a.uVal) + " (из GCI, стандартной ошибки порядка или ширины брекета)"
                       : "неопределённость валидации u_val = " + num(o.a.uVal);
    return s + " меньше самого смещения |E|: смещение измерено уверенно и при этом допустимо.\n";
}

static std::string statusCell(const Outcome& o, Footnotes& notes) {
    if (!resolvedBias(o)) return verdictSymbol(o.a.verdict);
    ++notes.count;
    notes.text += biasNote(o, notes.count);
    return "✅ смещение разрешено" + superscript(notes.count);
}

static std::string tableRow(const Outcome& o, const History& history, Footnotes& notes) {
    const Case& c = *o.c;
    std::string E = num(o.a.E);
    if (std::isfinite(o.a.relativeE) && o.ref.value != 0) E += " (" + num(100 * o.a.relativeE) + " %)";
    std::string order = "—";
    if (std::isfinite(o.r.observedOrder) || std::isfinite(o.r.theoreticalOrder)) order = num(o.r.observedOrder) + " / " + num(o.r.theoreticalOrder);
    const std::string z = std::isfinite(o.a.z) ? printfString("%+.1f", o.a.z) : "—";
    return "| " + cell(c.title) + " | " + sourceCell(c.ref) + " | " + referenceCell(o) + " | " + oursCell(o.r) + " | " + E + " | " + num(o.a.uVal) +
           " | " + z + " | " + order + " | " + statusCell(o, notes) + " | " + trend(history, o) + " | `" + c.id + "`<br>" +
           statusWord(o.status) + ", " + c.split + (c.test.empty() ? "" : "<br>`" + c.test + "`") + " |\n";
}

static std::string section(const char* title, const char* category, const std::vector<Outcome>& outcomes, const History& history,
                           Footnotes& notes) {
    std::string rows;
    notes.text.clear();
    for (const Outcome& o : outcomes)
        if (o.c->category == category) rows += tableRow(o, history, notes);
    if (rows.empty()) return "";
    return std::string("### ") + title +
           "\n\n| Случай | Эталон и источник | D ± u_D | S ± u_num ± u_input | E = S − D | u_val | z, σ | порядок набл. / теор. | Статус | История | id, слепота, тест |\n"
           "|---|---|---|---|---|---|---|---|---|---|---|\n" + rows + "\n" + (notes.text.empty() ? "" : notes.text + "\n");
}

// chi^2 of the validation pulls and the look-elsewhere reading (Significance.h).
static std::string significance(const std::vector<Outcome>& outcomes) {
    std::vector<std::pair<std::string, double>> pulls;
    int blinded = 0;
    for (const Outcome& o : outcomes) {
        if (o.c->category != "validation") continue;
        if (o.status == BlindStatus::Blinded) ++blinded;
        else pulls.push_back({o.c->id, o.a.z});
    }
    const SetSummary s = summarize(pulls);
    std::string t = "### Значимость набора валидации\n\n";
    if (s.n == 0) return t + "Нет открытых валидационных случаев с конечным z.\n\n";
    char buf[512];
    std::snprintf(buf, sizeof(buf), "χ² = %.2f при ndf = %d (χ²/ndf = %.2f), p = %.3g. Наибольший пул %+.1fσ (`%s`); поправка Бонферрони "
                  "на число случаев N = %d: глобальное p ≤ %.3g. Ожидаемое число |z| > 2 среди N = %d согласных случаев: %.2f, наблюдается %d.",
                  s.chi2, s.n, s.chi2 / s.n, s.pValue, s.largestPull, s.largestId.c_str(), s.n, s.bonferroniP, s.n, s.expectedOver2,
                  s.observedOver2);
    t += buf;
    if (blinded > 0) t += " Слепых валидационных случаев вне суммы: " + std::to_string(blinded) + ".";
    return t + "\n\n";
}

static std::string blindTable(const std::vector<Outcome>& outcomes) {
    std::string rows;
    for (const Outcome& o : outcomes)
        if (o.status != BlindStatus::Open || o.c->split == "holdout")
            rows += "| `" + o.c->id + "` | " + statusWord(o.status) + " | " + o.c->split + " | `" + o.c->commitment + "` |\n";
    if (rows.empty()) return "";
    return "### Слепой анализ\n\nЭталоны этих случаев запечатаны ([процедура](10-verification.md)); здесь только SHA-256 обязательства. "
           "Журнал снятия слепоты: [verification/unblinding-log.md](../verification/unblinding-log.md).\n\n"
           "| id | статус | выборка | SHA-256(value\\|uncertainty\\|salt) |\n|---|---|---|---|\n" + rows + "\n";
}

static std::string details(const std::vector<Outcome>& outcomes, const std::string& imagesDir) {
    std::string s = "### Числа за каждой строкой\n\n";
    for (const Outcome& o : outcomes) {
        const bool bias = resolvedBias(o);
        s += "- `" + o.c->id + "` " + (bias ? "✅ смещение разрешено —" : verdictSymbol(o.a.verdict)) + " " + o.a.wording +
             (bias ? ", но |E| в пределах допуска" : "") + ". " + o.r.detail + ".";
        if (!o.c->ref.note.empty()) s += " Эталон: " + o.c->ref.note + ".";
        if (o.r.convergence.size() >= 2) s += " [график](" + imagesDir + "/" + o.c->id + ".svg)";
        s += "\n";
    }
    return s + "\n";
}

static std::string header(const RunInfo& info, const std::vector<Outcome>& outcomes) {
    int counts[5] = {0, 0, 0, 0, 0}, biases = 0;
    for (const Outcome& o : outcomes) { ++counts[int(o.a.verdict)]; biases += resolvedBias(o); }
    const std::string ofWhich = biases ? " (из них " + std::to_string(biases) + " — «смещение разрешено»)" : "";
    char head[512];
    std::snprintf(head, sizeof(head), "Прогон `rf_verify --%s`: %s, коммит `%s`, %.0f с. Итог: %d ✅%s, %d 🟡, %d ❌, %d ⚠️, %d 🔒 слепых.",
                  info.mode.c_str(), info.utc.c_str(), info.commit.c_str(), info.seconds, counts[0], ofWhich.c_str(), counts[1], counts[2],
                  counts[3], counts[4]);
    return std::string(kBegin) + "\n## Реестр rf_verify (ASME V&V 20, слепой анализ)\n\n" + head + (info.note.empty() ? "" : " " + info.note) + "\n\n" +
           "Как читать ([метод](10-verification.md)): S — наше число, D — эталон, E = S − D, u_val = √(u_num² + u_input² + u_D²), "
           "пул z = E / u_val в σ: |z| < 2 согласие, 2–3 напряжение, 3–5 свидетельство ошибки, ≥ 5 ошибка установлена. Статус — по "
           "допуску случая на |E|, а не на |E| − u_val: большая неопределённость не превращает плохой ответ в зачёт. Статус и пул "
           "отвечают на разные вопросы: допуск — «достаточно ли близко для нашей цели», z — «отличимо ли расхождение от нашей "
           "неопределённости». Поэтому «✅ смещение разрешено» — |E| в пределах допуска, но |z| ≥ 2: неопределённость (для проверки "
           "кода — численная u_num: GCI, стандартная ошибка порядка или ширина брекета) меньше самого смещения, то есть смещение "
           "измерено, а не шум, и при этом допустимо; "
           "сноска под таблицей даёт E в % от эталона. Для проверки кода число — наблюдаемый порядок схемы против теоретического. "
           "«История»: ↑/↓ — изменение больше 3σ разброса между прогонами, ≈ — в пределах.\n\n";
}

std::string boardBlock(const RunInfo& info, const std::vector<Outcome>& outcomes, const History& history, const std::string& imagesDir) {
    std::string s = header(info, outcomes);
    Footnotes notes;
    s += section("Проверка кода (code verification)", "code-verification", outcomes, history, notes);
    s += section("Проверка решения (solution verification)", "solution-verification", outcomes, history, notes);
    s += section("Валидация (validation)", "validation", outcomes, history, notes);
    s += significance(outcomes);
    s += "![Пулы z всех случаев с полосами 1, 2, 3 σ](" + imagesDir + "/pulls.svg)\n\n";
    bool anyValidation = false;
    for (const Outcome& o : outcomes) anyValidation |= o.c->category == "validation" && o.status != BlindStatus::Blinded;
    if (anyValidation) s += "![E ± u_val валидационных случаев](" + imagesDir + "/validation.svg)\n\n";
    return s + blindTable(outcomes) + details(outcomes, imagesDir) + kEnd + "\n";
}

bool updateBoard(const std::string& boardPath, const std::string& block) {
    std::ifstream in(boardPath, std::ios::binary);
    if (!in) return false;
    std::stringstream buf;
    buf << in.rdbuf();
    std::string text = buf.str();
    in.close();
    const size_t b = text.find(kBegin), e = text.find(kEnd);
    if (b != std::string::npos && e != std::string::npos && e > b) {
        size_t after = e + std::string(kEnd).size();
        if (after < text.size() && text[after] == '\n') ++after;
        text = text.substr(0, b) + block + text.substr(after);
    } else { // first run: before the first section, and the rest becomes the hand-kept part
        const size_t first = text.find("\n## ");
        const size_t at = first == std::string::npos ? text.size() : first + 1;
        text = text.substr(0, at) + block + "\n## Строки, которые пока ведутся вручную\n\n"
               "Эти цифры печатают тесты `rf_tests`; в реестр `rf_verify` они ещё не перенесены.\n\n" + text.substr(at);
    }
    std::ofstream out(boardPath, std::ios::binary);
    out << text;
    return bool(out);
}

} // namespace rf::verify
