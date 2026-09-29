// Numbers in files (Format.h): standard streams imbued with the classic locale, so neither the
// program's C locale (printf, strtof) nor its global C++ locale can change the text.
#include "core/Format.h"

#include <iomanip>
#include <locale>
#include <sstream>

namespace rf {

bool parseNumber(const std::string& word, float& out) {
    std::istringstream in(word);
    in.imbue(std::locale::classic());
    in >> out;
    return !in.fail() && in.peek() == std::char_traits<char>::eof();
}

std::string numberText(float v) {
    std::ostringstream out;
    out.imbue(std::locale::classic());
    for (int digits = 6; digits <= 9; ++digits) {
        out.str({});
        out << std::setprecision(digits) << v; // the default notation: printf's %g
        float back = 0;
        if (parseNumber(out.str(), back) && back == v) break;
    }
    return out.str();
}

} // namespace rf
