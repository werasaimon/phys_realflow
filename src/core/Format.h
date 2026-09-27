#pragma once
// printf into a std::string, for readouts ("%.2f м/с") of the snapshot and the scenes.
#include <cstdarg>
#include <cstdio>
#include <string>

namespace rf {

inline std::string format(const char* f, ...) {
    char buf[256];
    va_list ap;
    va_start(ap, f);
    std::vsnprintf(buf, sizeof(buf), f, ap);
    va_end(ap);
    return buf;
}

} // namespace rf
