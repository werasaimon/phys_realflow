// The registry of verification and validation cases: built once from the case files, in the order
// the board shows them (code verification first, then solution verification, then validation).
#include "Benchmark.h"

#include <algorithm>

namespace rf::verify {

static int categoryRank(const std::string& category) {
    if (category == "code-verification") return 0;
    if (category == "solution-verification") return 1;
    return 2;
}

const std::vector<Case>& registry() {
    static const std::vector<Case> cases = [] {
        std::vector<Case> all;
        addRigidCases(all);
        addGasCases(all);
        addParticleCases(all);
        addRelativityCases(all);
        addPlasmaCases(all);
        std::stable_sort(all.begin(), all.end(),
                         [](const Case& a, const Case& b) { return categoryRank(a.category) < categoryRank(b.category); });
        return all;
    }();
    return cases;
}

} // namespace rf::verify
