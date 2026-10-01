#include "tests/TestRunner.h"
#include "tests/Tests.h"
int main() {
    run("motor work analytical and observation", testMotorWork);
    run("motor work rollback", testMotorWorkRollback);
    run("nut bolt assembly", testNutBoltAssembly);
    run("nut bolt cycle", testNutBoltCycle);
    run("chain drive cycle", testChainDriveCycle);
    if (const char* path=std::getenv("RF_JUNIT")) writeJUnit(path);
    return g_failures ? 1 : 0;
}
