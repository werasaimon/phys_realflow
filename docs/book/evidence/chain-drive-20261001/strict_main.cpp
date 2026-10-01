#include "tests/TestRunner.h"
#include "tests/Tests.h"
int main() {
    run("chain drive (strict FP)",testChainDriveCycle);
    if(const char* path=std::getenv("RF_JUNIT")) writeJUnit(path);
    return g_failures;
}
