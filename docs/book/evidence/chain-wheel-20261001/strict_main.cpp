#include "tests/TestRunner.h"
#include "tests/Tests.h"
int main() {
    run("chain wheel assembly (strict FP)",testChainWheelAssembly);
    run("chain wheel loaded, pulled and released (strict FP)",testChainWheelLoaded);
    if(const char* path=std::getenv("RF_JUNIT")) writeJUnit(path);
    return g_failures;
}
