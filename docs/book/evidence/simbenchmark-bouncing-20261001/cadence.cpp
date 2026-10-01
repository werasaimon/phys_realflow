// Diagnostic only: isolate integration/restitution error using the shared benchmark scene.
#define main bouncingBenchmarkEntry
#include "../../../../tools/SimBenchmarkBouncing.cpp"
#undef main
int main() {
    for (const std::string mode : {"fixed-discrete", "refined-discrete", "atomic-ccd"}) {
        BouncingSdk sdk(20, mode == "atomic-ccd");
        Work work;
        Metrics metrics;
        const double initial = observe(sdk.states(), 0, metrics, false);
        double impactGain = 0, freeFlightChange = 0, predictedFreeFlight = 0;
        int impacts = 0;
        for (int frame = 0; frame < 2000; ++frame) {
            const int parts = mode == "refined-discrete" && sdk.states()[0].position.y < 0.2f ? 16 : 1;
            const float h = 0.01f / float(parts);
            for (int part = 0; part < parts; ++part) {
                const auto before = sdk.states()[0];
                const double oldEnergy = observe(sdk.states(), work.time, metrics, false);
                if (!sdk.step(h, work)) return 2;
                const auto after = sdk.states()[0];
                const double energy = observe(sdk.states(), work.time, metrics, false);
                if (before.velocity.y < 0 && after.velocity.y > 0) {
                    impactGain += energy - oldEnergy;
                    ++impacts;
                } else {
                    freeFlightChange += energy - oldEnergy;
                    predictedFreeFlight -= 0.5 * count * mass * 9.81 * 9.81 * double(h) * h;
                }
            }
        }
        const double final = observe(sdk.states(), work.time, metrics, false);
        std::cout << "{\"mode\":\"" << mode << "\",\"initial_J\":" << number(initial)
                  << ",\"final_J\":" << number(final) << ",\"impacts\":" << impacts
                  << ",\"impact_gain_J\":" << number(impactGain)
                  << ",\"free_flight_change_J\":" << number(freeFlightChange)
                  << ",\"predicted_euler_change_J\":" << number(predictedFreeFlight)
                  << ",\"accepted_time_s\":" << number(work.time) << ",\"accepted\":" << work.accepted << ",\"rejected\":" << work.rejected << "}\n";
    }
}
