// An independent caller exercises failed-step recovery with all four broad-phase strategies.
#include "tests/TestRunner.h"
using namespace rf;
static void configure(RigidWorld& w, int kind) {
    if (kind == 0) w.setBroadPhase(std::make_unique<BruteForceBroadPhase>());
    if (kind == 1) w.setBroadPhase(std::make_unique<BvhBroadPhase>());
    if (kind == 2) w.setBroadPhase(std::make_unique<AABBTreeBroadPhase>(0.1f));
    if (kind == 3) w.setBroadPhase(std::make_unique<SweepAndPruneBroadPhase>(0.1f));
    w.params.sleeping = w.params.collideWithDomain = false;
    w.params.linearDamping = w.params.angularDamping = 0;
    w.addBox(Vector3(0), {0.01f,1,1}, Quaternion(),0,Vector3(1));
    w.addSphere({-1.5f,0,0},0.02f,1000,Vector3(1));
    w.addBox({5,-0.1f,0},{1,0.1f,1},Quaternion(),0,Vector3(1));
    w.addBox({5,0.1f,0},Vector3(0.1f),Quaternion(),1000,Vector3(1));
    for(int i=0;i<10;++i) w.step(0.001f);
    w.bodies()[1].vel={100,0,0};
    w.bodies()[3].force={2,3,4};
}
static uint64_t hash(const RigidWorld& w) {
    StateHash h;
    for(auto& b:w.bodies()) {h.add(b.pos);h.add(b.rot);h.add(b.vel);h.add(b.angVel);h.add(b.force);h.add(b.torque);}
    return h.h;
}
int main() {
    for(int kind=0;kind<4;++kind) {
        RigidWorld w,ref;configure(w,kind);configure(ref,kind);
        const auto before=hash(w);w.params.ccdMaxSubdivisions=1;
        CHECK(!w.tryStep(0.02f), "budget fixture must fail");
        CHECK(hash(w)==before, "failed step changed state");
        w.params.ccdMaxSubdivisions=16;
        for(int i=0;i<60;++i) {
            w.step(0.001f);ref.step(0.001f);
            CHECK(hash(w)==hash(ref), "broad phase %d mismatch after retry at %d",kind,i);
        }
        std::printf("%s: 60 subsequent states compared\n",w.broadPhase().name());
    }
    return g_failures;
}
