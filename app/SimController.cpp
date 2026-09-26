#include "SimController.h"

#include <chrono>

SimController::SimController() {
    sim_ = std::make_unique<rf::Simulation>();
    publish();
    thread_ = std::thread([this] { loop(); });
}

SimController::~SimController() {
    quit_ = true;
    cv_.notify_all();
    if (thread_.joinable()) thread_.join();
}

void SimController::post(Command cmd) {
    {
        std::lock_guard<std::mutex> lk(mtx_);
        queue_.push_back(std::move(cmd));
    }
    cv_.notify_all();
}

void SimController::setRunning(bool on) {
    running_ = on;
    cv_.notify_all();
}

void SimController::requestStep() {
    {
        std::lock_guard<std::mutex> lk(mtx_);
        ++stepRequests_;
    }
    cv_.notify_all();
}

std::shared_ptr<const rf::RenderSnapshot> SimController::snapshot() const {
    std::lock_guard<std::mutex> lk(snapMtx_);
    return snap_;
}

void SimController::publish() {
    auto s = std::make_shared<rf::RenderSnapshot>();
    sim_->fillSnapshot(*s);
    {
        std::lock_guard<std::mutex> lk(snapMtx_);
        snap_ = std::move(s);
    }
    ++serial_;
}

void SimController::loop() {
    using clock = std::chrono::steady_clock;
    while (!quit_) {
        std::deque<Command> cmds;
        bool doStep = false;
        {
            std::unique_lock<std::mutex> lk(mtx_);
            cv_.wait_for(lk, std::chrono::milliseconds(50),
                         [&] { return quit_ || !queue_.empty() || running_ || stepRequests_ > 0; });
            if (quit_) break;
            cmds.swap(queue_);
            if (running_) doStep = true;
            else if (stepRequests_ > 0) { --stepRequests_; doStep = true; }
        }
        for (auto& c : cmds) c(*sim_);
        if (!doStep && cmds.empty()) continue;

        auto t0 = clock::now();
        if (doStep) sim_->stepFrame();
        publish();

        // Optional real-time cap: do not run the SPH / rigid scenes faster than wall-clock.
        if (doStep && realtime_ && sim_->mode() != rf::SimMode::WindTunnel) {
            auto target = t0 + std::chrono::microseconds(int(sim_->frameDt * 1e6f));
            std::this_thread::sleep_until(target);
        }
    }
}
