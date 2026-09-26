#pragma once
// Runs the Simulation on a background thread. The GUI never touches the Simulation directly:
// it posts commands (executed between frames on the simulation thread) and reads immutable
// RenderSnapshots published after every frame.

#include "scene/Simulation.h"

#include <atomic>
#include <condition_variable>
#include <deque>
#include <functional>
#include <memory>
#include <mutex>
#include <thread>

class SimController {
public:
    using Command = std::function<void(rf::Simulation&)>;

    SimController();
    ~SimController();

    void post(Command cmd);
    void setRunning(bool on);
    bool isRunning() const { return running_; }
    void requestStep();
    void setRealtimeLimit(bool on) { realtime_ = on; }

    std::shared_ptr<const rf::RenderSnapshot> snapshot() const;
    uint64_t serial() const { return serial_; }

private:
    void loop();
    void publish();

    std::unique_ptr<rf::Simulation> sim_;
    std::thread thread_;
    std::mutex mtx_;
    std::condition_variable cv_;
    std::deque<Command> queue_;
    int stepRequests_ = 0;
    std::atomic<bool> running_{false}, quit_{false}, realtime_{true};

    mutable std::mutex snapMtx_;
    std::shared_ptr<const rf::RenderSnapshot> snap_;
    std::atomic<uint64_t> serial_{0};
};
