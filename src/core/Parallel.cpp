#include "core/Parallel.h"

#include <cstdlib>

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <mutex>
#include <thread>

#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__)
#include <immintrin.h>
#define RF_PAUSE() _mm_pause()
#else
#define RF_PAUSE() std::this_thread::yield()
#endif

namespace rf {

static thread_local bool tl_inPool = false;

// A job lives in one of two slots selected by generation parity. The shared ticket packs
// (generation << 32 | next chunk index); workers claim chunks with CAS, so a worker holding a
// stale copy of a job can never claim a chunk of a newer one.
struct ThreadPool::Impl {
    struct Slot {
        std::atomic<ChunkFn> fn{nullptr};
        std::atomic<void*> ctx{nullptr};
        std::atomic<uint32_t> chunks{0};
        std::atomic<int> remaining{0};
    };
    Slot slots[2];
    std::atomic<uint64_t> ticket{0};
    std::atomic<bool> quit{false};
    std::atomic<int> sleeping{0};
    std::mutex sleepMutex;
    std::condition_variable cv;
    std::mutex runMutex;
    std::vector<std::thread> threads;

    // Tries to execute one chunk of the current job. Returns false when there is nothing to do.
    bool tryWork() {
        uint64_t t = ticket.load(std::memory_order_acquire);
        uint32_t gen = uint32_t(t >> 32), idx = uint32_t(t);
        Slot& s = slots[gen & 1];
        ChunkFn fn = s.fn.load(std::memory_order_relaxed);
        void* ctx = s.ctx.load(std::memory_order_relaxed);
        uint32_t chunks = s.chunks.load(std::memory_order_relaxed);
        if (idx >= chunks) return false;
        if (!ticket.compare_exchange_weak(t, t + 1, std::memory_order_acq_rel)) return true; // retry
        tl_inPool = true;
        fn(ctx, int(idx));
        tl_inPool = false;
        s.remaining.fetch_sub(1, std::memory_order_acq_rel);
        return true;
    }

    void workerLoop() {
        using clock = std::chrono::steady_clock;
        while (!quit.load(std::memory_order_relaxed)) {
            if (tryWork()) continue;
            // Idle: spin, then yield, then sleep.
            auto idleStart = clock::now();
            int spins = 0;
            while (!quit.load(std::memory_order_relaxed)) {
                uint64_t t = ticket.load(std::memory_order_acquire);
                if (uint32_t(t) < slots[uint32_t(t >> 32) & 1].chunks) break;
                if (++spins < 4000) { RF_PAUSE(); continue; }
                if (clock::now() - idleStart < std::chrono::milliseconds(3)) { std::this_thread::yield(); continue; }
                std::unique_lock<std::mutex> lk(sleepMutex);
                sleeping.fetch_add(1);
                cv.wait_for(lk, std::chrono::milliseconds(20), [&] {
                    uint64_t t2 = ticket.load(std::memory_order_acquire);
                    return quit.load() || uint32_t(t2) < slots[uint32_t(t2 >> 32) & 1].chunks;
                });
                sleeping.fetch_sub(1);
                idleStart = clock::now();
                spins = 0;
            }
        }
    }
};

ThreadPool& ThreadPool::instance() {
    static ThreadPool pool;
    return pool;
}

ThreadPool::ThreadPool() : impl_(new Impl) {
    unsigned hw = std::thread::hardware_concurrency();
    // RF_THREADS=n limits the pool (profiling, benchmarks); default: all hardware threads.
    if (const char* env = std::getenv("RF_THREADS"); env && std::atoi(env) > 0) hw = unsigned(std::atoi(env));
    workers_ = hw > 1 ? hw - 1 : 0;
    for (unsigned i = 0; i < workers_; ++i) impl_->threads.emplace_back([this] { impl_->workerLoop(); });
}

ThreadPool::~ThreadPool() {
    impl_->quit = true;
    impl_->cv.notify_all();
    for (auto& t : impl_->threads) t.join();
    delete impl_;
}

void ThreadPool::run(int chunks, ChunkFn fn, void* ctx) {
    if (chunks <= 0) return;
    // Nested parallel calls from inside a job run serially (avoids deadlock on runMutex).
    if (workers_ == 0 || chunks == 1 || tl_inPool) {
        for (int c = 0; c < chunks; ++c) fn(ctx, c);
        return;
    }
    std::lock_guard<std::mutex> guard(impl_->runMutex);
    Impl& I = *impl_;
    uint64_t t = I.ticket.load(std::memory_order_acquire);
    uint32_t gen = uint32_t(t >> 32) + 1;
    Impl::Slot& s = I.slots[gen & 1];
    s.fn.store(fn, std::memory_order_relaxed);
    s.ctx.store(ctx, std::memory_order_relaxed);
    s.chunks.store(uint32_t(chunks), std::memory_order_relaxed);
    s.remaining.store(chunks, std::memory_order_relaxed);
    I.ticket.store(uint64_t(gen) << 32, std::memory_order_release);
    if (I.sleeping.load() > 0) {
        std::lock_guard<std::mutex> lk(I.sleepMutex);
        I.cv.notify_all();
    }
    while (I.tryWork()) {
    }
    while (s.remaining.load(std::memory_order_acquire) > 0) RF_PAUSE();
}

} // namespace rf
