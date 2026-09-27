#pragma once
// Minimal low-latency thread pool. MinGW's libgomp pays ~150 us per parallel region on Windows,
// which dominates iterative solvers (PCG does several regions per iteration), so the solvers use
// this pool instead: workers spin briefly between jobs and fall back to sleeping when idle.
//
// Determinism: the same input gives the same bits on 1, 2 or 64 threads - the rule of Box2D v3
// (E. Catto, "Determinism", box2d.org, 2024): a result that changes with the number of cores
// cannot be reproduced on another computer. Two rules make it so here:
//  - a reduction (parallelSum, parallelMax) cuts [0, n) into blocks that depend only on n and
//    minChunk - never on the number of threads - and adds the blocks' results in block order, so a
//    sum is rounded the same way on every machine;
//  - parallelFor / parallelRanges may cut by the number of threads, because there every index
//    writes only its own output: how the loop is cut cannot change a result.
// A solver that merges results of several workers sorts them by a stable key or writes them into
// slots addressed by index (the rigid narrow phase: one slot per body and per pair).
// tests/DeterminismTests.cpp checks the rule on the scenes with ThreadPool::setActiveThreads.

#include <algorithm>

namespace rf {

class ThreadPool {
public:
    using ChunkFn = void (*)(void* ctx, int chunk);

    static ThreadPool& instance();
    // The threads of the pool, the caller included. Per-worker scratch is sized by this number.
    int threadCount() const { return int(workers_) + 1; }
    // Which thread of the pool is calling: 0 for the thread that owns the simulation (and any
    // thread outside the pool), 1 .. threadCount() - 1 for the workers. Solvers keep one scratch
    // per worker and index it by this - a plain thread_local int, nothing with a destructor
    // (MinGW's TLS cleanup double-frees thread_local objects with destructors at thread exit).
    static int workerIndex();
    // Runs fn(ctx, c) for c in [0, chunks); the calling thread participates. Blocks until done.
    void run(int chunks, ChunkFn fn, void* ctx);

    // How many threads take part in run(): n in 1 .. threadCount(), 0 = all (the default). The
    // others stay asleep. For tests that compare 1, 2 and all threads; call it between jobs, never
    // from inside one. threadCount() does not change, so scratch sized by it stays valid.
    void setActiveThreads(int n);
    int activeThreads() const;

    ~ThreadPool();

private:
    ThreadPool();
    struct Impl;
    Impl* impl_;
    unsigned workers_ = 0;
};

// How many chunks a parallelFor over n items gets: a few per active thread, each at least
// minChunk items. Only for loops whose result does not depend on the cut (see the file head).
inline int parallelChunks(int n, int minChunk) {
    int threads = ThreadPool::instance().activeThreads();
    int byWork = std::max(1, n / std::max(1, minChunk));
    return std::max(1, std::min(threads * 4, byWork));
}

// The blocks of a reduction: at most kMaxReduceBlocks, each of at least minChunk items - a function
// of n and minChunk only, the same on every machine. 256 blocks keep 64 threads busy.
constexpr int kMaxReduceBlocks = 256;
inline int reduceBlocks(int n, int minChunk) { return std::max(1, std::min(kMaxReduceBlocks, n / std::max(1, minChunk))); }

// f(begin, end) over contiguous ranges of [0, n).
template <class F> void parallelRanges(int n, F&& f, int minChunk = 512) {
    if (n <= 0) return;
    int chunks = parallelChunks(n, minChunk);
    if (chunks == 1) { f(0, n); return; }
    struct Ctx { F* f; int n; int chunks; } ctx{&f, n, chunks};
    ThreadPool::instance().run(chunks, [](void* p, int c) {
        Ctx& x = *static_cast<Ctx*>(p);
        int b = int((long long)x.n * c / x.chunks), e = int((long long)x.n * (c + 1) / x.chunks);
        (*x.f)(b, e);
    }, &ctx);
}

// f(i) for i in [0, n).
template <class F> void parallelFor(int n, F&& f, int minChunk = 512) {
    parallelRanges(n, [&](int b, int e) { for (int i = b; i < e; ++i) f(i); }, minChunk);
}

// part[b] = f(begin, end) for block b of [0, n) cut into `blocks` equal blocks, on the pool.
template <class T, class F> void runReduceBlocks(int n, int blocks, F& f, T* part) {
    struct Ctx { F* f; int n; int blocks; T* part; } ctx{&f, n, blocks, part};
    ThreadPool::instance().run(blocks, [](void* p, int c) {
        Ctx& x = *static_cast<Ctx*>(p);
        const int b = int((long long)x.n * c / x.blocks), e = int((long long)x.n * (c + 1) / x.blocks);
        x.part[c] = (*x.f)(b, e);
    }, &ctx);
}

// Sum of the f(begin, end) partial results (T must support +=), block by block in block order.
// The partial results live on the stack: a reduction allocates nothing.
template <class T, class F> T parallelSum(int n, F&& f, int minChunk = 512) {
    if (n <= 0) return T();
    const int blocks = reduceBlocks(n, minChunk);
    if (blocks == 1) return f(0, n);
    T part[kMaxReduceBlocks];
    runReduceBlocks(n, blocks, f, part);
    T s = part[0];
    for (int c = 1; c < blocks; ++c) s += part[c];
    return s;
}

// The largest of init and the f(begin, end) partial results, on the same blocks as parallelSum
// (a maximum does not round, but a NaN would win or lose depending on the order).
template <class T, class F> T parallelMax(int n, T init, F&& f, int minChunk = 512) {
    if (n <= 0) return init;
    const int blocks = reduceBlocks(n, minChunk);
    if (blocks == 1) return std::max(init, f(0, n));
    T part[kMaxReduceBlocks];
    runReduceBlocks(n, blocks, f, part);
    T m = init;
    for (int c = 0; c < blocks; ++c) m = std::max(m, part[c]);
    return m;
}

} // namespace rf
