#pragma once
// Minimal low-latency thread pool. MinGW's libgomp pays ~150 us per parallel region on Windows,
// which dominates iterative solvers (PCG does several regions per iteration), so the solvers use
// this pool instead: workers spin briefly between jobs and fall back to sleeping when idle.

#include <algorithm>
#include <vector>

namespace rf {

class ThreadPool {
public:
    using ChunkFn = void (*)(void* ctx, int chunk);

    static ThreadPool& instance();
    int threadCount() const { return int(workers_) + 1; }
    // Runs fn(ctx, c) for c in [0, chunks); the calling thread participates. Blocks until done.
    void run(int chunks, ChunkFn fn, void* ctx);

    ~ThreadPool();

private:
    ThreadPool();
    struct Impl;
    Impl* impl_;
    unsigned workers_ = 0;
};

inline int parallelChunks(int n, int minChunk) {
    int threads = ThreadPool::instance().threadCount();
    int byWork = std::max(1, n / std::max(1, minChunk));
    return std::max(1, std::min(threads * 4, byWork));
}

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

// Sum of f(begin, end) partial results (T must support +=).
template <class T, class F> T parallelSum(int n, F&& f, int minChunk = 512) {
    if (n <= 0) return T();
    int chunks = parallelChunks(n, minChunk);
    std::vector<T> part(chunks, T());
    struct Ctx { F* f; int n; int chunks; T* part; } ctx{&f, n, chunks, part.data()};
    if (chunks == 1) return f(0, n);
    ThreadPool::instance().run(chunks, [](void* p, int c) {
        Ctx& x = *static_cast<Ctx*>(p);
        int b = int((long long)x.n * c / x.chunks), e = int((long long)x.n * (c + 1) / x.chunks);
        x.part[c] = (*x.f)(b, e);
    }, &ctx);
    T s = part[0];
    for (int c = 1; c < chunks; ++c) s += part[c];
    return s;
}

template <class T, class F> T parallelMax(int n, T init, F&& f, int minChunk = 512) {
    if (n <= 0) return init;
    int chunks = parallelChunks(n, minChunk);
    std::vector<T> part(chunks, init);
    struct Ctx { F* f; int n; int chunks; T* part; } ctx{&f, n, chunks, part.data()};
    if (chunks == 1) return std::max(init, f(0, n));
    ThreadPool::instance().run(chunks, [](void* p, int c) {
        Ctx& x = *static_cast<Ctx*>(p);
        int b = int((long long)x.n * c / x.chunks), e = int((long long)x.n * (c + 1) / x.chunks);
        x.part[c] = (*x.f)(b, e);
    }, &ctx);
    T m = init;
    for (const T& v : part) m = std::max(m, v);
    return m;
}

} // namespace rf
