#pragma once

// ---------------------------------------------------------------------------
// Lightweight, compile-gated phase profiler for the parallel workers.
//
// Enable by building with -DANALOGS_PROFILE (e.g.
//   PKG_CPPFLAGS="-DANALOGS_PROFILE" R CMD INSTALL .
// or set it in ~/.R/Makevars temporarily). When the flag is absent, every
// macro below expands to nothing and there is ZERO runtime cost: no timers,
// no registration, no clock calls.
//
// Design notes:
//   - Each worker thread accumulates nanosecond totals and call counts into a
//     thread_local ProfileTimer. This keeps the hot path lock-free.
//   - On first use, each thread's timer registers a pointer to itself into a
//     shared registry (guarded by a mutex). Registration happens once per
//     thread, not per focal, so it stays off the hot path.
//   - report() locks the registry and sums across all threads, printing a
//     compact table. This is what you paste back to drive optimization
//     sequencing.
//
// Phases (per focal):
//   GATHER - lattice candidate collection (query / knn_query)
//   EXACT  - exact geo/environment distance pass + selection over candidates
//   AGG    - aggregation / finalize / solve_ridge (AggWorker only)
// ---------------------------------------------------------------------------

#ifdef ANALOGS_PROFILE

#include <cstdint>
#include <vector>
#include <mutex>
#include <chrono>
#include <Rcpp.h>

namespace analogs {

enum class ProfilePhase : int {
      GATHER = 0,
            EXACT  = 1,
            AGG    = 2,
            N_PHASES = 3
};

struct ProfileTimer {
      uint64_t ns[3];      // accumulated nanoseconds per phase
      uint64_t calls[3];   // number of scoped intervals per phase

      ProfileTimer() {
            for (int p = 0; p < 3; ++p) { ns[p] = 0; calls[p] = 0; }
            registry_register(this);
      }

      inline void add(int phase, uint64_t dt) {
            ns[phase]    += dt;
            calls[phase] += 1;
      }

      // --- shared registry across threads ------------------------------------
      static std::mutex& registry_mutex() {
            static std::mutex m;
            return m;
      }
      static std::vector<ProfileTimer*>& registry() {
            static std::vector<ProfileTimer*> v;
            return v;
      }
      static void registry_register(ProfileTimer* t) {
            std::lock_guard<std::mutex> lk(registry_mutex());
            registry().push_back(t);
      }

      // Zero all registered timers (call before a parallelFor).
      static void reset_all() {
            std::lock_guard<std::mutex> lk(registry_mutex());
            for (ProfileTimer* t : registry()) {
                  for (int p = 0; p < 3; ++p) { t->ns[p] = 0; t->calls[p] = 0; }
            }
      }

      // Sum across threads and print a compact table (call after a parallelFor).
      static void report(const char* label) {
            std::lock_guard<std::mutex> lk(registry_mutex());

            uint64_t tot_ns[3]    = {0, 0, 0};
            uint64_t tot_calls[3] = {0, 0, 0};
            int n_threads_active = 0;

            for (ProfileTimer* t : registry()) {
                  bool active = false;
                  for (int p = 0; p < 3; ++p) {
                        tot_ns[p]    += t->ns[p];
                        tot_calls[p] += t->calls[p];
                        if (t->calls[p] > 0) active = true;
                  }
                  if (active) ++n_threads_active;
            }

            const uint64_t grand_ns = tot_ns[0] + tot_ns[1] + tot_ns[2];
            const double grand_ms   = grand_ns / 1.0e6;

            const char* names[3] = {"GATHER", "EXACT ", "AGG   "};

            Rcpp::Rcout << "\n[analogs profile] " << label
                        << "  (threads active: " << n_threads_active << ")\n";
            Rcpp::Rcout << "  phase    total_ms      %tracked   calls        mean_ns/call\n";
            for (int p = 0; p < 3; ++p) {
                  const double ms  = tot_ns[p] / 1.0e6;
                  const double pct = grand_ns > 0
                  ? (100.0 * static_cast<double>(tot_ns[p]) / static_cast<double>(grand_ns))
                        : 0.0;
                  const double mean_ns = tot_calls[p] > 0
                  ? (static_cast<double>(tot_ns[p]) / static_cast<double>(tot_calls[p]))
                        : 0.0;
                  Rcpp::Rcout << "  " << names[p]
                              << "  " << ms
                              << "        " << pct << "%"
                              << "     " << tot_calls[p]
                              << "        " << mean_ns << "\n";
            }
            Rcpp::Rcout << "  total tracked: " << grand_ms << " ms\n\n";
      }
};

// Thread-local profiler accessor. Using a function-local static (Meyers
// singleton) rather than a namespace-scope `extern thread_local` avoids the
// cross-translation-unit TLS-wrapper symbol that fails to resolve under
// macOS/Clang flat-namespace dlopen (as R uses). No separate definition in
// any .cpp is required; each TU gets the same per-thread instance via this
// inline accessor. The ProfileTimer constructor registers itself into the
// shared registry on first use per thread.
inline ProfileTimer& profiler() {
      static thread_local ProfileTimer tls_profiler;
      return tls_profiler;
}

// RAII scope guard: times from construction to destruction, adds to phase.
struct ProfileScope {
      int phase;
      std::chrono::steady_clock::time_point t0;
      explicit ProfileScope(int phase_)
            : phase(phase_), t0(std::chrono::steady_clock::now()) {}
      ~ProfileScope() {
            const auto t1 = std::chrono::steady_clock::now();
            const uint64_t dt = static_cast<uint64_t>(
                  std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count());
            profiler().add(phase, dt);
      }
};

} // namespace analogs

// Concatenation helpers so we can generate unique scope-guard variable names.
#define ANALOGS_PROF_CAT_(a, b) a##b
#define ANALOGS_PROF_CAT(a, b) ANALOGS_PROF_CAT_(a, b)

// Bracket a lexical scope. Usage: { ANALOGS_PROFILE_SCOPE(GATHER); ...work... }
#define ANALOGS_PROFILE_SCOPE(PHASE)                                \
::analogs::ProfileScope ANALOGS_PROF_CAT(_analogs_prof_, __LINE__)( \
            static_cast<int>(::analogs::ProfilePhase::PHASE))

#define ANALOGS_PROFILE_RESET() ::analogs::ProfileTimer::reset_all()
#define ANALOGS_PROFILE_REPORT(LABEL) ::analogs::ProfileTimer::report(LABEL)

#else  // !ANALOGS_PROFILE  -- all macros compile to nothing

#define ANALOGS_PROFILE_SCOPE(PHASE) ((void)0)
#define ANALOGS_PROFILE_RESET() ((void)0)
#define ANALOGS_PROFILE_REPORT(LABEL) ((void)0)

#endif // ANALOGS_PROFILE
