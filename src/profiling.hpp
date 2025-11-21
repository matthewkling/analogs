#pragma once

#include <chrono>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>

namespace analogs {

class ProfileTimer {
public:
      using Clock = std::chrono::high_resolution_clock;
      using TimePoint = Clock::time_point;
      using Duration = std::chrono::duration<double, std::milli>;

      struct Event {
            std::string name;
            double duration_ms;
            size_t count;

            Event() : name(""), duration_ms(0.0), count(0) {}
            Event(const std::string& n, double d, size_t c = 1)
                  : name(n), duration_ms(d), count(c) {}
      };

      ProfileTimer() : enabled_(false) {}

      void enable() { enabled_ = true; }
      void disable() { enabled_ = false; }
      bool is_enabled() const { return enabled_; }

      void start(const std::string& name) {
            if (!enabled_) return;
            timers_[name] = Clock::now();
      }

      void stop(const std::string& name) {
            if (!enabled_) return;

            auto end = Clock::now();
            auto it = timers_.find(name);
            if (it == timers_.end()) return;

            double ms = Duration(end - it->second).count();

            auto& event = events_[name];
            event.name = name;
            event.duration_ms += ms;
            event.count += 1;

            timers_.erase(it);
      }

      void record(const std::string& name, double value_ms) {
            if (!enabled_) return;

            auto& event = events_[name];
            event.name = name;
            event.duration_ms += value_ms;
            event.count += 1;
      }

      void increment_counter(const std::string& name, size_t delta = 1) {
            if (!enabled_) return;
            counters_[name] += delta;
      }

      std::map<std::string, Event> get_events() const {
            return events_;
      }

      std::map<std::string, size_t> get_counters() const {
            return counters_;
      }

      void clear() {
            events_.clear();
            counters_.clear();
            timers_.clear();
      }

      // RAII helper for scoped timing
      class ScopedTimer {
      public:
            ScopedTimer(ProfileTimer& pt, const std::string& name)
                  : timer_(pt), name_(name) {
                  timer_.start(name_);
            }
            ~ScopedTimer() {
                  timer_.stop(name_);
            }
      private:
            ProfileTimer& timer_;
            std::string name_;
      };

private:
      bool enabled_;
      std::map<std::string, TimePoint> timers_;
      std::map<std::string, Event> events_;
      std::map<std::string, size_t> counters_;
};

// Thread-local profiler instance
extern thread_local ProfileTimer g_profiler;

// Convenience macros
#define PROFILE_START(name) if (analogs::g_profiler.is_enabled()) analogs::g_profiler.start(name)
#define PROFILE_STOP(name) if (analogs::g_profiler.is_enabled()) analogs::g_profiler.stop(name)
#define PROFILE_SCOPE(name) analogs::ProfileTimer::ScopedTimer _prof_scope_##__LINE__(analogs::g_profiler, name)
#define PROFILE_COUNT(name, delta) if (analogs::g_profiler.is_enabled()) analogs::g_profiler.increment_counter(name, delta)

// Per-worker profiling aggregator
class WorkerProfileAggregator {
public:
      void add_worker_profile(const ProfileTimer& worker_prof) {
            auto events = worker_prof.get_events();
            for (const auto& kv : events) {
                  auto& agg = aggregated_[kv.first];
                  agg.name = kv.first;
                  agg.duration_ms += kv.second.duration_ms;
                  agg.count += kv.second.count;
            }

            auto counters = worker_prof.get_counters();
            for (const auto& kv : counters) {
                  counters_[kv.first] += kv.second;
            }
      }

      std::map<std::string, ProfileTimer::Event> get_aggregated() const {
            return aggregated_;
      }

      std::map<std::string, size_t> get_counters() const {
            return counters_;
      }

private:
      std::map<std::string, ProfileTimer::Event> aggregated_;
      std::map<std::string, size_t> counters_;
};

} // namespace analogs
