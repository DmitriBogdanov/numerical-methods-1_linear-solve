#pragma once

#include <chrono>



// # StaticTimer #
// Small class used for time measurements, fully static aka does not require creating local instances
struct StaticTimer {
	using Clock = std::chrono::steady_clock;
	using Milliseconds = std::chrono::milliseconds;

	inline static void start() {
		_start_timepoint = Clock::now();
	}

	inline static unsigned long long elapsed() {
		return std::chrono::duration_cast<Milliseconds>(Clock::now() - _start_timepoint).count();
			// time since last StaticTimer::start() call in ms 
	}

private:
	inline static Clock::time_point _start_timepoint = StaticTimer::Clock::now();
		// 'inline static' requires C++17
};