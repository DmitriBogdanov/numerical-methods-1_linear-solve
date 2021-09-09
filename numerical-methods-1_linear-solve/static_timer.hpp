#pragma once

#include <chrono>



// Ќебольшой класс дл€ замера времени
struct StaticTimer {
	using Clock = std::chrono::steady_clock;
	using Milliseconds = std::chrono::milliseconds;

	inline static void start() {
		_start_timepoint = Clock::now();
	}

	inline static unsigned long long elapsed() {
		return std::chrono::duration_cast<Milliseconds>(Clock::now() - _start_timepoint).count();
			// врем€ в мс с последнего вызова StaticTimer::start() 
	}

private:
	inline static Clock::time_point _start_timepoint = StaticTimer::Clock::now();
		// 'inline static' требует C++17
};