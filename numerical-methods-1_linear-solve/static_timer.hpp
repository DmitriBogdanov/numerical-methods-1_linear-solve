#pragma once

#include <chrono>



// ��������� ����� ��� ������ �������
struct StaticTimer {
	using Clock = std::chrono::steady_clock;
	using Milliseconds = std::chrono::milliseconds;

	inline static void start() {
		_start_timepoint = Clock::now();
	}

	inline static unsigned long long elapsed() {
		return std::chrono::duration_cast<Milliseconds>(Clock::now() - _start_timepoint).count();
			// ����� � �� � ���������� ������ StaticTimer::start() 
	}

private:
	inline static Clock::time_point _start_timepoint = StaticTimer::Clock::now();
		// 'inline static' ������� C++17
};