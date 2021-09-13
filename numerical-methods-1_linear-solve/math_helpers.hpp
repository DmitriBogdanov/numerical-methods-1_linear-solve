#pragma once

#include <algorithm>



template<typename T>
constexpr bool isZero(T value) {
	constexpr T EPSILON = static_cast<T>(1e-10);

	return std::abs(value) < EPSILON;
}

template<typename T>
constexpr T sqr(T value) {
	return value * value;
}
