#pragma once
#include <atomic>
#include <limits>

#include "enumerate_irreducible_meanders.hpp"
#include "submeander.h"

// Compatibility with the original example, with the corrected [n][k] indexing.
inline std::atomic<unsigned long long> irreducible_meander_numbers[max_n + 1][max_n + 1]{};

inline void brute_force() {
	const auto numbers = enumerate_irreducible_meanders(max_n);
	for (int n = 0; n <= max_n; ++n)
		for (int k = 0; n + 2 * k <= max_n; ++k)
			if (numbers[n][k] > std::numeric_limits<unsigned long long>::max())
				throw std::overflow_error("Use enumerate_irreducible_meanders for counts exceeding 64 bits");
	for (int n = 0; n <= max_n; ++n)
		for (int k = 0; k <= max_n; ++k)
			irreducible_meander_numbers[n][k] = n + 2 * k <= max_n
				? numbers[n][k].convert_to<unsigned long long>() : 0;
}
