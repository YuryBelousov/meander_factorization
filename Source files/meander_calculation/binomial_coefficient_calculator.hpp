#pragma once
#include <algorithm>
#include <cstddef>
#include <unordered_map>
#include <utility>

/*
This file contains an auxiliary technical class for calculating binomial coefficients
(as well as storing those already calculated).

Template parameter:
	- coef_type: the numeric type used for internal calculations (e.g., long long). It should support
	arithmetic operations (+, *, assignment) and construction from size_t.
*/


template <typename coef_type>
class binomial_coefficient_calculator {
private:
	struct index_hash {
		size_t operator()(const std::pair<size_t, size_t>& index) const noexcept {
			return index.first ^ (index.second + size_t(0x9e3779b9) + (index.first << 6) + (index.first >> 2));
		}
	};
	std::unordered_map<std::pair<size_t, size_t>, coef_type, index_hash> binomial_coef;

public:
	binomial_coefficient_calculator() = default;

	const coef_type& get_binomial_coef(size_t n, size_t k) {
		if (k > n) {
			static const coef_type zero(0);
			return zero;
		}
		k = std::min(k, n - k);
		auto index = std::make_pair(n, k);
		auto [it, inserted] = binomial_coef.try_emplace(index, coef_type{});
		auto& value = it->second; // References remain valid when recursive calls rehash the cache.
		if (inserted) {
			try {
				switch (k) {
				case 0:	value = coef_type(1);	break;
				case 1:	value = coef_type(n);	break;
				case 2:	value = coef_type(n / 2) * coef_type(n % 2 ? n : n - 1);	break;
				default:
					value = get_binomial_coef(n - 1, k) + get_binomial_coef(n - 1, k - 1);	break;
				}
			}
			catch (...) {
				binomial_coef.erase(index);
				throw;
			}
		}
		return value;
	}
};