#pragma once
#include <unordered_map>
#include <vector>

/*
This file contains an auxiliary technical class for calculating binomial coefficients
(as well as storing those already calculated).

Template parameter:
	- number_type: the numeric type used for internal calculations (e.g., long long). It should support
	basic arithmetic operations (+, -, *, assignment) and allow implicit conversion from int.
*/


template <typename coef_type>
class binomial_coefficient_calculator {
private:
	std::unordered_map<uint32_t, coef_type> binomial_coef;

private:
	static uint32_t pack(size_t n, size_t k) noexcept {
		constexpr int bit_shift = 16;
		return (uint32_t(n) << bit_shift) | uint32_t(k);
	};

public:
	binomial_coefficient_calculator() = default;

	const coef_type& get_binomial_coef(size_t n, size_t k) noexcept {
		k = std::min(k, n - k);
		auto index = pack(n, k);
		auto [it, inserted] = binomial_coef.try_emplace(index, coef_type{});

		if (inserted)
			switch (k) {
			case 0:	it->second = coef_type(1);	break;
			case 1:	it->second = coef_type(n);	break;
			case 2:	it->second = coef_type(n) * coef_type(n - 1) / coef_type(2);	break;
			default:
				it->second = get_binomial_coef(n - 1, k) + get_binomial_coef(n - 1, k - 1);	break;
			}
		return it->second;
	}
};