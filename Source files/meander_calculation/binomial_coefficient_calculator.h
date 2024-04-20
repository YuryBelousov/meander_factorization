#pragma once
#include <map>

#include "boost/multiprecision/cpp_int.hpp"
#include "meander_calculation/integer_encoder.h"

/*
This file contains an auxiliary technical class for calculating binomial coefficients 
(as well as storing those already calculated).
*/

using big_int = boost::multiprecision::uint1024_t;

class binomial_coefficient_calculator {
private:
	std::map<int_pair_encoded, big_int> binomial_coef;

public:
	binomial_coefficient_calculator() = default;

	const big_int& get_binomial_coef(int n, int k) noexcept {
		k = std::min(k, n - k);
		auto index = int_encoder::encode(n, k);
		if (binomial_coef.contains(index))
			return binomial_coef[index];
		big_int new_value = 1;
		switch (k) {
		case 0:
			break;
		case 1:
			new_value = n;
			break;
		case 2:
			new_value = n * (n - 1) / 2;
			break;
		default:
			new_value = get_binomial_coef(n - 1, k) + get_binomial_coef(n - 1, k - 1);
			break;
		}
		binomial_coef[index] = new_value;
		return binomial_coef[index];
	}
};