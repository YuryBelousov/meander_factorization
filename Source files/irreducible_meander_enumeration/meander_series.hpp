#pragma once
#include <algorithm>
#include <array>
#include <stdexcept>
#include <vector>

#include "boost/multiprecision/cpp_int.hpp"

namespace irreducible_meander_enumeration {

// Compositional inversion needs signed integers, possibly larger than the final counts.
using integer = boost::multiprecision::cpp_int;

class meander_series {
public:
	struct term {
		int n, k;
		integer value;
	};

private:
	int weight_limit, transverse_limit, width;
	std::vector<integer> coefficients;

public:
	// Retain x^n t^k with n + 2k <= max_n and n <= max_transverse.
	explicit meander_series(int max_n, int max_transverse = -1)
		: weight_limit(max_n), transverse_limit(max_transverse < 0 ? max_n : std::min(max_n, max_transverse)),
		width(max_n / 2 + 1) {
		if (max_n < 0) throw std::invalid_argument("negative series cutoff");
		coefficients.resize((transverse_limit + 1) * width);
	}

	int limit() const { return weight_limit; }
	int max_transverse() const { return transverse_limit; }
	bool contains(int n, int k) const {
		return n >= 0 && n <= transverse_limit && k >= 0 && n + 2 * k <= weight_limit;
	}
	const integer& operator()(int n, int k) const {
		static const integer zero = 0;
		return contains(n, k) ? coefficients[n * width + k] : zero;
	}
	integer& operator()(int n, int k) {
		if (!contains(n, k)) throw std::out_of_range("series coefficient outside cutoff");
		return coefficients[n * width + k];
	}
	std::vector<term> terms() const {
		std::vector<term> result;
		for (int n = 0; n <= transverse_limit; ++n)
			for (int k = 0; n + 2 * k <= weight_limit; ++k)
				if ((*this)(n, k) != 0) result.push_back({n, k, (*this)(n, k)});
		return result;
	}
	std::vector<std::array<int, 2>> orders_by_degree() const {
		std::vector<std::array<int, 2>> result;
		for (int degree = 0; degree <= weight_limit; ++degree)
			for (int k = 0; k <= degree && degree + k <= weight_limit; ++k)
				if (degree - k <= transverse_limit) result.push_back({degree - k, k});
		return result;
	}
};

// The cubic recurrence excludes the two snakes of total order one.
// Maintain S^2 as S is calculated, avoiding repeated cubic convolutions.
inline meander_series iterated_snakes(int max_n, int max_transverse = -1) {
	meander_series snakes(max_n, max_transverse), square(max_n, max_transverse);
	std::vector<meander_series::term> established;
	constexpr int quadratic[][3] = {
		{0, 0, 2}, {0, 1, -8}, {0, 2, 7}, {1, 0, -5}, {1, 1, 5}
	};
	constexpr int linear[][3] = {
		{0, 1, -3}, {0, 2, -5}, {0, 3, 11}, {1, 0, 1},
		{1, 1, -16}, {1, 2, 18}, {2, 0, -7}, {2, 1, 7}
	};
	constexpr int constant[][3] = {
		{0, 2, -1}, {0, 3, 2}, {0, 4, 5}, {1, 1, -4}, {1, 2, -3},
		{1, 3, 13}, {2, 0, -1}, {2, 1, -8}, {2, 2, 11}, {3, 0, -3}, {3, 1, 3}
	};
	const auto& known_snakes = snakes;
	const auto& known_square = square;
	auto cubic_term = [&](int n, int k) {
		integer value = 0;
		for (const auto& term : established)
			value += term.value * known_square(n - term.n, k - term.k);
		return value;
	};
	for (const auto [n, k] : snakes.orders_by_degree()) {
		if (n + k < 2) continue;
		integer value = cubic_term(n, k) - cubic_term(n, k - 1);
		for (const auto& [i, j, coefficient] : quadratic)
			value -= coefficient * known_square(n - i, k - j);
		for (const auto& [i, j, coefficient] : linear)
			value -= coefficient * known_snakes(n - i, k - j);
		for (const auto& [i, j, coefficient] : constant)
			if (n == i && k == j) value -= coefficient;
		if (value < 0) throw std::logic_error("negative snake coefficient");
		snakes(n, k) = value;
		if (value == 0) continue;
		for (const auto& term : established)
			if (square.contains(n + term.n, k + term.k))
				square(n + term.n, k + term.k) += 2 * value * term.value;
		if (square.contains(2 * n, 2 * k)) square(2 * n, 2 * k) += value * value;
		established.push_back({n, k, value});
	}
	return snakes;
}

// Find F with F(P, Q) = x + t, where P = x + ... is odd and Q = t + ... is even in x.
// Eliminate in degree n + k: substitution is unit triangular in this order.
// Weighted order n + 2k would fail, since Q can contain x^2 at the same weight as t.
inline meander_series inverse_sum(const meander_series& P, const meander_series& Q) {
	if (P.limit() != Q.limit() || P.max_transverse() != Q.max_transverse())
		throw std::invalid_argument("series cutoffs differ");
	if (P.limit() < 2 || P(1, 0) != 1 || Q(0, 1) != 1)
		throw std::invalid_argument("composition map has the wrong linear part");
	const auto P_terms = P.terms(), Q_terms = Q.terms();
	for (const auto& term : P_terms)
		if (term.n % 2 != 1) throw std::invalid_argument("P is not odd in x");
	for (const auto& term : Q_terms)
		if (term.n % 2 != 0 || term.n + 2 * term.k < 2)
			throw std::invalid_argument("Q has invalid parity or valuation");

	const int max_n = P.limit(), max_transverse = P.max_transverse(), width = max_n / 2 + 1;
	meander_series result(max_n, max_transverse), residual(max_n, max_transverse);
	residual(1, 0) = residual(0, 1) = 1;
	std::vector<meander_series> images((max_transverse + 1) * width, result);
	auto image = [&](int n, int k) -> meander_series& { return images[n * width + k]; };
	image(0, 0)(0, 0) = 1;
	for (const auto [n, k] : result.orders_by_degree()) {
		if (n == 0 && k == 0) continue;
		// Obtain P^n Q^k from an earlier image with one multiplication.
		const auto& previous = n > 0 ? image(n - 1, k) : image(0, k - 1);
		const auto& factor = n > 0 ? P_terms : Q_terms;
		auto& current = image(n, k);
		for (const auto& u : previous.terms())
			for (const auto& v : factor)
				if (current.contains(u.n + v.n, u.k + v.k))
					current(u.n + v.n, u.k + v.k) += u.value * v.value;
		const integer coefficient = residual(n, k);
		if (coefficient == 0) continue;
		result(n, k) = coefficient;
		for (const auto& term : current.terms())
			residual(term.n, term.k) -= coefficient * term.value;
	}
	if (!residual.terms().empty()) throw std::logic_error("composition inverse has residual");
	return result;
}

// Canonical factorization gives M = H + S odot H and H = x + t + I odot M,
// where odot substitutes odd/even parts in x. Thus I is the difference of
// inverse sums for (x + S_odd, t + S_even) and (M_odd, M_even).
inline meander_series derive_irreducibles(const meander_series& all) {
	if (all(0, 0) != 0) throw std::invalid_argument("M_(0,0) must be zero");
	const int max_n = all.limit(), max_transverse = all.max_transverse();
	meander_series P_snake(max_n, max_transverse), Q_snake(max_n, max_transverse),
		P_all(max_n, max_transverse), Q_all(max_n, max_transverse);
	P_snake(1, 0) = Q_snake(0, 1) = 1;
	for (const auto& term : iterated_snakes(max_n, max_transverse).terms())
		(term.n % 2 ? P_snake : Q_snake)(term.n, term.k) += term.value;
	for (const auto& term : all.terms())
		(term.n % 2 ? P_all : Q_all)(term.n, term.k) = term.value;
	auto result = inverse_sum(P_snake, Q_snake);
	const auto all_inverse = inverse_sum(P_all, Q_all);
	for (const auto [n, k] : result.orders_by_degree()) {
		result(n, k) -= all_inverse(n, k);
		if (result(n, k) < 0)
			throw std::logic_error("negative irreducible coefficient: check the all-meander input");
	}
	return result;
}

namespace series_detail {

using polynomial = std::vector<integer>;
using rows = std::array<polynomial, 5>;

inline polynomial product(const polynomial& A, const polynomial& B) {
	const int size = static_cast<int>(A.size());
	polynomial result(size);
	for (int i = 0; i < size; ++i)
		if (A[i] != 0)
			for (int j = 0; j < size - i; ++j)
				if (B[j] != 0) result[i + j] += A[i] * B[j];
	return result;
}

inline polynomial reciprocal(const polynomial& A) {
	if (A[0] != 1) throw std::invalid_argument("unit constant required");
	polynomial result(A.size());
	result[0] = 1;
	for (int i = 1; i < static_cast<int>(A.size()); ++i)
		for (int j = 1; j <= i; ++j) result[i] -= A[j] * result[i - j];
	return result;
}

inline polynomial power(const polynomial& A, int exponent) {
	polynomial result(A.size());
	result[0] = 1;
	for (int i = 0; i < exponent; ++i) result = product(result, A);
	return result;
}

template <size_t size>
inline polynomial evaluate(const int (&coefficients)[size], const polynomial& A) {
	polynomial result(A.size());
	result[0] = coefficients[size - 1];
	for (size_t i = size - 1; i-- > 0;) {
		result = product(result, A);
		result[0] += coefficients[i];
	}
	return result;
}

// Cyclic rerooting gives row 2 = row 1' / 4 and row 4 = row 3' / 8.
inline void differentiate_odd_rows(rows& result, int max_k) {
	for (int n : {2, 4})
		for (int k = 0; k <= max_k; ++k) {
			integer numerator = (k + 1) * result[n - 1][k + 1];
			if (numerator % (2 * n) != 0) throw std::logic_error("nonintegral row derivative");
			result[n][k] = numerator / (2 * n);
		}
}

inline unsigned totient(unsigned n) {
	unsigned result = n;
	for (unsigned p = 2; p * p <= n; ++p)
		if (n % p == 0) {
			while (n % p == 0) n /= p;
			result -= result / p;
		}
	if (n > 1) result -= result / n;
	return result;
}

} // namespace series_detail

// Exact ALL-meander rows M_0,...,M_4. Closed transverse rows 2 and 4 are
// C_(n,k) = (n + k) * M_(n-1,k) / n, so their transfer states can be omitted.
inline series_detail::rows all_meander_rows(int max_k) {
	if (max_k < 0) throw std::invalid_argument("negative row cutoff");
	const int size = max_k + 2; // Differentiation needs one extra odd-row term.
	series_detail::rows result;
	for (auto& row : result) row.resize(size);
	for (int k = 1; k < size; ++k) result[0][k] = 1;

	// M_1(t) = 1 / ((1 - 2t) sqrt(1 - 4t)).
	integer central_binomial = 1;
	result[1][0] = 1;
	for (int k = 1; k < size; ++k) {
		central_binomial *= 2 * (2 * k - 1);
		central_binomial /= k;
		result[1][k] = 2 * result[1][k - 1] + central_binomial;
	}

	// For w = t + w^2, M_3 = 2Q(w) / ((1 - 2w)^9 (1 - 2w + 2w^2)^5).
	series_detail::polynomial w(size), A(size), B(size);
	w[1] = 1;
	for (int k = 2; k < size; ++k)
		for (int j = 1; j < k; ++j) w[k] += w[j] * w[k - j];
	constexpr int Q[] = {1, -16, 120, -556, 1767, -4036, 6716, -7976, 6159,
		-1860, -2058, 3236, -2070, 600, 16, -48, 6};
	const auto w2 = series_detail::product(w, w);
	A[0] = B[0] = 1;
	for (int k = 1; k < size; ++k) {
		A[k] = -2 * w[k];
		B[k] = -2 * w[k] + 2 * w2[k];
	}
	const auto denominator = series_detail::product(series_detail::power(A, 9), series_detail::power(B, 5));
	result[3] = series_detail::product(series_detail::evaluate(Q, w), series_detail::reciprocal(denominator));
	for (auto& coefficient : result[3]) coefficient *= 2;
	series_detail::differentiate_odd_rows(result, max_k);
	return result;
}

inline series_detail::rows irreducible_meander_rows(int max_k) {
	if (max_k < 0) throw std::invalid_argument("negative row cutoff");
	const int size = max_k + 2;
	series_detail::polynomial T(size), A(size), B(size);
	// (1 + t)(T - T^2) = t; all subsequent operations are integral.
	for (int k = 1; k < size; ++k) {
		T[k] = k == 1 ? integer(1) : -T[k - 1];
		for (int j = 1; j < k; ++j) T[k] += T[j] * T[k - j];
		for (int j = 1; j < k - 1; ++j) T[k] += T[j] * T[k - 1 - j];
	}
	series_detail::rows result;
	for (auto& row : result) row.resize(size);
	// I_1 = 2T^4; I_3 = T^5 P(T) / ((1 - 2T)^5 (1 - 2T + 2T^2)).
	result[1] = series_detail::power(T, 4);
	for (auto& coefficient : result[1]) coefficient *= 2;
	constexpr int P[] = {16, -144, 576, -1308, 1728, -992, -752, 2100, -2000, 992, -224, 12};
	const auto T2 = series_detail::product(T, T);
	A[0] = B[0] = 1;
	for (int k = 1; k < size; ++k) {
		A[k] = -2 * T[k];
		B[k] = -2 * T[k] + 2 * T2[k];
	}
	const auto numerator = series_detail::product(series_detail::power(T, 5), series_detail::evaluate(P, T));
	const auto denominator = series_detail::product(series_detail::power(A, 5), B);
	result[3] = series_detail::product(numerator, series_detail::reciprocal(denominator));
	series_detail::differentiate_odd_rows(result, max_k);
	return result;
}

// Complete the unenumerated boundary from fixed rows, totient formulas and
// 2n I_(n,k) = (k + 1) I_(n-1,k+1). Compare wherever inversion also supplies a value.
inline meander_series complete_irreducibles(const meander_series& irreducibles) {
	const int max_n = irreducibles.limit();
	const auto rows = irreducible_meander_rows(max_n / 2);
	meander_series result(max_n);
	for (int n = 0; n <= max_n; ++n)
		for (int k = 0; result.contains(n, k); ++k) {
			if (n < 1 || k < 3) result(n, k) = 0;
			else if (n <= 4) result(n, k) = rows[n][k];
			else if (k == 3)
				result(n, k) = n % 2 == 0 ? series_detail::totient(n / 2 + 4) - 2 : 0;
			else if (n % 2 == 1 && k == 4)
				result(n, k) = integer((n + 1) / 2) * (series_detail::totient((n + 1) / 2 + 4) - 2);
			else if (n % 2 == 0) {
				if (!irreducibles.contains(n, k)) throw std::logic_error("missing even boundary coefficient");
				result(n, k) = irreducibles(n, k);
			}
			else {
				if (!irreducibles.contains(n + 1, k - 1)) throw std::logic_error("missing odd boundary source");
				integer numerator = 2 * (n + 1) * irreducibles(n + 1, k - 1);
				if (numerator % k != 0) throw std::logic_error("nonintegral cyclic conversion");
				result(n, k) = numerator / k;
			}
			if (irreducibles.contains(n, k) && result(n, k) != irreducibles(n, k))
				throw std::logic_error("irreducible formula mismatch");
		}
	return result;
}

} // namespace irreducible_meander_enumeration
