#pragma once
#include <vector>
#include <array>
#include <algorithm>
#include <optional>

#include "meander_calculation/binomial_coefficient_calculator.hpp"

/*
This file defines a class for calculating various classes of meanders (excluding irreducible meanders,
which are handled separately in submeander.h). The calculations are performed using truncated generating
functions. For mathematical details, see the corresponding paper.

Template parameter:
	- number_type: the numeric type used for internal calculations (e.g., long long). It should support
	basic arithmetic operations (+, -, *, assignment) and allow implicit conversion from int.

The main functions are :
	-read_irreducible_meander_numbers(file) :
	Reads the numbers of irreducible meanders from a file.
	Note: Reading starts from n = 1 and k = 3, as there are no irreducible meanders for n < 1 or k < 3.

	- calculate_all(n, k) :
	Computes the numbers of meanders, iterated snakes, and irreducible meanders up to the given order(n, k).
	Note: During calculation, only the known part of the irreducible meanders is used. Therefore, results
	for higher orders will be incorrect.

	- get_ * something * (n, k) :
	Returns the number of corresponding meanders of order(n, k).
	Note: If the requested number has not been calculated yet, the function recalculates all numbers up to
	the given order.

Recommendations(to avoid recomputations) :
	1) If you need all numbers up to a given order, call calculate_all first.
	2) If you need specific numbers up to order(n, k), call the appropriate get_*...*(n, k) function.
	After that, you can call the getters with smaller parameters without recomputation.
*/

template <typename number_type>
class meander_calculator {
	using container_for_opt_numbers = std::vector<std::vector<std::optional<number_type>>>;

	using numbers_1D_container = std::vector<number_type>;
	using numbers_2D_container = std::vector<numbers_1D_container>;
	using numbers_3D_container = std::vector<numbers_2D_container>;

private:
	binomial_coefficient_calculator <number_type> binom_calculator;

	container_for_opt_numbers irreducible_meanders;
	container_for_opt_numbers iterated_irreducible_meanders;
	container_for_opt_numbers iterated_snakes;

	container_for_opt_numbers meander_numbers;

private:
	// Ensures that numbers[n][k] is accessible. Expands the container if necessary.
	void ensure_capacity(container_for_opt_numbers& numbers, size_t n, size_t k, bool is_jagged) noexcept {
		if (is_jagged) {
			if (numbers.size() <= n)
				numbers.resize(1 + std::max(n, 2 * numbers.size()));
			for (size_t i = 0; i <= n; ++i) {
				auto& row = numbers[i];
				size_t need_columns = k + (n - i) / 2;
				if (row.size() <= need_columns)
					row.resize(1 + std::max(need_columns, 2 * row.size()), std::nullopt);
			}
		}
		else {
			if (numbers.size() <= n)
				numbers.resize(1 + std::max(n, 2 * numbers.size()), std::vector<std::optional<number_type>>(k + 1, std::nullopt));
			for (auto& row : numbers) {
				if (row.size() <= k)
					row.resize(1 + k, std::nullopt);
			}
		}
	}

	// Calculates sum_{i,j} A[i][j]*A[n-i][k-j]. Assumes all numbers up to (n, k) have been precomputed.
	number_type quadric_term(const container_for_opt_numbers& numbers, size_t n, size_t k) const noexcept {
		number_type sum(0);
		for (size_t i = 0; i <= n; ++i) {
			const size_t j_min = (i < 2 ? 2 - i : 0);	// numbers[i][j], where i+j<2 are always 0.
			const size_t j_max = (i == n ? k : k + 1);	// Do not use numbers[n][k] itself.
			for (size_t j = j_min; j < j_max; ++j)
				sum += numbers[i][j].value() * numbers[n - i][k - j].value();
		}
		return sum;
	}

	// Calculates sum_{i1+i2+i3=n, j1+j2+j3=k} A[i1][j1]*A[i2][j2]*A[i3][j3]. Assumes all numbers up to (n, k) have been precomputed.
	number_type cubic_term(const container_for_opt_numbers& numbers, size_t n, size_t k) const noexcept {
		number_type sum = 0;
		for (size_t i1 = 0; i1 <= n; ++i1) {
			const size_t j1_min = (i1 < 2 ? 2 - i1 : 0);	// numbers[i1][j1], where i1+j1<2 are always 0.
			const size_t j1_max = (i1 == n ? k : k + 1);	// Do not use numbers[n][k] itself.
			for (size_t j1 = j1_min; j1 < j1_max; ++j1)
				for (size_t i2 = 0; i2 <= n - i1; ++i2) {
					const size_t j2_min = (i2 < 2 ? 2 - i2 : 0);
					const size_t j2_max = (i2 == n && j1 == 0 ? k : k + 1 - j1);
					for (size_t j2 = j2_min; j2 < j2_max; ++j2)
						sum += numbers[i1][j1].value() * numbers[i2][j2].value() * numbers[n - i1 - i2][k - j1 - j2].value();
				}
		}
		return sum;
	}

	// Calculates the truncated product of two bivariate polynomials A(x, t) and B(x, t).
	// Assumes A and B have the same truncation.
	numbers_2D_container product(const numbers_2D_container& A, const numbers_2D_container& B) const noexcept {
		numbers_2D_container result(A.size());
		for (size_t i = 0; i < A.size(); ++i)
			result[i] = numbers_1D_container(A[i].size(), 0);
		for (size_t i1 = 0; i1 < A.size(); ++i1)
			for (size_t j1 = 0; j1 < A[i1].size(); ++j1)
				if (A[i1][j1] != 0)
					for (size_t i2 = 0; i1 + i2 < B.size(); ++i2)
						for (size_t j2 = 0; j1 + j2 < B[i1 + i2].size() && j2 < B[i2].size(); ++j2)
							if (B[i2][j2] != 0)
								result[i1 + i2][j1 + j2] += A[i1][j1] * B[i2][j2];
		return result;
	}

	// Computes truncated powers of a bivariate polynomial A(x, t) up to max_power.
	numbers_3D_container get_powers(const numbers_2D_container& A, size_t max_power) const noexcept {
		numbers_3D_container A_powers(1 + max_power);
		A_powers[0] = numbers_2D_container(1, numbers_1D_container(1, 1));
		A_powers[1] = A;
		for (size_t power = 2; power <= max_power; ++power)
			A_powers[power] = product(A_powers[power - 1], A);
		return A_powers;
	}

	// Calculates a truncated version of A(x + B1(x, t), t + B2(x, t)), 
	// where B1(x, t) is the odd part and B2(x, t) is the even part of B(x, t).
	// Assumes A and B have the same truncation.
	numbers_2D_container compose(const numbers_2D_container& A, const numbers_2D_container& B) const noexcept {
		numbers_2D_container B1(A.size());
		numbers_2D_container B2(A.size());

		// Create odd and even parts of B.
		for (size_t i = 0; i < A.size(); ++i) {
			B1[i].resize(A[i].size(), 0);
			B2[i].resize(A[i].size(), 0);
			for (size_t j = 0; j < A[i].size(); ++j) {
				if (number_type b_ij = B[i][j]; b_ij != 0)
					(i % 2 == 0 ? B2[i][j] : B1[i][j]) = b_ij;
			}
		}

		if (A.size() >= 2) B1[1][0] = 1;		// x + B1(x,t).
		if (A[0].size() >= 2) B2[0][1] = 1;		// t + B2(x,t).

		numbers_3D_container B1_powers = get_powers(B1, A.size());
		numbers_3D_container B2_powers = get_powers(B2, A[0].size());

		numbers_2D_container result(A.size());
		for (size_t i = 0; i < A.size(); ++i)
			result[i] = numbers_1D_container(A[i].size(), 0);

		for (size_t i = 0; i < A.size(); ++i) {
			for (size_t j = 0; j < A[i].size(); ++j)
				if (number_type c = A[i][j]; c != 0) {
					// Add contribution of A[i][j] * (B1^i * B2^j) to the result.
					numbers_2D_container term = (i == 0 ? B2_powers[j] : (j == 0 ? B1_powers[i] : product(B1_powers[i], B2_powers[j])));
					for (size_t x = 0; x < A.size(); ++x)
						for (size_t t = 0; t < A[x].size(); ++t)
							result[x][t] += c * term[x][t];
				}
		}
		return result;
	}

	// Calculates the number of iterated snakes up to a given order (n, k).
	void calculate_iterated_snakes_up_to(size_t n, size_t k, bool is_jagged) noexcept {
		// Solves the third-degree recurrence relation for the iterated snake generating function.
		ensure_capacity(iterated_snakes, n, k, is_jagged);
		for (size_t i = 0; i <= n; ++i) {
			number_type prev_cubic = 0;
			const size_t j_min = (i < 2 ? 2 - i : 0);
			const size_t j_max = (is_jagged ? k + (n - i) / 2 : k);
			for (size_t j = j_min; j <= j_max; ++j) {
				if (iterated_snakes[i][j].has_value()) {
					prev_cubic = cubic_term(iterated_snakes, i, j);
					continue;
				}

				// 1) Cubic_term part: IS(x,t)^3 * (1 + t)
				number_type S3_0 = cubic_term(iterated_snakes, i, j),
					S3_1 = (j == 0 ? 0 : prev_cubic),
					S2 = 0, S1 = 0;

				prev_cubic = S3_0;	// Save previously calculated cubic term.

				// 2) Quadratic: IS(x,t)^2 * B(x,t)
				static constexpr std::array<std::array<int, 3>, 5> B{ { {0,0,2}, {0,1,-8}, {0,2,7}, {1,0,-5}, {1,1,5} } };
				for (const auto& monomial : B) {
					auto [deg_x, deg_t, coef] = monomial;
					if (i >= deg_x && j >= deg_t)
						S2 += coef * quadric_term(iterated_snakes, i - deg_x, j - deg_t);
				}

				// 3) Linear: IS(x,t) * C(x,t)
				static constexpr std::array<std::array<int, 3>, 8> C{
				{ {0,1, -3}, {0,2, -5}, {0,3, 11}, {1,0,  1}, {1,1,-16}, {1,2, 18}, {2,0, -7}, {2,1,  7} } };
				for (const auto& monomial : C) {
					auto [deg_x, deg_t, coef] = monomial;
					if (i < deg_x || j < deg_t || (i == deg_x && j == deg_t)) continue;
					S1 += coef * iterated_snakes[i - deg_x][j - deg_t].value();
				}

				// 4) Inhomogeneity: D_{n,k}
				constexpr int shift = 10;
				number_type S0(0);
				switch ((i << shift) | j) {
				case (0 << shift) | 2:	S0 = -1;	break;
				case (0 << shift) | 3:	S0 = 2;		break;
				case (0 << shift) | 4:	S0 = 5;		break;
				case (1 << shift) | 1:	S0 = -4;	break;
				case (1 << shift) | 2:	S0 = -3;	break;
				case (1 << shift) | 3:	S0 = 13;	break;
				case (2 << shift) | 0:	S0 = -1;	break;
				case (2 << shift) | 1:	S0 = -8;	break;
				case (2 << shift) | 2:	S0 = 11;	break;
				case (3 << shift) | 0:	S0 = -3;	break;
				case (3 << shift) | 1:	S0 = 3;		break;
				default:							break;
				}

				// 5) Solve the master recurrence equation: 
				iterated_snakes[i][j] = S3_0 - S3_1 - S2 - S1 - S0;
			}
		}
	}

	// Calculates the number of iterated irreducible meanders up to given complexity.
	void calculate_iterated_irreducibe_meanders_up_to(size_t n, size_t k) noexcept {
		// For the calculation it uses an equation for generating function: II(x,t) = A(x + II1(x,t), t + II2(x,t)).
		// Where A(x,t) is a generating function for the irreducible meanders.
		ensure_capacity(iterated_irreducible_meanders, n, k, true);
		if (!iterated_irreducible_meanders[n][k].has_value()) {
			numbers_2D_container res_A(n + 1), A(n + 1);

			// Prepare truncated version of A(x,t) and II(x,t)
			for (size_t i = 0; i <= n; ++i) {
				A[i].resize(1 + k + (n - i) / 2, 0);
				res_A[i].resize(1 + k + (n - i) / 2, 0);
				for (size_t j = 0; j < A[i].size(); ++j)
					A[i][j] = get_irreducible_meanders_number(i, j);
			}

			// Perform convolution until convergence.
			while (true) {
				auto temp = compose(A, res_A);
				if (temp == res_A) break;
				res_A = std::move(temp);
			}

			// Save result.
			for (size_t i = 0; i < res_A.size(); ++i)
				for (size_t j = 0; j < res_A[i].size(); ++j)
					iterated_irreducible_meanders[i][j] = std::move(res_A[i][j]);
		}
	}

public:
	meander_calculator() noexcept : iterated_snakes{ {0,0}, {0,4} }, iterated_irreducible_meanders{ {0,0}, {0,0} } {};

	// Calculates the number of meanders up to given order.
	void calculate_all(size_t n, size_t k) noexcept {
		// For the calculation it uses an truncated generating functions.
		calculate_iterated_irreducibe_meanders_up_to(n, k);
		calculate_iterated_snakes_up_to(n, k, true);

		numbers_2D_container iter_snakes_started(n + 1), iter_snakes_started_conv(n + 1),
			iter_irr_started(n + 1), iter_irr_started_conv(n + 1);

		// Prepare truncated version of all ingredients.
		for (size_t i = 0; i <= n; ++i) {
			iter_snakes_started[i].resize(1 + k + (n - i) / 2, 0);
			iter_snakes_started_conv[i].resize(1 + k + (n - i) / 2, 0);
			iter_irr_started[i].resize(1 + k + (n - i) / 2, 0);
			iter_irr_started_conv[i].resize(1 + k + (n - i) / 2, 0);
			for (size_t j = 0; j <= k + (n - i) / 2; ++j) {
				iter_snakes_started[i][j] = iter_snakes_started_conv[i][j] = get_iterated_snakes_number(i, j);
				iter_irr_started[i][j] = iter_irr_started_conv[i][j] = get_iterated_irreducible_meanders_number(i, j);
			}
		}

		// Perform convolution until convergence.
		while (true) {
			auto temp = compose(iter_snakes_started, iter_irr_started_conv);
			if (temp == iter_snakes_started_conv) break;
			iter_snakes_started_conv = std::move(temp);
			iter_irr_started_conv = compose(iter_irr_started, iter_snakes_started_conv);
		}

		// Save result.
		meander_numbers.resize(iter_snakes_started.size());
		for (size_t i = 0; i < iter_snakes_started.size(); ++i) {
			meander_numbers[i].resize(iter_snakes_started[i].size());
			for (size_t j = 0; j < iter_snakes_started[i].size(); ++j)
				if (i + j == 1) 
					meander_numbers[i][j] = 1;
				else
					meander_numbers[i][j] = iter_snakes_started_conv[i][j] + iter_irr_started_conv[i][j];
		}
	}

	number_type get_meander_number(size_t n, size_t k) noexcept {
		ensure_capacity(meander_numbers, n, k, false);
		if (!meander_numbers[n][k].has_value())
			calculate_all(n, k);
		return meander_numbers[n][k].value();
	}

	number_type get_snakes_number(size_t n, size_t k) noexcept {
		if (n + k < 2) return 0;
		return (n % 2 == 0 ? 1 : 2) * binom_calculator.get_binomial_coef(n + k, k);
	}

	number_type get_iterated_snakes_number(size_t n, size_t k) noexcept {
		ensure_capacity(iterated_snakes, n, k, false);
		if (!iterated_snakes[n][k].has_value())
			calculate_iterated_snakes_up_to(n, k, false);
		return iterated_snakes[n][k].value();
	}

	number_type get_irreducible_meanders_number(size_t n, size_t k) const noexcept {
		if (irreducible_meanders.size() <= n || irreducible_meanders[0].size() <= k)
			return 0;
		return irreducible_meanders[n][k].value_or(0);
	}

	number_type get_iterated_irreducible_meanders_number(size_t n, size_t k) noexcept {
		ensure_capacity(iterated_irreducible_meanders, n, k, false);
		if (!iterated_irreducible_meanders[n][k].has_value())
			calculate_iterated_irreducibe_meanders_up_to(n, k);
		return iterated_irreducible_meanders[n][k].value();
	}

	// Read number of irreducible meanders from a given file. Note: it starts reading from n = 1 and k = 3.
	void read_irreducible_meander_numbers(std::istream& file_with_meanders) {
		int max_n = 0, max_k = 0;
		file_with_meanders >> max_n >> max_k;
		ensure_capacity(irreducible_meanders, max_n + 1, max_k + 1, false);
		for (int n = 1; n <= max_n; ++n)
			for (int k = 3; k <= max_k; ++k) {
				number_type irr_meand_num = 0;
				file_with_meanders >> irr_meand_num;
				if (irr_meand_num != 0)
					irreducible_meanders[n][k] = irr_meand_num;
			}
	}
};