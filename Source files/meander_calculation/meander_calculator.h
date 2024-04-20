#pragma once
#include <vector>
#include <algorithm>

#include "meander_calculation/two_dimensional_partitions.h"
#include "meander_calculation/binomial_coefficient_calculator.h"

/*
This file contains the class for calculating the numbers of different classes of meanders 
(except for the irreducible ones, for which the files `submeander.h` and `brute_force_functions.h`
are used). The calculations naively follow the formulas given in the corresponding article.
*/


using namespace int_encoder;

class meander_calculator {
	using get_some_meander_number = big_int (meander_calculator::* const)(int, int );
	using container_for_numbers = std::unordered_map<int_pair_encoded, big_int>;

private:	
	container_for_numbers irreducible_meanders;			
	container_for_numbers iterated_irreducible_meanders;
	container_for_numbers iterated_snakes;

	container_for_numbers iter_snakes_started;
	container_for_numbers iter_irreducible_started;

	two_dimensional_partition_calculator partition_calculator;
	binomial_coefficient_calculator binom_calculator;

public:	
	meander_calculator() {
		iterated_snakes[encode(0, 0)] = 
			iterated_snakes[encode(0, 1)] = 
			iterated_snakes[encode(1, 0)] = 0;
	}

private:
	big_int tree_calculation(container_for_numbers& calculating_type, get_some_meander_number root_type, get_some_meander_number child_type, int n, int k) noexcept {
		auto index = encode(n, k);
		if (calculating_type.contains(index))
			return calculating_type[index];

		big_int res = 0;
		for (int child_n = 0; child_n <= n; ++child_n) {
			for (int child_k = 0; child_k <= k; ++child_k) {
				for (const auto& mu : partition_calculator.get_all_two_dim_partition(child_n, child_k)) {
					int odd_childs = partition_parity(mu),
						even_childs = mu.size() - odd_childs,
						root_n = n - child_n + odd_childs,
						root_k = k - child_k + even_childs;

					if (odd_childs > root_n || even_childs > root_k)
						continue;
					auto local_res = (this->*root_type)(root_n, root_k);
					if (local_res.is_zero()) 
						continue;

					int n_residual = root_n,
						k_residual = root_k,
						cur_order = -1,
						s_cur_order = 0;
					for (const auto& child_order : mu) {
						if (cur_order == -1)
							cur_order = child_order;
						if (cur_order == child_order)
							++s_cur_order;
						else {
							if (decode_first_coord(cur_order) % 2 == 1) {
								local_res *= binom_calculator.get_binomial_coef(n_residual, s_cur_order);
								n_residual -= s_cur_order;
							}
							else {
								local_res *= binom_calculator.get_binomial_coef(k_residual, s_cur_order);
								k_residual -= s_cur_order;
							}
							cur_order = child_order;
							s_cur_order = 1;
						}
						auto temp = (this->*child_type)(decode_first_coord(cur_order), decode_second_coord(cur_order));
						local_res *= temp;
						if (local_res.is_zero()) 
							break;
					}
					if (!local_res.is_zero()) {
						if (decode_first_coord(cur_order) % 2 == 1)
							local_res *= binom_calculator.get_binomial_coef(n_residual, s_cur_order);
						else
							local_res *= binom_calculator.get_binomial_coef(k_residual, s_cur_order);
						res += local_res;
					}
				}
			}
		}
		calculating_type[index] = res + (this->*root_type)(n, k);
		return calculating_type[index];
	}

	big_int get_iter_snakes_started_number(int n, int k) noexcept {
		return tree_calculation(iter_snakes_started, 
			&meander_calculator::get_iterated_snakes_number, 
			&meander_calculator::get_iter_irreducible_started_number,
			n, k);
	}

	big_int get_iter_irreducible_started_number(int n, int k) noexcept{
		return tree_calculation(iter_irreducible_started,
			&meander_calculator::get_iterated_ireducible_meanders_number,
			&meander_calculator::get_iter_snakes_started_number,
			n, k);
	}

public:
	big_int get_iterated_ireducible_meanders_number(int n, int k) noexcept {
		if (k < 3 || n < 1) 
			return 0;
		
		auto index = encode(n, k);
		if (iterated_irreducible_meanders.contains(index))
			return iterated_irreducible_meanders[index];
	

		big_int res = 0;
		for (int child_n = 0; child_n < n; ++child_n) {
			for (int child_k = 0; child_k < k-2; ++child_k) {
				for (const auto& mu : partition_calculator.get_all_two_dim_partition(child_n, child_k)) {
					int even_childs = partition_parity(mu),
						odd_childs = mu.size() - even_childs,
						root_n = n - child_n + odd_childs - mu.size(),
						root_k = k - child_k + even_childs - 3*mu.size();

					if (odd_childs > root_n || even_childs > root_k)
						continue;
					auto local_res = get_ireducible_meanders_number(root_n, root_k);
					if (local_res.is_zero())
						continue;

					int n_residual = root_n,
						k_residual = root_k,
						cur_order = -1,
						s_cur_order = 0;
					for (const auto& child_order : mu) {
						if (cur_order == -1)
							cur_order = child_order;
						if (cur_order == child_order)
							++s_cur_order;
						else {
							if (decode_first_coord(cur_order) % 2 == 0) {
								local_res *= binom_calculator.get_binomial_coef(n_residual, s_cur_order);
								n_residual -= s_cur_order;
							}
							else {
								local_res *= binom_calculator.get_binomial_coef(k_residual, s_cur_order);
								k_residual -= s_cur_order;
							}
							cur_order = child_order;
							s_cur_order = 1;
						}
						auto temp = get_iterated_ireducible_meanders_number(decode_first_coord(cur_order) + 1, decode_second_coord(cur_order) + 3);
						local_res *= temp;
						if (local_res.is_zero())
							break;
					}
					if (!local_res.is_zero()) {
						if (decode_first_coord(cur_order) % 2 == 0)
							local_res *= binom_calculator.get_binomial_coef(n_residual, s_cur_order);
						else
							local_res *= binom_calculator.get_binomial_coef(k_residual, s_cur_order);
						res += local_res;
					}
				}
			}
		}
		iterated_irreducible_meanders[index] = res + get_ireducible_meanders_number(n, k);
		return iterated_irreducible_meanders[index];
	}

	big_int get_iterated_snakes_number(int n, int k) noexcept {
		auto index = encode(n, k);
		if (iterated_snakes.contains(index))
			return iterated_snakes[index];

		int delta = n % 2;
		big_int res = 0;
		for (int root_n = 2 - delta; root_n <= n; root_n += 2)
			for (int root_k = 0; root_k <= k; ++root_k) {
				if (root_n == 1 && root_k == 0 
					|| n == root_n && k == root_k)
					continue;
				for (const auto& mu : partition_calculator.get_all_two_dim_partition((n - root_n) / 2, k - root_k)) {
					if (mu.size() > root_n)
						continue;
					big_int local_res = binom_calculator.get_binomial_coef(root_n + root_k, root_k);
					int n_residual = root_n,
						cur_order = -1,
						s_cur_order = 0;
					for (const auto& child_order : mu) {
						if (cur_order == -1)
							cur_order = child_order;
						if (cur_order == child_order)
							++s_cur_order;
						else {
							cur_order = child_order;
							local_res *= binom_calculator.get_binomial_coef(n_residual, s_cur_order);
							n_residual -= s_cur_order;
							s_cur_order = 1;
						}
						local_res *= get_iterated_snakes_number(2 * decode_first_coord(cur_order) + 1, decode_second_coord(cur_order)) / 2;
					}
					local_res *= binom_calculator.get_binomial_coef(n_residual, s_cur_order);
					res += local_res;
				}
			}
		iterated_snakes[index] = (1 + delta) * (res + binom_calculator.get_binomial_coef(n + k, k));
		return iterated_snakes[index];
	}

	big_int get_ireducible_meanders_number(int n, int k) noexcept {
		auto index = encode(n, k);
		if (irreducible_meanders.contains(index))
			return irreducible_meanders[index];
		return 0;
	}

	big_int get_meander_number(int n, int k) noexcept {
		if (n + k == 1)
			return 1;
		return get_iter_snakes_started_number(n, k) + get_iter_irreducible_started_number(n, k);
	}

	void read_irreducible_meander_numbers(std::istream& file_with_meanders) {
		int max_n = 0, max_k = 0;
		file_with_meanders >> max_n >> max_k;
		for (int n = 1; n <= max_n; ++n)
			for (int k = 3; k <= max_k; ++k) {
				big_int irr_meand_num = 0;
				file_with_meanders >> irr_meand_num;
				if (irr_meand_num != 0) 
					irreducible_meanders[encode(n, k)] = irr_meand_num;
			}
	}
};
