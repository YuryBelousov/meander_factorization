#pragma once
#include <map>
#include <set>
#include <vector>
#include <cmath>

#include "meander_calculation/integer_encoder.h"

/*
This file contains the class needed to work with two-dimensional partitions.
The main function is get_all_two_dim_partition(n, k) which returns a vector of all partitions of (n, k).


				ALGORITHM DESCRIPTION
First, note that we want to find all partitions of the pair (n, k), i.e., all sets {(a_1, b_1), ..., (a_q, b_q)} 
with the following properties:
1) sum a_i = n;
2) sum b_i = k;
3) |a_i| + |b_i| > 0 for all i.
Let f(n, k, a, b) be the function that returns all the partitions of (n, k) where each element is not greater
than (a, b). Then f(n, k, a, b) can be constructed recursively. For each pair (c, d) such that (c,d) <= (a, b) we 
take the union of {(c,d)} and each partition from f(n-c, k-d, c,d). Here we have to make sure, that f(n-c, k-d, c,d)
is not empty (i.e., there are partitions of (n-c, k-d) where every element is greater than (c,d)). But the only case
where this set is empty is when c == 0 and n != 0. 

To speed up the computation within the algorithm work, we encode a pair of integers into a single number 
(see integer_encoder.h).
*/

using namespace int_encoder;

using two_dim_partition = std::vector<int_pair_encoded>;
using two_dim_partition_container = std::vector<two_dim_partition>;

//This function returns number of odd elements (i.e., of type (2n+1, k)) in a given partition
int partition_parity(const two_dim_partition& partition) noexcept {
	int res = 0;
	for (const auto& elt : partition) 
		res += (decode_first_coord(elt) % 2);
	return res;
}

bool is_partition_of_key(int_pair_encoded key, const two_dim_partition& partition) noexcept {
	int_pair_encoded sum = key;
	for (const auto& elt : partition) 
		sum = code_of_difference(sum, elt);
	return sum == encode(0,0);
}

class two_dimensional_partition_calculator {
private:
	//This container is used to store partitions that have already been calculated (see algorithm description)
	std::map<int_pair_encoded, two_dim_partition_container> all_two_dim_partitions;

private:
	//Add new_element to all partition
	two_dim_partition_container add_element_to_partitions(int_pair_encoded new_element, two_dim_partition_container partitions) const noexcept {
		if (partitions.empty())
			return { {new_element} };
		for (auto& partition : partitions)
			partition.push_back(new_element);
		return partitions;
	}

	//Find all possible values to reduce the decrement that are less than or equal to reduction_upper_boundary.
	std::vector<int_pair_encoded> all_possible_reduction_elements(int_pair_encoded decrement, int_pair_encoded reduction_upper_boundary) const noexcept { 
		auto rub_x = decode_first_coord(reduction_upper_boundary),
			rub_y = decode_second_coord(reduction_upper_boundary),
			d_x = decode_first_coord(decrement),
			d_y = decode_second_coord(decrement);
		std::vector<int_pair_encoded> results;
		for (int i = std::min(rub_x, d_x); i >= 0; --i) {
			for (int j = d_y; j >= 0; --j) {
				/*
				To make sure that (i, j) is a possible element for reduction, we need to check that
				1) (i,j) <= reduction_upper_boundary,
				2) decrement can be reduced to zero by elements not greater than (i, j).
				*/
				if ((i < rub_x || (i == rub_x && j <= rub_y))
					&& (i != 0 || (i == 0 && j != 0 && d_x == 0 ))) {
					results.push_back(encode(i, j));
				}
			}
		}
		return results;
	}

	const two_dim_partition_container& get_all_two_dim_partition_recursion(int_pair_encoded key_element, int_pair_encoded reduction_upper_boundary) noexcept {
		//Check if we have already found all the partitions
		auto key_and_boundary_coded = encode_2(key_element, reduction_upper_boundary);
		if (all_two_dim_partitions.contains(key_and_boundary_coded))
			return all_two_dim_partitions[key_and_boundary_coded];
		all_two_dim_partitions[key_and_boundary_coded];
		for (const auto& reducer : all_possible_reduction_elements(key_element, reduction_upper_boundary)) {
			auto& ref = all_two_dim_partitions[key_and_boundary_coded];
			auto temp = add_element_to_partitions(reducer,
				get_all_two_dim_partition_recursion(code_of_difference(key_element, reducer), reducer));
			ref.insert(ref.cend(),
				std::make_move_iterator(temp.cbegin()),
				std::make_move_iterator(temp.cend()));
		}
		return all_two_dim_partitions[key_and_boundary_coded];
	}


public:
	two_dimensional_partition_calculator() {
		all_two_dim_partitions[encode(0, 0)] = {};
	}

	const two_dim_partition_container& get_all_two_dim_partition(int n, int k, bool use_allowed = false) noexcept {
		return get_all_two_dim_partition_recursion(encode(n, k), encode(n, k));
	}

};