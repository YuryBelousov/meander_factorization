#pragma once
#include <thread>
#include <atomic>
#include <set>
#include <mutex>

#include "irreducible_meander_enumeration/submeander.h"

/*
This file contains functions to perform a complete enumeration of irreducible meanders.
The algorithm is straightforward, its description can be found in the file submeander.h.

The implementation of the algorithm is as follows. The `brute_force` function distributes
the work among the threads in such a way as to minimise the maximum running time among
all the threads (to approximate the running time, we use the estimate obtained from prior
computations involving meanders of complexity 30). Each thread executes the function
`brute_forse_several_steps` that sequentially enumerates all almost irreducible meanders 
starting from a given permutation. For each step of the enumeration, the 
`brute_force_one_step` function performs a depth-first search.
*/


constexpr int number_of_threads = 7;

std::atomic<unsigned long long> irreducible_meander_numbers[max_n + 1][max_n + 1];
std::mutex cout_mutex;

void brute_forse_one_step(submeander& m) noexcept {
	int av_points[max_n / 2];
	std::memcpy(av_points, m.get_avaliable_points(), sizeof(av_points));
	for (int i = 0; i < max_n / 2; ++i) {
		if (av_points[i] == meanless_value)
			break;
		m.go_throught_point(av_points[i]);
		if (m.finished()) {
			int n = m.get_n();
			int k = m.get_cups();
			++irreducible_meander_numbers[n - 2 * k - 1][k - 3];
		}
		else
			brute_forse_one_step(m);
		m.step_back();
	}
}


void brute_forse_several_steps(const std::set<std::vector<int>>& code_starts) noexcept {
	auto start = std::chrono::steady_clock::now();

	for (const auto& points : code_starts) {
		submeander m(points);
		brute_forse_one_step(m);
	}
	std::lock_guard<std::mutex> guard(cout_mutex);
	for (const auto& points : code_starts) {
		std::cout << "{ ";
		for (auto pt : points)
			std::cout << pt << ", ";
		std::cout << "}, ";
	}
	std::cout << "\n";
}


void brute_forse() noexcept {
	std::thread threads[number_of_threads];

	std::set<std::vector<int>> starting_points[number_of_threads];
	for (int i = 4; i < max_n - 2; i += 2)
		starting_points[0].insert({ i });
	for (int i = 3; i < max_n - 2; i += 2)
		starting_points[0].insert({ 2, i });

	//This function approximates the number of corresponding meanders.
	//The approximation is based on the calculation for max_n=30. 
	//For n=38 it does not work too well. Apparently, we should use some other way to divide the work between threads.
	auto weight_approximation = [](const std::vector<int>& points, int max_ind) -> double {
		double t = static_cast<double>(points[0] - 2) / static_cast<double>(max_ind - 4);
		double w = 0.2 + t * (-1.08 + t * (2.88 + t * (-3.46 + t * 1.59)));
		if (points.size() != 1) {
			t = static_cast<double>(points[1] - 3) / static_cast<double>(max_ind - 3);
			w *= 0.27 + t * (-1.75 + t * (4.81 + t * (-5.69 + t * 2.46)));
		}
		return w;
		};

	//Try to reduce the maximal weight among all threads
	bool can_reduce = true;
	while (can_reduce) {
		can_reduce = false;
		double max_sum_value = -1.,
			min_sum_value = 10.;
		int max_ind = -1,
			min_ind = -1;
		for (int i = 0; i < number_of_threads; ++i) {
			double temp_value = 0.;
			for (const auto& e : starting_points[i])
				temp_value += weight_approximation(e, max_n);
			if (max_sum_value < temp_value) {
				max_sum_value = temp_value;
				max_ind = i;
			}
			if (min_sum_value > temp_value) {
				min_sum_value = temp_value;
				min_ind = i;
			}
		}
		double min_value = max_n;
		std::vector<int> min_element = {};
		for (const auto& elt_in_max : starting_points[max_ind]) {
			if (min_value > weight_approximation(elt_in_max, max_n)) {
				min_value = weight_approximation(elt_in_max, max_n);
				min_element = elt_in_max;
			}
		}
		if (min_sum_value + min_value < max_sum_value) {
			can_reduce = true;
			starting_points[min_ind].insert(min_element);
			starting_points[max_ind].erase(min_element);
		}
	}

	for (int i = 0; i < number_of_threads; ++i)
		threads[i] = std::thread(brute_forse_several_steps, std::cref(starting_points[i]));

	for (int i = 0; i < number_of_threads; ++i)
		threads[i].join();
}