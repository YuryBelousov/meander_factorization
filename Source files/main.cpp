#include <iostream>
#include <fstream>
#include <string>
#include <time.h>

#include "boost/multiprecision/cpp_int.hpp"

#include "meander_calculation/meander_calculator.hpp"
#include "irreducible_meander_enumeration/brute_force_functions.h"

using big_int = boost::multiprecision::uint1024_t;

std::string time_passed_since(const std::clock_t& start) {
	double elapsed = double(std::clock() - start) / CLOCKS_PER_SEC;
	std::ostringstream temp_stream;
	temp_stream << std::fixed << std::setprecision(3);

	if (elapsed < 60.0)
		temp_stream << elapsed << " sec";
	else if (elapsed < 3600.0)
		temp_stream << int(elapsed / 60.) << " min " << (elapsed - int(elapsed / 60) * 60) << " sec";
	else if (elapsed < 86400.0)
		temp_stream << int(elapsed / 3600.) << " hour " << ((elapsed - int(elapsed / 3600) * 3600) / 60) << " min";
	else
		temp_stream << int(elapsed / 86400.) << " days " << ((elapsed - int(elapsed / 86400) * 86400) / 3600) << " hours";
	return temp_stream.str();
}


int main() {
	auto start = std::clock();

	// Calculate irreducible meander numbers up to n+2k <= max_n, see submeander.h
	brute_force();
	std::cout << "It takes " << time_passed_since(start)
		<< " to enumerate all irreducible meander with n + 2k < " 
		<< max_n << std::endl;
	start = std::clock();

	// Create meander calculator.
	meander_calculator<big_int> mnd_calc;
	std::ifstream irr_meander_stream("files with numbers/irreducible_meanders_38.txt");
	mnd_calc.read_irreducible_meander_numbers(irr_meander_stream);

	// Calculate non-singular iterated snakes.
	const size_t N_nonsing_iter_snakes = 500;
	mnd_calc.get_iterated_snakes_number(N_nonsing_iter_snakes, 0);
	std::cout 
		<< "It takes " << time_passed_since(start)
		<< " to enumerate all non-singular iterated snakes of total order less than " 
		<< N_nonsing_iter_snakes << std::endl;
	start = std::clock();

	// Calculate iterated snakes.
	const size_t N_iter_snakes = 100, K_iter_snakes = 30;
	mnd_calc.get_iterated_snakes_number(N_iter_snakes, K_iter_snakes);
	std::cout 
		<< "It takes " << time_passed_since(start)
		<< " to enumerate all iterated snakes of order less than (" << N_iter_snakes 
		<< ", " << K_iter_snakes << ")" << std::endl;
	start = std::clock();

	// Calculate all meanders. Note: only the known part of the irreducible meanders is used (up to n+2k = 38). 
	// Therefore, results for higher orders will be incorrect.
	const size_t pre_calculated_irreducible = 38;
	mnd_calc.calculate_all(pre_calculated_irreducible, 0);
	std::cout 
		<< "It takes " << time_passed_since(start)
		<< " to enumerate all non-singular meanders of total order up to " << pre_calculated_irreducible
		<< std::endl;

	// Save results.
	std::ofstream file_with_irreducible_meanders_out("files with numbers/irreducible_meanders_test.txt");
	for (size_t n = 0; n <= max_n; ++n)
		for (size_t k = 0; k <= max_n; ++k)
		file_with_irreducible_meanders_out << irreducible_meander_numbers[n][k] << (k == max_n ? '\n' : ' ');
	file_with_irreducible_meanders_out.close();

	std::ofstream file_with_non_singular_iter_snakes("files with numbers/iterated_snakes_non_singular_" + std::to_string(N_nonsing_iter_snakes) + ".txt");
	for (size_t n = 0; n <= N_nonsing_iter_snakes; ++n)
		file_with_non_singular_iter_snakes << mnd_calc.get_iterated_snakes_number(n, 0) << std::endl;
	file_with_non_singular_iter_snakes.close();

	std::ofstream file_with_iter_snakes("files with numbers/iterated_snakes_" + std::to_string(N_iter_snakes) + "_" + std::to_string(K_iter_snakes) + ".txt");
	for (size_t n = 0; n <= N_iter_snakes; ++n)
		for (size_t k = 0; k <= K_iter_snakes; ++k)
			file_with_iter_snakes << mnd_calc.get_iterated_snakes_number(n, k) << (k == K_iter_snakes ? '\n' : ' ');
	file_with_iter_snakes.close();

	std::ofstream file_with_meander_numbers("files with numbers/meander_numbers.txt");
	for (size_t n = 0; n <= pre_calculated_irreducible; ++n) {
		const size_t k_max = (pre_calculated_irreducible - n) / 2;
		for (size_t k = 0; k <= k_max; ++k)
			file_with_meander_numbers << mnd_calc.get_meander_number(n, k) << (k == k_max ? '\n' : ' ');
	}
	file_with_meander_numbers.close();

	return 0;
}