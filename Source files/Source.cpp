#include <iostream>
#include <fstream>
#include <string>
#include <chrono>

#include "meander_calculation/meander_calculator.h"
#include "irreducible_meander_enumeration/brute_force_functions.h"


// Technical function. It returns a string containing the time passed.
std::string get_time_passed(const std::chrono::steady_clock::time_point& start) noexcept{
	auto end = std::chrono::steady_clock::now();
	auto duration_in_millisec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	double time = duration_in_millisec.count()/1000.;
	std::stringstream stream;
	stream << std::fixed << std::setprecision(2);
	std::string time_type(" sec.");
	if (time > 60) {
		time /= 60.;
		time_type = " min.";
		if (time > 60) {
			time /= 60.;
			time_type = " hours.";
			if (time > 24) {
				time /= 24.;
				time_type = " days.";
			}
		}
	}
	stream << time;
	return std::string(stream.str() + time_type);
}

int main() {
	auto start = std::chrono::steady_clock::now();

	//Calculate irreducible meander numbers up to n+2k <= max_n, see submeander.h
	brute_forse();
	std::cout << get_time_passed(start) << std::endl;
	std::ofstream file_with_irr_mnd_out("files with numbers/irreducible_meanders.txt");
	for (int n = 1; n < max_n - 5; ++n) {
		for (int k = 3; k < max_n / 2; ++k)
			file_with_irr_mnd_out << irreducible_meander_numbers[n - 1][k - 3] << " ";
		file_with_irr_mnd_out << std::endl;
	}
	file_with_irr_mnd_out.close();

	
	meander_calculator mnd;

	//Calculate non-singular iterated snakes numbers up to total order 100
	auto one_step = std::chrono::steady_clock::now();

	std::ofstream file_with_non_singular_iter_snakes("files with numbers/iterated_snakes_non_singular.txt");
	for (int n = 0; n <= 100; ++n) {
		std::cout << n << ") ";
		file_with_non_singular_iter_snakes << mnd.get_iterated_snakes_number(n, 0) << std::endl;
		std::cout << get_time_passed(one_step) << std::endl;
		one_step = std::chrono::steady_clock::now();
	}	
	file_with_non_singular_iter_snakes.close();

	//Calculate iterated snakes numbers up to n <= 30, k <= 15
	std::ofstream file_with_iter_snakes("files with numbers/iterated_snakes.txt");
	for (int n = 0; n <= 30; ++n) {
		std::cout << n << ") ";
		for (int k = 0; k <= 15; ++k) {			
			file_with_iter_snakes << mnd.get_iterated_snakes_number(n, k) << " ";
			std::cout << get_time_passed(one_step) << " ";
			one_step = std::chrono::steady_clock::now();
		}
		std::cout << std::endl;
		file_with_iter_snakes << std::endl;
	}
	file_with_iter_snakes.close();

	return 0;
}