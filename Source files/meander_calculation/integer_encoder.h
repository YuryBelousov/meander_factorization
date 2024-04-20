#pragma once
#include <utility>

/*
This file contains auxiliary technical functions that are used to encode a pair of integers into 
a single integer. We use bit shifts for this purpose.
*/

using int_pair = std::pair<int, int>;
using int_pair_encoded = int;

namespace int_encoder {
	constexpr int shift = 7;

	constexpr int_pair_encoded encode(int x, int y) noexcept{
		return (x << shift) | y;
	}

	constexpr int_pair_encoded encode_2(const int_pair_encoded& code_x, const int_pair_encoded& code_y) noexcept {
		return ((code_x << shift) << shift) | code_y;
	}

	constexpr int decode_first_coord(const int_pair_encoded& code) noexcept {
		return code >> shift;
	}

	constexpr int decode_second_coord(const int_pair_encoded& code) noexcept {
		return code & ((1 << shift) - 1);
	}

	constexpr int_pair_encoded decode_2_first_coord(const int_pair_encoded& code) noexcept {
		return (code >> shift) >> shift;
	}

	constexpr int_pair_encoded decode_2_second_coord(const int_pair_encoded& code) noexcept {
		return code & (((1 << shift) << shift) - 1);
	}

	constexpr int_pair decode(const int_pair_encoded& code) noexcept {
		return int_pair(decode_first_coord(code), decode_second_coord(code));
	}

	constexpr int_pair decode_2(const int_pair_encoded& coded_code) noexcept {
		return int_pair(decode_2_first_coord(coded_code), decode_2_second_coord(coded_code));
	}

	constexpr int_pair_encoded code_of_difference(const int_pair_encoded& code_1, const int_pair_encoded& code_2) noexcept {
		auto x = decode_first_coord(code_1) - decode_first_coord(code_2);
		auto y = decode_second_coord(code_1) - decode_second_coord(code_2);
		return encode(x, y);
	}
};
