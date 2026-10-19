#pragma once
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <stdexcept>

#include "boost/multiprecision/cpp_int.hpp"

namespace irreducible_meander_enumeration {

#if defined(__SIZEOF_INT128__)
using count_type = unsigned __int128;
#else
using count_type = boost::multiprecision::uint128_t;
#endif
using key_type = std::uint64_t;
using word_type = std::uint64_t;
using count_table = std::array<std::array<count_type, 28>, 55>;

inline std::uint64_t hash(key_type key) {
	key ^= key >> 30;
	key *= 0xbf58476d1ce4e5b9ULL;
	key ^= key >> 27;
	key *= 0x94d049bb133111ebULL;
	return key ^ (key >> 31);
}

// Frontier ends form a balanced bracket word (opening = 1), followed by a 1 sentinel.
// The key stores this word in 53 bits, the river split in 6 bits, and k in 5 bits.
inline word_type word(key_type key) { return key & ((word_type(1) << 53) - 1); }
inline int split(key_type key) { return int((key >> 53) & 63); }
inline int tangencies(key_type key) { return int(key >> 59); }
inline int length(word_type bits) { return 63 - std::countl_zero(bits); }
inline int symbol(word_type bits, int i) { return int((bits >> i) & 1); }

inline word_type reverse_bits(word_type bits) {
	bits = ((bits >> 1) & 0x5555555555555555ULL) | ((bits & 0x5555555555555555ULL) << 1);
	bits = ((bits >> 2) & 0x3333333333333333ULL) | ((bits & 0x3333333333333333ULL) << 2);
	bits = ((bits >> 4) & 0x0f0f0f0f0f0f0f0fULL) | ((bits & 0x0f0f0f0f0f0f0f0fULL) << 4);
	bits = ((bits >> 8) & 0x00ff00ff00ff00ffULL) | ((bits & 0x00ff00ff00ff00ffULL) << 8);
	bits = ((bits >> 16) & 0x0000ffff0000ffffULL) | ((bits & 0x0000ffff0000ffffULL) << 16);
	return (bits >> 32) | (bits << 32);
}

inline key_type pack_key(word_type bits, int h, int s, int k, bool reflection) {
	key_type key = bits | (word_type(s) << 53);
	if (reflection && h) {
		const word_type sentinel = word_type(1) << h;
		const word_type mask = sentinel - 1;
		const word_type reflected = ((reverse_bits(bits & mask) >> (64 - h)) ^ mask) | sentinel;
		key = std::min(key, reflected | (word_type(h - s) << 53));
	}
	return key | (word_type(k) << 59);
}

inline word_type replace(word_type bits, int i, int value) {
	const word_type mask = word_type(1) << i;
	return (bits & ~mask) | (word_type(value) << i);
}

inline word_type insert_pair(word_type bits, int i) {
	const word_type mask = (word_type(1) << i) - 1;
	return (bits & mask) | ((bits & ~mask) << 2) | (word_type(1) << i);
}

inline word_type erase_pair(word_type bits, int i) {
	const word_type mask = (word_type(1) << i) - 1;
	return (bits & mask) | ((bits >> (i + 2)) << i);
}

inline int partner(word_type bits, int i, int h) {
	const int direction = symbol(bits, i) ? 1 : -1;
	int level = 1;
	for (int j = i + direction; j >= 0 && j < h; j += direction) {
		level += direction * (2 * symbol(bits, j) - 1);
		if (!level) return j;
	}
	throw std::logic_error("Invalid frontier brackets");
}

// A matched pair would close a separate component. Only the final pair may close,
// which is handled by the enumeration; otherwise splice the two remaining mates.
inline bool join(word_type& bits, int i, int h) {
	const int a = symbol(bits, i), b = symbol(bits, i + 1);
	if (a == 1 && b == 0) return false;
	if (a == b) bits = replace(bits, partner(bits, a ? i + 1 : i, h), a);
	bits = erase_pair(bits, i);
	return true;
}

} // namespace irreducible_meander_enumeration
