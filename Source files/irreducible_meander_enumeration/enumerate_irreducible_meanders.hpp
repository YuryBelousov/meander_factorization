#pragma once
#include <iostream>
#include "frontier_storage.hpp"
#include "meander_series.hpp"

namespace irreducible_meander_enumeration {

// Each step adds one contact: v = n+k. Thus k in the key also determines n.
template <class emit_function, class complete_function>
inline void extend_frontier(key_type key, int vertex, int max_weight, int max_transverse,
	bool reflection, emit_function&& emit, complete_function&& complete) {
	const int next_vertex = vertex + 1;
	auto add = [&](word_type bits, int h, int s, int k) {
		const int n = next_vertex - k, weight = next_vertex + k;
		if (n > max_transverse || weight > max_weight) return;
		// Closing costs at least max(s, h-s) in weight. Only crossings can save
		// weight over tangencies, and retained rows have at least n=6.
		int remaining = std::max({ s, h - s, h - (max_transverse - n), 6 - n });
		if (h > 2 && 2 * s == h && symbol(bits, s - 1) == 1 && symbol(bits, s) == 0)
			remaining = std::max(remaining, h / 2 + 2);
		if (remaining <= max_weight - weight)
			emit(pack_key(bits, h, s, k, reflection));
	};
	auto close_or_add = [&](word_type bits, int h, int s, int k, int index) {
		const int n = next_vertex - k;
		if (n > max_transverse || next_vertex + k > max_weight) return;
		if (h == 2) {
			if (n >= 6) complete(n, k);
		}
		else if (join(bits, index, h)) add(bits, h - 2, s, k);
	};
	const word_type bits = word(key);
	const int h = length(bits), s = split(key), k = tangencies(key);
	if (h + 2 <= 52) add(insert_pair(bits, s), h + 2, s + 1, k);
	if (s > 0) add(bits, h, s - 1, k);
	if (s < h) add(bits, h, s + 1, k);
	if (s > 0 && s < h) close_or_add(bits, h, s - 1, k, s - 1);
	if (h + 2 <= 52) {
		const word_type child = insert_pair(bits, s);
		add(child, h + 2, s, k + 1);
		add(child, h + 2, s + 2, k + 1);
	}
	// These two tangencies have the same frontier but are distinct configurations.
	if (s > 0) add(bits, h, s, k + 1);
	if (s < h) add(bits, h, s, k + 1);
	if (s >= 2) close_or_add(bits, h, s - 2, k + 1, s - 2);
	if (h - s >= 2) close_or_add(bits, h, s, k + 1, s);
}

inline meander_series enumerate_all_meanders(int max_weight, int max_transverse,
	const enumeration_options& options, std::vector<layer_statistics>& statistics) {
	if (max_weight < 2 || max_weight > 53 || max_transverse < 0 || max_transverse > std::min(max_weight, 52))
		throw std::invalid_argument("Enumeration requires 2 <= cutoff <= 53 and 0 <= transverse cutoff <= min(cutoff, 52)");
	meander_series result(max_weight, max_transverse);
	const auto low_rows = all_meander_rows(max_weight / 2);
	for (int n = 0; n <= std::min(4, max_transverse); ++n)
		for (int k = 0; result.contains(n, k); ++k)
			result(n, k) = low_rows[n][k];
	if (max_transverse <= 4) return result;

	const int closed_weight = max_weight + 2 - max_weight % 2;
	const int closed_transverse = max_transverse + max_transverse % 2;
	count_table closed{};
	scratch_files scratch(options.scratch);
	const auto initial = scratch.file(0, 0, 8);
	{
		std::ofstream out(initial, std::ios::binary);
		const key_type key = pack_key(1, 0, 0, 0, options.reflection);
		const std::uint64_t one = 1;
		out.write(reinterpret_cast<const char*>(&key), 8);
		out.write(reinterpret_cast<const char*>(&one), 8);
		out.close();
		if (!out) throw std::runtime_error("Cannot write initial frontier");
	}
	std::vector<source_file> source{ { initial, 1, 8 } };
	unsigned partition_hint = 1;
	for (int vertex = 1; vertex <= (closed_weight + closed_transverse) / 2 && !source.empty(); ++vertex) {
		const auto start = std::chrono::steady_clock::now();
		std::uint64_t source_states = 0, source_bytes = 0;
		std::vector<job> jobs;
		for (size_t i = 0; i < source.size(); ++i) {
			source_states += source[i].records;
			source_bytes += source[i].bytes();
			for (std::uint64_t first = 0; first < source[i].records; first += job_records)
				jobs.push_back({ i, first, std::min(job_records, source[i].records - first) });
		}
		auto next = build_layer(source, jobs, scratch, vertex, options,
			[&](key_type key, auto&& emit, auto&& complete) {
				extend_frontier(key, vertex - 1, closed_weight, closed_transverse,
					options.reflection, emit, complete);
			}, partition_hint);
		if (!options.partitions) {
			partition_hint = next.partitions;
			if (partition_hint > 1 && next.memory_peak < options.ram_bytes / 4) partition_hint /= 2;
		}
		for (int n = 6; n <= closed_transverse; n += 2)
			for (int k = 0; n + 2 * k <= closed_weight; ++k)
				closed[n][k] = checked_add(closed[n][k], next.counts[n][k]);
		next.vertex = vertex;
		next.source_states = source_states;
		next.source_bytes = source_bytes;
		for (const auto& file : source) scratch.remove(file.path);
		source = std::move(next.files);
		next.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
		statistics.push_back(next);
		if (options.verbose)
			std::cerr << "Contact " << vertex << ": " << next.states << " states, "
				<< double(next.memory_peak) / 1048576 << " MiB, " << next.seconds << " sec\n";
	}
	// Exact open/closed rank identities. Rows 0..4 have already been supplied by formulas.
	for (int n = 5; n <= max_transverse; ++n)
		for (int k = 0; result.contains(n, k); ++k) {
			const count_type count = n % 2 ? closed[n + 1][k] : closed[n][k + 1];
			integer value = (integer(std::uint64_t(count >> 64)) << 64) + std::uint64_t(count);
			const int numerator = n % 2 ? n + 1 : k + 1;
			const int denominator = n % 2 ? n + k + 1 : 2 * (n + k + 1);
			value *= numerator;
			if (value % denominator != 0) throw std::logic_error("Nonintegral open/closed conversion");
			result(n, k) = value / denominator;
		}
	return result;
}

} // namespace irreducible_meander_enumeration

// Calculate the complete table [n][k] for n+2k <= max_n (0..53).
// Cutoffs 52 and 53 share the same transfer; cyclic identities supply the odd boundary.
// memory_limit counts frontier arrays and worker buffers, in bytes. Temporary files
// are created in a unique subdirectory of directory and removed when the call ends.
inline std::vector<std::vector<boost::multiprecision::cpp_int>> enumerate_irreducible_meanders(
	int max_n, size_t memory_limit = size_t(32) << 30,
	unsigned number_of_threads = std::clamp(std::thread::hardware_concurrency(), 1u, 256u),
	const std::filesystem::path& directory = "meander_scratch") {
	using namespace irreducible_meander_enumeration;
	if (max_n < 0 || max_n > 53) throw std::invalid_argument("Enumeration cutoff must be between 0 and 53");
	if (!number_of_threads || number_of_threads > 256 || memory_limit < 1048576)
		throw std::invalid_argument("Enumeration requires 1..256 threads and at least 1 MiB of memory");
	const int transverse = max_n < 14 ? 0 : std::max(4, max_n - 8 - max_n % 2);
	meander_series irreducibles(max_n, transverse);
	if (max_n >= 14) {
		enumeration_options options;
		options.ram_bytes = memory_limit;
		options.threads = number_of_threads;
		options.scratch = directory;
		std::vector<layer_statistics> statistics;
		irreducibles = derive_irreducibles(enumerate_all_meanders(max_n, transverse, options, statistics));
	}
	const auto complete = complete_irreducibles(irreducibles);
	std::vector<std::vector<integer>> result(max_n + 1, std::vector<integer>(max_n / 2 + 1));
	for (int n = 0; n <= max_n; ++n)
		for (int k = 0; complete.contains(n, k); ++k)
			result[n][k] = complete(n, k);
	return result;
}
