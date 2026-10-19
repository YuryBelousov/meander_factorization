#pragma once
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

#include "frontier.hpp"

namespace irreducible_meander_enumeration {

namespace fs = std::filesystem;
constexpr unsigned shard_count = 64;
constexpr std::size_t buffer_bytes = 65536;
constexpr unsigned batch_size = 64;
constexpr std::uint64_t job_records = 16384;

struct budget_exceeded : std::runtime_error {
	budget_exceeded() : std::runtime_error("destination memory budget exceeded") {}
};

// Disk records contain one or two 64-bit limbs, never a numeric object's layout.
inline count_type read_count(const char* data, unsigned width) {
	std::uint64_t low;
	std::memcpy(&low, data, 8);
	if (width == 8) return low;
	std::uint64_t high;
	std::memcpy(&high, data + 8, 8);
	return (count_type(high) << 64) | low;
}

inline void write_count(char* data, count_type value, unsigned width) {
	const std::uint64_t low = std::uint64_t(value);
	std::memcpy(data, &low, 8);
	if (width == 16) {
		const std::uint64_t high = std::uint64_t(value >> 64);
		std::memcpy(data + 8, &high, 8);
	}
}

// The limit includes old arrays during rehash/promotion and explicit worker
// buffers. Small C++ control objects, stream buffers and OS file cache are extra.
class memory_budget {
	std::size_t limit;
	std::atomic_size_t used_bytes{0}, peak_bytes{0};
public:
	explicit memory_budget(std::size_t limit) : limit(limit) {}
	void acquire(std::size_t bytes) {
		auto old = used_bytes.load();
		do {
			if (bytes > limit || old > limit - bytes) throw budget_exceeded();
		} while (!used_bytes.compare_exchange_weak(old, old + bytes));
		auto peak = peak_bytes.load();
		while (peak < old + bytes && !peak_bytes.compare_exchange_weak(peak, old + bytes)) {}
	}
	void release(std::size_t bytes) { used_bytes.fetch_sub(bytes); }
	std::size_t used() const { return used_bytes.load(); }
	std::size_t peak() const { return peak_bytes.load(); }
};

struct memory_reservation {
	memory_budget& budget;
	std::size_t bytes;
	memory_reservation(memory_budget& b, std::size_t n) : budget(b), bytes(n) { budget.acquire(n); }
	~memory_reservation() { budget.release(bytes); }
};

inline count_type checked_add(count_type a, count_type b) {
	const count_type sum = a + b;
	if (sum < a) throw std::overflow_error("128-bit frontier count overflow");
	return sum;
}

// Separate key/count arrays use 16 bytes per slot until a carry requires
// a high limb. Promotion is exact and applies to one shard at a time.
class compact_map {
	memory_budget& budget;
	std::mutex& growth_mutex;
	std::unique_ptr<std::uint64_t[]> keys, lows, highs;
	std::size_t capacity = 0, entry_count = 0;

	void grow() {
		std::lock_guard lock(growth_mutex); // only one shard rehashes at a time
		const std::size_t next = capacity ? 2 * capacity : 16;
		const auto bytes = next * (highs ? 24 : 16);
		budget.acquire(bytes);
		std::unique_ptr<std::uint64_t[]> new_keys, new_lows, new_highs;
		try {
			new_keys.reset(new std::uint64_t[next]{});
			new_lows.reset(new std::uint64_t[next]);
			if (highs) new_highs.reset(new std::uint64_t[next]{});
		}
		catch (...) { budget.release(bytes); throw; }
		for (std::size_t i = 0; i < capacity; ++i) if (keys[i]) {
			std::size_t j = hash(keys[i]) & (next - 1);
			while (new_keys[j]) j = (j + 1) & (next - 1);
			new_keys[j] = keys[i]; new_lows[j] = lows[i];
			if (new_highs) new_highs[j] = highs[i];
		}
		const auto old_bytes = memory_bytes();
		keys = std::move(new_keys); lows = std::move(new_lows); highs = std::move(new_highs);
		capacity = next;
		budget.release(old_bytes);
	}
	void promote() {
		if (highs) return;
		std::lock_guard lock(growth_mutex);
		const auto bytes = capacity * 8;
		budget.acquire(bytes);
		try { highs.reset(new std::uint64_t[capacity]{}); }
		catch (...) { budget.release(bytes); throw; }
	}
public:
	compact_map(memory_budget& budget, std::mutex& growth) : budget(budget), growth_mutex(growth) {}
	~compact_map() { budget.release(memory_bytes()); }
	compact_map(const compact_map&) = delete;
	compact_map& operator=(const compact_map&) = delete;
	std::size_t size() const { return entry_count; }
	std::size_t memory_bytes() const { return capacity * (highs ? 24 : 16); }
	unsigned count_bytes() const { return highs ? 16 : 8; }
	void add(key_type key, std::uint64_t hash, count_type value) {
		if (!key || !value) throw std::logic_error("zero key/count in compact map");
		if (!capacity) grow();
		auto index = std::size_t(hash) & (capacity - 1);
		while (keys[index]) {
			if (keys[index] == key) {
				const count_type old = count_type(lows[index]) | (highs ? count_type(highs[index]) << 64 : 0);
				const auto sum = checked_add(old, value);
				if (sum >> 64) promote();
				lows[index] = std::uint64_t(sum);
				if (highs) highs[index] = std::uint64_t(sum >> 64);
				return;
			}
			index = (index + 1) & (capacity - 1);
		}
		if ((entry_count + 1) * 4 > capacity * 3) {
			grow();
			index = std::size_t(hash) & (capacity - 1);
			while (keys[index]) index = (index + 1) & (capacity - 1);
		}
		if (value >> 64) promote();
		keys[index] = key; lows[index] = std::uint64_t(value);
		if (highs) highs[index] = std::uint64_t(value >> 64);
		++entry_count;
	}
	template <class visitor_type> void for_each(visitor_type&& visit) const {
		for (std::size_t i = 0; i < capacity; ++i) if (keys[i])
			visit(keys[i], count_type(lows[i]) | (highs ? count_type(highs[i]) << 64 : 0));
	}
};

struct enumeration_options {
	fs::path scratch = "disk_scratch";
	std::size_t ram_bytes = std::size_t(32) << 30;
	unsigned partitions = 0; // zero selects adaptive hash buckets
	unsigned threads = std::clamp(std::thread::hardware_concurrency(), 1u, 256u);
	bool reflection = true, verbose = false;
};

struct layer_statistics {
	int vertex = 0;
	std::uint64_t source_states = 0, source_bytes = 0;
	std::uint64_t states = 0, bytes = 0, transitions = 0, generated = 0, scanned = 0;
	std::uint64_t narrow_states = 0, wide_states = 0;
	std::uint64_t bucket_records = 0, read_records = 0, written_records = 0;
	std::uint64_t read_bytes = 0, written_bytes = 0, scratch_peak_bytes = 0;
	std::size_t memory_peak = 0;
	unsigned workers = 0, partitions = 1, retries = 0;
	double seconds = 0, generation_seconds = 0, reduction_seconds = 0;
};

// Only generated files in a unique child directory are ever removed. An
// existing scratch directory and all its other contents remain untouched.
class scratch_files {
	std::vector<fs::path> owned_paths;
	fs::path run_directory;
public:
	explicit scratch_files(const fs::path& parent) {
		fs::create_directories(parent);
		const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
		for (unsigned suffix = 0;; ++suffix) {
			run_directory = parent / ("meanders-" + std::to_string(stamp) + "-" + std::to_string(suffix));
			if (fs::create_directory(run_directory)) break;
		}
	}
	~scratch_files() {
		std::error_code ec;
		for (const auto& path : owned_paths) fs::remove(path, ec);
		fs::remove(run_directory, ec); // removes only if empty, never recursive
	}
	fs::path file(int vertex, unsigned partition, unsigned count_bytes) {
		auto path = run_directory / ("v" + std::to_string(vertex) + "-p" + std::to_string(partition)
						   + "-c" + std::to_string(count_bytes) + ".bin");
		owned_paths.push_back(path);
		return path;
	}
	void remove(const fs::path& path) { fs::remove(path); }
};

struct source_file {
	fs::path path;
	std::uint64_t records = 0;
	unsigned count_bytes = 8;
	std::uint64_t bytes() const { return records * (8 + count_bytes); }
};
struct job { std::size_t file; std::uint64_t first, count; };

template <class visitor_type>
void read_records(const source_file& file, const job& current_job, std::vector<char>& buffer, visitor_type&& visit) {
	if (file.count_bytes != 8 && file.count_bytes != 16) throw std::logic_error("bad counter width");
	const unsigned record_bytes = 8 + file.count_bytes;
	std::ifstream input(file.path, std::ios::binary);
	if (!input) throw std::runtime_error("cannot open a source layer");
	input.seekg(std::streamoff(current_job.first * record_bytes));
	if (!input) throw std::runtime_error("cannot seek in a source layer");
	for (auto left = current_job.count; left;) {
		const auto count = std::min<std::uint64_t>(left, buffer.size() / record_bytes);
		if (!count) throw std::logic_error("record buffer too small");
		input.read(buffer.data(), std::streamsize(count * record_bytes));
		if (!input) throw std::runtime_error("source layer truncated or read failed");
		for (std::size_t i = 0; i < count; ++i) {
			key_type key;
			const char* p = buffer.data() + i * record_bytes;
			std::memcpy(&key, p, 8);
			const count_type ways = read_count(p + 8, file.count_bytes);
			if (!key || !ways) throw std::runtime_error("invalid zero record in source layer");
			visit(key, ways);
		}
		left -= count;
	}
}

struct pending_entry { key_type key; std::uint64_t hash; count_type ways; };

struct layer_result : layer_statistics {
	count_table counts{};
	std::vector<source_file> files;
};

struct io_accounting {
	std::uint64_t read_records = 0, written_records = 0, read_bytes = 0, written_bytes = 0;
	std::uint64_t live = 0, peak = 0;
	explicit io_accounting(std::uint64_t source_bytes) : live(source_bytes), peak(source_bytes) {}
	void read(std::uint64_t records, unsigned width) {
		read_records += records; read_bytes += records * (8 + width);
	}
	void write(std::uint64_t records, unsigned width) {
		written_records += records; const auto bytes = records * (8 + width);
		written_bytes += bytes; live += bytes; peak = std::max(peak, live);
	}
	void remove(const source_file& file, scratch_files& scratch) {
		scratch.remove(file.path); live -= file.bytes();
	}
};

struct file_factory {
	scratch_files& scratch;
	int vertex;
	unsigned next = 0;
	fs::path make(unsigned width) { return scratch.file(vertex, next++, width); }
};

struct bucket_file { source_file file; unsigned depth = 0; };

// Only buffer_bytes of explicit output buffers is live, regardless of the
// number of buckets. The C streams are unbuffered, so there is no hidden 64 KiB
// stream buffer per bucket. A 64-slot LRU also bounds open file handles.
class bucket_writer {
	struct buffer_slot {
		unsigned bucket = 0;
		std::FILE* file = nullptr;
		std::size_t used = 0;
		std::uint64_t touched = 0;
	};
	std::vector<source_file> files;
	std::array<buffer_slot, shard_count> slots{};
	std::vector<char> buffer = std::vector<char>(buffer_bytes);
	unsigned slot_count, partition_count, depth;
	std::size_t slot_bytes;
	std::uint64_t access_clock = 0;
	io_accounting& accounting;
	std::mutex mutex;
	bool finished = false;

	void flush(unsigned i) {
		auto& slot = slots[i];
		if (!slot.used) return;
		if (std::fwrite(buffer.data() + i * slot_bytes, 1, slot.used, slot.file) != slot.used)
			throw std::runtime_error("bucket write failed (disk full or I/O error)");
		accounting.write(slot.used / 24, 16); slot.used = 0;
	}
	void close(unsigned i) {
		auto& slot = slots[i];
		if (!slot.file) return;
		flush(i);
		auto* file = slot.file; slot.file = nullptr;
		if (std::fclose(file)) throw std::runtime_error("bucket close failed (disk full or I/O error)");
	}
	unsigned obtain(unsigned bucket) {
		unsigned candidate = 0;
		for (unsigned i = 0; i < slot_count; ++i) {
			if (slots[i].file && slots[i].bucket == bucket) return i;
			if (!slots[i].file) { candidate = i; break; }
			if (slots[i].touched < slots[candidate].touched) candidate = i;
		}
		close(candidate);
		auto& slot = slots[candidate];
		const auto& path = files[bucket].path;
#if defined(_MSC_VER)
		_wfopen_s(&slot.file, path.c_str(), L"ab");
#elif defined(_WIN32)
		slot.file = _wfopen(path.c_str(), L"ab");
#else
		slot.file = std::fopen(path.c_str(), "ab");
#endif
		if (!slot.file) throw std::runtime_error("cannot open destination bucket");
		if (std::setvbuf(slot.file, nullptr, _IONBF, 0)) {
			std::fclose(slot.file); slot.file = nullptr;
			throw std::runtime_error("cannot disable destination bucket buffering");
		}
		slot.bucket = bucket;
		return candidate;
	}
	void append_unlocked(key_type key, count_type ways, std::uint64_t hash) {
		if (!key || !ways) throw std::logic_error("zero key/count in bucket");
		const unsigned bucket = (hash >> (32 + depth)) & (partition_count - 1);
		const unsigned i = obtain(bucket);
		auto& slot = slots[i]; slot.touched = ++access_clock;
		if (slot.used + 24 > slot_bytes) flush(i);
		char* p = buffer.data() + i * slot_bytes + slot.used;
		std::memcpy(p, &key, 8); write_count(p + 8, ways, 16);
		slot.used += 24; ++files[bucket].records;
	}
public:
	bucket_writer(file_factory& factory, unsigned partitions, unsigned depth, io_accounting& accounting)
		: slot_count(std::min(partitions, shard_count)), partition_count(partitions), depth(depth),
		  slot_bytes(buffer_bytes / slot_count / 24 * 24), accounting(accounting) {
		files.reserve(partitions);
		for (unsigned i = 0; i < partitions; ++i) files.push_back({factory.make(16), 0, 16});
	}
	~bucket_writer() {
		for (auto& slot : slots) if (slot.file) std::fclose(slot.file);
	}
	void append(key_type key, count_type ways) {
		std::lock_guard lock(mutex); append_unlocked(key, ways, hash(key));
	}
	void append_batch(const pending_entry* items, std::size_t count) {
		std::lock_guard lock(mutex);
		for (std::size_t i = 0; i < count; ++i) append_unlocked(items[i].key, items[i].ways, items[i].hash);
	}
	std::vector<bucket_file> finish() {
		std::lock_guard lock(mutex);
		if (finished) throw std::logic_error("bucket writer finished twice");
		for (unsigned i = 0; i < slot_count; ++i) close(i);
		std::vector<bucket_file> result;
		const unsigned child_depth = depth + std::countr_zero(partition_count);
		for (const auto& file : files) if (file.records) {
			if (fs::file_size(file.path) != file.bytes()) throw std::runtime_error("bucket size mismatch");
			result.push_back({file, child_depth});
		}
		finished = true;
		return result;
	}
};

inline void write_tables(const std::array<std::unique_ptr<compact_map>, shard_count>& tables,
	file_factory& factory, io_accounting& accounting, layer_result& result) {
	std::vector<char> buffer(buffer_bytes);
	for (unsigned width : {8u, 16u}) {
		std::uint64_t records = 0;
		for (const auto& table : tables) if (table->count_bytes() == width) records += table->size();
		if (!records) continue;
		const auto path = factory.make(width);
		std::ofstream out(path, std::ios::binary | std::ios::trunc);
		if (!out) throw std::runtime_error("cannot create destination layer");
		const unsigned record_bytes = 8 + width;
		std::size_t used = 0;
		auto flush = [&] {
			if (!used) return;
			out.write(buffer.data(), std::streamsize(used));
			if (!out) throw std::runtime_error("destination write failed (disk full or I/O error)");
			accounting.write(used / record_bytes, width); used = 0;
		};
		for (const auto& table : tables) if (table->count_bytes() == width)
			table->for_each([&](key_type key, count_type ways) {
				if (used + record_bytes > buffer.size()) flush();
				std::memcpy(buffer.data() + used, &key, 8);
				write_count(buffer.data() + used + 8, ways, width); used += record_bytes;
			});
		flush(); out.close();
		if (!out) throw std::runtime_error("destination close failed (disk full or I/O error)");
		if (fs::file_size(path) != records * record_bytes) throw std::runtime_error("destination size mismatch");
		result.files.push_back({path, records, width}); result.states += records; result.bytes += records * record_bytes;
		(width == 8 ? result.narrow_states : result.wide_states) += records;
	}
}

inline unsigned workers_for(const enumeration_options& options, std::size_t jobs, std::size_t worker_bytes) {
	return std::min<unsigned>(unsigned(std::min<std::size_t>(options.threads, jobs)),
		unsigned(std::max<std::size_t>(1, options.ram_bytes / (4 * worker_bytes))));
}

// The only source-extension pass. In-memory batches already committed before a
// failing insertion remain in the tables; the failing item and later items are
// appended once after the tables have been spilled with every shard frozen.
// Ordinary batches only lock their own shard: no global lock is contended in
// the common case where the entire destination layer fits in memory.
template <class extension_type>
std::vector<bucket_file> generate(const std::vector<source_file>& source, const std::vector<job>& jobs,
	file_factory& factory, io_accounting& accounting, const enumeration_options& options,
	extension_type&& extend, layer_result& result, unsigned initial_partitions) {
	if (jobs.empty()) return {};
	const std::size_t worker_bytes = buffer_bytes + shard_count * batch_size * sizeof(pending_entry) + sizeof(count_table);
	const unsigned workers = workers_for(options, jobs.size(), worker_bytes);
	result.workers = workers;
	memory_budget budget(options.ram_bytes);
	memory_reservation buffers(budget, workers * worker_bytes + buffer_bytes);
	std::mutex growth;
	std::array<std::unique_ptr<compact_map>, shard_count> tables;
	std::array<std::mutex, shard_count> mutexes;
	for (auto& table : tables) table = std::make_unique<compact_map>(budget, growth);
	std::mutex spill_mutex;
	std::atomic_bool spilled{false};
	std::unique_ptr<bucket_writer> writer;
	const unsigned partitions = options.partitions ? options.partitions : initial_partitions;
	if (partitions > 1) {
		writer = std::make_unique<bucket_writer>(factory, partitions, 0, accounting);
		result.partitions = partitions;
		spilled.store(true, std::memory_order_release);
	}
	struct worker_result { count_table counts{}; std::uint64_t scanned = 0, generated = 0, read_bytes = 0; };
	std::vector<worker_result> local(workers);
	std::atomic_size_t next_job{0};
	std::atomic_bool cancelled{false};
	std::exception_ptr failure;
	std::mutex error_mutex;
	{
		std::vector<std::jthread> threads;
		for (unsigned worker = 0; worker < workers; ++worker) threads.emplace_back([&, worker] {
			try {
				std::vector<char> buffer(buffer_bytes);
				std::array<std::vector<pending_entry>, shard_count> batches;
				for (auto& batch : batches) batch.reserve(batch_size);
				auto flush = [&](unsigned shard) {
					auto& batch = batches[shard];
					if (batch.empty()) return;
					std::size_t first = 0;
					bool spill = false;
					if (spilled.load(std::memory_order_acquire)) {
						writer->append_batch(batch.data(), batch.size());
					}
					else {
						std::lock_guard shard_lock(mutexes[shard]);
						// A spill may have frozen this shard while the worker
						// waited for its lock. Publication follows the complete
						// dump/reset, so the writer is ready when this is true.
						if (spilled.load(std::memory_order_acquire)) {
							writer->append_batch(batch.data(), batch.size());
						}
						else {
							try {
								for (; first < batch.size(); ++first)
									tables[shard]->add(batch[first].key, batch[first].hash, batch[first].ways);
							}
							catch (const budget_exceeded&) { spill = true; }
						}
					}
					if (spill) {
						if (options.partitions == 1)
							throw std::runtime_error("fixed destination partition count is insufficient for memory budget");
						std::lock_guard spill_lock(spill_mutex);
						if (!spilled.load(std::memory_order_acquire)) {
							// Every ordinary flush holds at most one shard
							// mutex. Taking all of them in ascending order here
							// cannot form a lock-order cycle, even if several
							// workers hit the budget at the same time.
							std::array<std::unique_lock<std::mutex>, shard_count> frozen;
							for (unsigned i = 0; i < shard_count; ++i)
								frozen[i] = std::unique_lock<std::mutex>(mutexes[i]);
							writer = std::make_unique<bucket_writer>(factory, 2, 0, accounting);
							for (const auto& table : tables)
								table->for_each([&](key_type key, count_type ways) { writer->append(key, ways); });
							for (auto& table : tables) table.reset();
							result.partitions = 2; ++result.retries;
							spilled.store(true, std::memory_order_release);
						}
						writer->append_batch(batch.data() + first, batch.size() - first);
					}
					batch.clear();
				};
				auto& own = local[worker];
				for (std::size_t i; !cancelled && (i = next_job++) < jobs.size();) {
					const auto& job = jobs[i];
					own.read_bytes += job.count * (8 + source[job.file].count_bytes);
					read_records(source[job.file], job, buffer, [&](key_type key, count_type ways) {
						++own.scanned;
						extend(key, [&](key_type child) {
							++own.generated;
							const auto hash_value = hash(child);
							const unsigned shard = (hash_value >> 48) & (shard_count - 1);
							batches[shard].push_back({child, hash_value, ways});
							if (batches[shard].size() == batch_size) flush(shard);
						}, [&](int n, int k) {
							if (n < 0 || std::size_t(n) >= own.counts.size() || k < 0
								|| std::size_t(k) >= own.counts[n].size())
								throw std::logic_error("completion outside count table");
							++own.generated;
							own.counts[n][k] = checked_add(own.counts[n][k], ways);
						});
					});
				}
				if (!cancelled) for (unsigned shard = 0; shard < shard_count; ++shard) flush(shard);
			}
			catch (...) {
				cancelled = true;
				std::lock_guard lock(error_mutex);
				if (!failure) failure = std::current_exception();
			}
		});
	}
	result.memory_peak = budget.peak();
	if (failure) std::rethrow_exception(failure);
	for (const auto& own : local) {
		for (std::size_t n = 0; n < result.counts.size(); ++n)
			for (std::size_t k = 0; k < result.counts[n].size(); ++k)
				result.counts[n][k] = checked_add(result.counts[n][k], own.counts[n][k]);
		result.scanned += own.scanned; result.generated += own.generated;
		accounting.read_records += own.scanned; accounting.read_bytes += own.read_bytes;
	}
	result.transitions = result.generated;
	if (writer) {
		auto buckets = writer->finish();
		for (const auto& bucket : buckets) result.bucket_records += bucket.file.records;
		return buckets;
	}
	write_tables(tables, factory, accounting, result);
	return {};
}

inline bool reduce(const bucket_file& bucket, file_factory& factory, io_accounting& accounting,
	const enumeration_options& options, layer_result& result) {
	std::vector<job> jobs;
	for (std::uint64_t first = 0; first < bucket.file.records; first += job_records)
		jobs.push_back({0, first, std::min(job_records, bucket.file.records - first)});
	const std::size_t worker_bytes = buffer_bytes + shard_count * batch_size * sizeof(pending_entry);
	const unsigned workers = workers_for(options, jobs.size(), worker_bytes);
	result.workers = std::max(result.workers, workers);
	memory_budget budget(options.ram_bytes);
	memory_reservation buffers(budget, workers * worker_bytes + buffer_bytes);
	std::mutex growth;
	std::array<std::unique_ptr<compact_map>, shard_count> tables;
	std::array<std::mutex, shard_count> mutexes;
	for (auto& table : tables) table = std::make_unique<compact_map>(budget, growth);
	std::vector<std::uint64_t> read(workers);
	std::atomic_size_t next_job{0};
	std::atomic_bool cancelled{false};
	std::exception_ptr failure;
	std::mutex error_mutex;
	{
		std::vector<std::jthread> threads;
		for (unsigned worker = 0; worker < workers; ++worker) threads.emplace_back([&, worker] {
			try {
				std::vector<char> buffer(buffer_bytes);
				std::array<std::vector<pending_entry>, shard_count> batches;
				for (auto& batch : batches) batch.reserve(batch_size);
				auto flush = [&](unsigned shard) {
					if (batches[shard].empty()) return;
					std::lock_guard lock(mutexes[shard]);
					for (const auto& item : batches[shard]) tables[shard]->add(item.key, item.hash, item.ways);
					batches[shard].clear();
				};
				for (std::size_t i; !cancelled && (i = next_job++) < jobs.size();) {
					// read_records reads a complete input buffer before visits;
					// count physical reads, including a buffer whose reduction
					// later fails partway through its visitor callbacks.
					const auto& job = jobs[i];
					std::uint64_t visited = 0;
					constexpr auto per_buffer = buffer_bytes / 24;
					read_records(bucket.file, job, buffer, [&](key_type key, count_type ways) {
						if (visited % per_buffer == 0)
							read[worker] += std::min<std::uint64_t>(per_buffer, job.count - visited);
						++visited;
						const auto hash_value = hash(key);
						const unsigned shard = (hash_value >> 48) & (shard_count - 1);
						batches[shard].push_back({key, hash_value, ways});
						if (batches[shard].size() == batch_size) flush(shard);
					});
				}
				if (!cancelled) for (unsigned shard = 0; shard < shard_count; ++shard) flush(shard);
			}
			catch (...) {
				cancelled = true;
				std::lock_guard lock(error_mutex);
				if (!failure) failure = std::current_exception();
			}
		});
	}
	for (auto count : read) accounting.read(count, 16);
	result.memory_peak = std::max(result.memory_peak, budget.peak());
	if (failure) {
		try { std::rethrow_exception(failure); }
		catch (const budget_exceeded&) { return false; }
	}
	write_tables(tables, factory, accounting, result);
	return true;
}

inline std::vector<bucket_file> split_bucket(const bucket_file& bucket, file_factory& factory,
	io_accounting& accounting, const enumeration_options& options,
	layer_result& result) {
	memory_budget budget(options.ram_bytes);
	memory_reservation buffers(budget, 2 * buffer_bytes);
	std::vector<char> buffer(buffer_bytes);
	bucket_writer writer(factory, 2, bucket.depth, accounting);
	read_records(bucket.file, {0, 0, bucket.file.records}, buffer,
		[&](key_type key, count_type ways) { writer.append(key, ways); });
	accounting.read(bucket.file.records, 16);
	result.memory_peak = std::max(result.memory_peak, budget.peak());
	return writer.finish();
}

template <class extension_type>
layer_result build_layer(const std::vector<source_file>& source, const std::vector<job>& jobs,
						scratch_files& scratch, int vertex, const enumeration_options& options, extension_type&& extend,
						unsigned initial_partitions = 1) {
	if (!options.threads || options.threads > 256 || options.ram_bytes < 1048576)
		throw std::invalid_argument("bucket transfer needs 1..256 threads and at least 1 MiB budget");
	if (options.partitions > 4096 || (options.partitions && !std::has_single_bit(options.partitions)))
		throw std::invalid_argument("partitions must be zero (auto) or a power of two up to 4096");
	if (!initial_partitions || initial_partitions > 4096 || !std::has_single_bit(initial_partitions))
		throw std::invalid_argument("initial bucket hint must be a power of two up to 4096");
	std::uint64_t source_bytes = 0;
	for (const auto& file : source) source_bytes += file.bytes();
	io_accounting accounting(source_bytes);
	file_factory factory{scratch, vertex};
	layer_result result;
	const auto start = std::chrono::steady_clock::now();
	auto pending = generate(source, jobs, factory, accounting, options, std::forward<extension_type>(extend), result,
							initial_partitions);
	const auto generated = std::chrono::steady_clock::now();
	result.generation_seconds = std::chrono::duration<double>(generated - start).count();
	const bool bucketized = !pending.empty();
	std::uint64_t initial_leaves = pending.size();
	while (!pending.empty()) {
		auto bucket = std::move(pending.back()); pending.pop_back();
		if (reduce(bucket, factory, accounting, options, result)) {
			accounting.remove(bucket.file, scratch);
		}
		else {
			if (options.partitions)
				throw std::runtime_error("fixed destination partition count is insufficient for memory budget");
			if (bucket.depth == 12)
				throw std::runtime_error("destination bucket still exceeds memory budget at 4096 hash partitions");
			auto children = split_bucket(bucket, factory, accounting, options, result);
			accounting.remove(bucket.file, scratch);
			initial_leaves += children.size() - 1;
			pending.insert(pending.end(), std::make_move_iterator(children.begin()), std::make_move_iterator(children.end()));
			++result.retries;
			result.partitions = std::max(result.partitions, 1u << (bucket.depth + 1));
		}
	}
	if (bucketized) result.partitions = std::max(result.partitions, unsigned(initial_leaves));
	result.reduction_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - generated).count();
	result.read_records = accounting.read_records; result.written_records = accounting.written_records;
	result.read_bytes = accounting.read_bytes; result.written_bytes = accounting.written_bytes;
	result.scratch_peak_bytes = accounting.peak;
	return result;
}

} // namespace irreducible_meander_enumeration
