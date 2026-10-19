# Meander factorization

C++20 companion code for [Prime Factorization of Meanders](https://arxiv.org/abs/2112.10289). It calculates:

| Number | Counts |
|---|---|
| $I_{n,k}$ | Irreducible meanders |
| $S_{n,k}$ | Iterated snakes |
| $M_{n,k}$ | All singular meanders |

Here, $n$ is the number of transverse crossings and $k$ is the number of tangencies. The repository includes irreducible counts for $n+2k\le51$.

[ALGORITHM.md](ALGORITHM.md) describes the enumeration and generating-function calculations.

## Build and run

Use a C++20 compiler on a 64-bit machine with Boost.Multiprecision. With GCC, adjusting the Boost include path:

```sh
g++ -std=c++20 -O3 -DNDEBUG -pthread -I"Source files" -I/path/to/boost "Source files/main.cpp" -o meanders
./meanders
```

Run from the repository root. The example prints elapsed times and saves results in `Files with numbers`. Calculation limits are set in `Source files/main.cpp`.

## Irreducible enumeration

Include `irreducible_meander_enumeration/enumerate_irreducible_meanders.hpp`:

```cpp
const auto numbers = enumerate_irreducible_meanders(32);
std::cout << numbers[6][4] << std::endl;
```

The result contains exact `boost::multiprecision::cpp_int` values indexed by `[n][k]`, complete for `n + 2k <= N`. Entries outside the cutoff are left zero.

The full call is `enumerate_irreducible_meanders(N, memory_limit, threads, directory)`.

- `N`: Inclusive cutoff for `n + 2k`.
- `memory_limit`: Memory budget in bytes. Default: 32 GiB. Process overhead and the operating system's file cache are additional.
- `threads`: Number of worker threads, from 1 to 256. Default: the hardware thread count, limited to this range.
- `directory`: Directory for temporary working files. Default: `meander_scratch`. Files belong to the current call and are removed when it ends. For large calculations, use a local SSD with sufficient free space.

## Numerical data and calculations

[Files with numbers/irreducible_meanders_51.txt](Files%20with%20numbers/irreducible_meanders_51.txt) contains the complete irreducible table through $n+2k\le51$. Its header gives the rectangular limits `45 25`; subsequent rows contain values for $n=1,\ldots,45$ and $k=3,\ldots,25$. Zeros outside $n+2k\le51$ are padding and do not represent calculated coefficients.

Use `meander_calculator<number_type>` to calculate meander numbers from this table. `read_irreducible_meander_numbers` loads the data, `calculate_all(n, k)` calculates numbers up to the given order, and the getters return individual values. General meander counts beyond the supplied irreducible data are incomplete. The example uses `boost::multiprecision::uint1024_t` and loads the table through 51.

## Measured performance

The table below shows measurements for a full irreducible enumeration through $n+2k\le51$ on an AMD Ryzen 9 9950X. It was built with GCC 16.1.0 and `-O3 -march=native -DNDEBUG`, using automatic partitioning, 32 threads and a 56 GiB memory budget.

| Elapsed time | Peak managed memory | Peak scratch space |
|---|---|---|
| 6 h 39 min 19 s | 55.5 GiB | 951.5 GiB |

The output from this run is in [Files with numbers/irreducible_meanders_51.txt](Files%20with%20numbers/irreducible_meanders_51.txt).

## Reading the code

Start with [Source files/main.cpp](Source%20files/main.cpp).

| File | Purpose |
|---|---|
| [enumerate_irreducible_meanders.hpp](Source%20files/irreducible_meander_enumeration/enumerate_irreducible_meanders.hpp) | Enumeration and conversion of closed counts to open counts |
| [frontier.hpp](Source%20files/irreducible_meander_enumeration/frontier.hpp) | Frontier states, connectivity, and reflection |
| [frontier_storage.hpp](Source%20files/irreducible_meander_enumeration/frontier_storage.hpp) | Shared state tables and temporary disk storage |
| [meander_series.hpp](Source%20files/irreducible_meander_enumeration/meander_series.hpp) | Exact inversion, fixed rows, and boundary coefficients |
| [meander_calculator.hpp](Source%20files/meander_calculation/meander_calculator.hpp) | Generating-function calculations for meanders and iterated snakes |
