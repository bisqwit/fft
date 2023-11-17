#include <cstddef>
#include <initializer_list>

/* Finds the smallest factor ≥ 1 of the given integer. If N is prime, returns N. */
std::size_t smallest_factor(std::size_t N);

/* Attempts to return a product of factors of N as close to sqrt(N) as possible. */
std::size_t half_factors(std::size_t N);

/* Builds a list of unique factors of N in target[], up to limit. Returns number of primes written. */
std::size_t get_factors(std::size_t N, std::size_t* target, std::size_t limit);

/* Returns true if N can be represented as a product of integers 2, 3 and 5, and 7. */
bool factors_only(std::size_t N, const std::size_t* tab, std::size_t nfac);
std::size_t small_factors_only(std::size_t N, std::size_t upto=7);

std::size_t biggest_factor(std::size_t N);
