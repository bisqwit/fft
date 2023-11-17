#include <set>
#include <cmath>
#include <mutex>
#include <vector>
#include <type_traits>
#include <unordered_map>
#include "factor.hh"
static class lore /* This structure caches prime-related information. */
{
    std::unordered_map<std::size_t, std::size_t> smallest_factors{ {0,0}, {1,1}, {2,2} };
    std::unordered_map<std::size_t, std::size_t> half_factors{};
    std::vector<signed char> is_prime{1,1,1}; // 0=not known, -1=no, 1=yes
    std::size_t              lastprime{2};    // extent of contiguous knowledge
    std::set<std::size_t>    primes{2};       // list of known primes
    std::mutex lock{};
public:
    auto get_lock() { return std::unique_lock<std::mutex>(lock); }
    template<typename T>
    std::size_t find(std::size_t N, T&&) requires std::is_same_v<std::remove_reference_t<T>,std::unique_lock<std::mutex>>
    {
        if(auto i = smallest_factors.find(N); i != smallest_factors.end()) return i->second;
        std::size_t solution = N;
        for(auto p: primes)
        {
            if(N % p == 0) { solution = p; break; }
            if(p*p >= N) { break; }
        }
        for(std::size_t p = lastprime|1; p*p <= N; p+=2) if(N % p == 0) { solution = p; break; }
        if(is_prime.size() <= N) is_prime.resize(N+1);
        for(is_prime[N] = (solution==N) ? 1 : -1; lastprime+1 < is_prime.size() && is_prime[lastprime+1] != 0; ) ++lastprime;
        if(solution == N) primes.insert(solution);
        return smallest_factors.emplace(N, solution).first->second;
    }
    std::size_t find_half(std::size_t N, std::unique_lock<std::mutex>&&)
    {
        if(auto i = half_factors.find(N); i != half_factors.end()) return i->second;
        std::size_t result = N, a = std::sqrt(N);
        for(a |= 1; a > 1; a -= (a>3?2:1)) { if(N % a == 0) { result = a; break; } }
        half_factors.emplace(N, result);
        return result;
    }
    lore() { auto l = get_lock(); for(std::size_t a=0; a<1048576*2; ++a) { find(a, l); if(a<31) find(1<<a, l); } }
} prime_lore;

std::size_t smallest_factor(std::size_t N)
{
    return prime_lore.find(N, prime_lore.get_lock());
}
std::size_t get_factors(std::size_t N, std::size_t* target, std::size_t limit)
{
    auto lock = prime_lore.get_lock();
    for(std::size_t orig=N,count=0,prev=0,f; ; N/=f)
        if(f=prime_lore.find(N, lock); count >= limit || f <= 1 || f == orig) { return count; }
        else if(f != prev) { target[count++] = prev = f; }
}
std::size_t small_factors_only(std::size_t N, std::size_t upto)
{
    auto lock = prime_lore.get_lock();
    for(size_t f; ; N/=f)
        { f = prime_lore.find(N, lock); if(f < 2 || f > upto) return f == 1; }
}

bool factors_only(std::size_t N, const std::size_t* tab, std::size_t nfac)
{
    for(std::size_t f=0; N && f<nfac; ++f) for(auto v = tab[f]; N % v == 0; N /= v) { }
    return N == 1;
}
std::size_t half_factors(std::size_t N)
{
    return prime_lore.find_half(N, prime_lore.get_lock());
}
std::size_t biggest_factor(std::size_t N)
{
    std::size_t array[16], n = get_factors(N, array, 16);
    return n ? array[n-1] : N;
}
