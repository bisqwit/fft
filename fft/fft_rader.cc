#include <mutex>
#include <algorithm>
#include <unordered_map>
#include "fft_rader.hh"
#include "fft_any.hh"
#include "factor.hh"

/* Compute n^m mod p, where m >= 0 and p > 0. */
static std::size_t powermod(std::size_t n, std::size_t m, std::size_t p)
{
    if(m == 0) return 1;
    if(m % 2 == 0) { n = powermod(n, m/2, p); return (n*n) % p; }
    return (n * powermod(n, m-1, p)) % p;
}

static std::pair<std::size_t, std::size_t> find_generator(std::size_t prime)
{
    // log2 of the product of first 16 primes ≈ 64.82. If size_t ≤ 64 bits,
    const std::size_t maxfac = 16;  // then size 16 array is enough.
    if(prime == 2) return { 1, 1 }; // prime = N , pm1 = N - 1.
    std::size_t factors[maxfac], pm1 = prime-1, g = 2; // Smallest possible g = 2.
    std::size_t count = get_factors(pm1, factors, maxfac);
    for(std::size_t m=0; m<count; )
        if(powermod(g, pm1 / factors[m], prime) == 1)
            { ++g; m = 0; }
        else
            { ++m; }
    return { g, powermod(g, prime-2, prime) };
}

struct RaderConfiguration
{
    cvector omega;
    std::size_t g;
    std::size_t ginv;
};
#include <assert.h>
static std::unordered_map<std::size_t, RaderConfiguration> data;
static std::mutex lock;
static const RaderConfiguration& configure_rader(std::size_t N)
{
    std::unique_lock<std::mutex> lk(lock);
    if(auto i = data.find(N); i != data.end()) return i->second;
    lk.unlock();

    auto [g, ginv] = find_generator(N); // (g*ginv) % N  is 1.
    assert((g*ginv)%N==1);

    exp_vector exps = ExpCircle(N);

    cvector omega(N-1);
    for(std::size_t gpower=1, i=0; i < N-1; ++i, gpower = (gpower*ginv)%N)
    {
        omega[i] = exps(gpower,N);
    }

    RaderConfiguration cfg;
    cfg.g     = g;
    cfg.ginv  = ginv;
    cfg.omega = FFT_any_fast(omega);
    for(auto& f: cfg.omega) f /= float(N-1);
    lk.lock();
    return data.emplace(N, std::move(cfg)).first->second;
}

void FFT_rader(const complex* input, complex* output, std::size_t N,
               std::size_t instride, std::size_t outstride)
{
    //fprintf(stderr, "N=%zu, ", N); fflush(stderr);
    const auto& cfg = configure_rader(N);
    const std::size_t g = cfg.g, ginv = cfg.ginv;
    //printf("g=%zu, ginv=%zu, instride=%zu, outstride=%zu\n", g,ginv,instride,outstride);

    cvector buf(N-1);

    // Permutation of input
    const complex in0 = input[0];
    for(std::size_t gpower=1, k=0; k<N-1; ++k, gpower = (gpower*g)%N)
        buf[k] = input[gpower*instride];

    // FFT
    exp_vector exps = ExpCircle(N-1);
    FFT_any_fast(&buf[0], output + outstride, exps, N-1, 1, outstride);

    // Set output DC component:
    output[0] = in0 + output[outstride];

    // Multiply by omega:
    if(outstride == 1)
        for(std::size_t k=0; k<N-1; ++k)
            output[k+1] = std::conj( output[k+1] * cfg.omega[k] + (k==0 ? in0 : 0));
    else
        for(std::size_t k=0; k<N-1; ++k)
            output[(k+1)*outstride] = std::conj( output[(k+1)*outstride] * cfg.omega[k] + (k==0 ? in0 : 0));
    // this will add input[0] to all of the outputs after the ifft
    //output[outstride] += std::conj(in0);

    // IFFT
    FFT_any_fast(output+outstride, &buf[0], exps, N-1, outstride, 1);

    // Inverse permutation
    for(std::size_t gpower=1, k=0; k<N-1; ++k, gpower = (gpower*ginv)%N)
        output[gpower*outstride] = std::conj(buf[k]);
}

cvector FFT_rader(const cvector& input)
{
    const auto N = input.size();
    cvector result(N);
    FFT_rader(&input[0], &result[0], N, 1, 1);
    return result;
}
