#include <mutex>
#include <limits>
#include <algorithm>
#include <unordered_map>
#include "fft_bluestein.hh"
#include "fft_any.hh"
#include "factor.hh"

struct BluesteinConfiguration
{
    cvector w, omega;
    std::size_t nb; // convolution size
};

static std::unordered_map<std::size_t, BluesteinConfiguration> data;
static std::mutex lock;
static const BluesteinConfiguration& configure_bluestein(std::size_t N, unsigned char maxfac, std::size_t m)
{
    constexpr unsigned B = std::numeric_limits<unsigned char>::digits, W = sizeof(std::size_t)*B;
    std::size_t key = N | (std::size_t(m?0:maxfac) << (W-B));
    std::unique_lock<std::mutex> lk(lock);
    if(auto i = data.find(key); i != data.end()) return i->second;

    std::size_t nb = 2*N - 1; // Minimal nb
    if(m)
        nb = m;
    else if(maxfac == 2)
        for(unsigned a=1; a<W; a<<=1) { nb = 1 + ((nb-1) | ((nb-1) >> a)); }
    else
        while(!small_factors_only(nb, maxfac)) { ++nb; } // Choose a hopefully-fast n'

    cvector w(N);
    exp_vector exps = ExpCircle(N*2);
    for(std::size_t k=0; k<N; ++k) { w[k] = exps(k*k, N*2); }

    cvector omega = w;
    for(auto& f: omega) f /= float(nb);
    omega.resize(nb);
    for(std::size_t k=1; k<N; ++k) omega[nb-k] = omega[k];

    lk.unlock();
    omega = FFT_any_fast(omega);
    lk.lock();

    BluesteinConfiguration cfg;
    cfg.nb    = nb;
    cfg.w     = std::move(w);
    cfg.omega = std::move(omega);

    return data.emplace(key, std::move(cfg)).first->second;
}

void FFT_bluestein(const complex* input, complex* output, std::size_t N,
                   std::size_t instride, std::size_t outstride,
                   unsigned char maxfac, std::size_t m)
{
    const auto& cfg = configure_bluestein(N, maxfac, m);
    const std::size_t nb = cfg.nb;
    cvector buf(nb*2);
    complex *buf1 = &buf[0], *buf2 = &buf[nb];

    // Multiply input by Bluestein sequence
    if(instride == 1)
        for(std::size_t k=0; k<N; ++k)
            buf1[k] = input[k] * cfg.w[k];
    else
        for(std::size_t k=0; k<N; ++k)
            buf1[k] = input[k*instride] * cfg.w[k];

    // Convolution: FFT
    exp_vector exps = ExpCircle(nb);
    FFT_any_fast(buf1, buf2, exps, nb, 1,1);

    // Convolution: pointwise multiplication
    for(std::size_t k=0; k<nb; ++k)
        buf2[k] = std::conj(buf2[k]) * cfg.omega[k];

    // Convolution: IFFT (which is calculated by regular FFT)
    FFT_any_fast(buf2, buf1, exps, nb, 1,1);

    // Multiply output by Bluestein sequence
    if(outstride == 1)
        for(std::size_t k=0; k<N; ++k)
            output[k] = std::conj(buf1[k]) * cfg.w[k];
    else
        for(std::size_t k=0; k<N; ++k)
            output[k*outstride] = std::conj(buf1[k]) * cfg.w[k];
}

cvector FFT_bluestein(const cvector& input)
{
    const auto N = input.size();
    cvector result(N);
    FFT_bluestein(&input[0], &result[0], N, 1, 1,   7);
    return result;
}
cvector FFT_bluestein_pow2(const cvector& input)
{
    const auto N = input.size();
    cvector result(N);
    FFT_bluestein(&input[0], &result[0], N, 1, 1,   2);
    return result;
}
cvector FFT_bluestein_fac2_3(const cvector& input)
{
    const auto N = input.size();
    cvector result(N);
    FFT_bluestein(&input[0], &result[0], N, 1, 1,   3);
    return result;
}
