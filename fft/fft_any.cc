#include "dft.hh"
#include "fft_tukey.hh"
#include "fft_rader.hh"
#include "fft_radix2.hh"
#include "fft_bluestein.hh"
#include "factor.hh"

#include <unordered_map>
#include <unordered_set>
#include <mutex>
#include <bit>

/* o(name, r1tab,r2val,  rec1val,rec1mul, rec2val,rec2mul, mul,add,temp,straightness, condition) */
#define LIST_METHODS(o) \
    /* DFT: All memory access is straight */ \
    o(DFT,        z,0,       0,0,   0,0,   N<5 ? 4*(N==3) : (N*N), N*(N-(N<5)), 0, 1.0, true) \
    /* Radix2: Two interleaved N/2 FFTs, three straight N/2 loops */ \
    o(Radix2,     z,0,       N/2,2, 0,0,   N/2, N, 0,     3/5.0,   smallest==2) \
    /* Rader: Two straight N-1 FFTs, two N-1 messed up loops, one N-1 straight loop, */ \
    o(Rader,      z,0,       N-1,2, 0,0,   N-1, 0, 1,     2/5.0,   smallest==N) \
    /* Bluestein: Two straight nb FFTs, two straight N loops, one straight nb loop */ \
    o(Bluestein,  M,0,r1,2, 0,0,    N*2+r1, 0, 3,  1.0,     r1 && N<original*3 && !resolving.contains(r1)) \
    /* Rec1: r2-interleaved r1-FFT times r2, one r2 straight loop times N */ \
    o(Tukey_Rec1, fac,N/r1, r1,r2, 0,0,   N*r2,N*r2,1,   (N*1.0/(N+r2)),       smallest!=N) \
    /* Rec2: r2-interleaved r1-FFT times r2, one r2 straight loop times r1, r1-interleaved r2-FFT times r1 */ \
    o(Tukey_Rec2, fac,N/r1, r1,r2, r2,r1, N,0,1,         (r1*1.0/(r2+r1+r1)),  smallest!=N)

#define o(name,T,R,c,d,E,F,m,a,t,s,condition) name,
enum class AnyMethod: unsigned { LIST_METHODS(o) };
#undef o

struct AnyConfiguration
{
    std::size_t r1, r2, mult,add,temp;
    float straight;
    AnyMethod method;
};
static std::unordered_map<std::size_t, AnyConfiguration> data;
static std::unordered_map<std::size_t, std::size_t> cache;
static std::unordered_set<std::size_t> resolving;
static std::size_t original;
static std::mutex lock;
static AnyConfiguration& configure_any(std::size_t N)
{
    std::unique_lock<std::mutex> lk(lock);
    if(auto i = data.find(N); i != data.end()) return i->second;

    if(N % 2 == 0) //(N & (N-1)) == 0)
    {
        // For powers of two and even numbers, just give Radix2
        AnyConfiguration cfg{};
        cfg.method = AnyMethod::Radix2;
        return data.emplace(N, std::move(cfg)).first->second;
    }
    if(std::size_t r = smallest_factor(N); r != N && r < 20)
    {
        AnyConfiguration cfg{};
        cfg.r1 = r;
        cfg.r2 = N/r;
        cfg.method = AnyMethod::Tukey_Rec2;
        return data.emplace(N, std::move(cfg)).first->second;
    }

    if(resolving.empty()) original = N;
    resolving.insert(N);

    std::size_t ftab[16], fsize = get_factors(N, ftab, 16);
    std::vector<std::size_t> M, fac, z(1);
#if 0
    for(std::size_t nb=2*N-1; ; ++nb)
    {
        std::size_t m = 0;
        for(std::size_t t = nb;;)
        {
            if(auto i = cache.find(t); i != cache.end()) { m |= i->second; break; }
            else if(t % 2 == 0) { m |= 1; t /= 2; }
            else if(t % 3 == 0) { m |= 2; t /= 3; }
            else if(t % 5 == 0) { m |= 4; t /= 5; }
            else if(t % 7 == 0) { m |= 8; t /= 7; }
            else { m = cache.emplace(nb, m+256*t).first->second; break; }
        }
        if(unsigned facmask = m&255, b = std::popcount(facmask); m < 2*256 && b >= 1 && b <= 2)
        {
            M.push_back(nb);
            if(facmask == 1) break; // Up to next power of two
        }
    }
#else
    std::size_t nb = 2*N-1; // Round up to a power of two
    nb--;
    nb |= (nb>>1);
    nb |= (nb>>2);
    nb |= (nb>>4);
    nb |= (nb>>8);
    nb |= (nb>>16);
    nb |= (nb>>32);
    nb++;
    M.push_back(nb);
#endif
    std::size_t smallest = smallest_factor(N), best = ~std::size_t{};
    for(std::size_t f=0; f< std::size_t(1<<fsize); ++f)
        if(std::size_t product = 1; true) //std::popcount(f) <= 8)
        {
            for(std::size_t b=0; b<fsize; ++b) { if(f & (1<<b)) { product *= ftab[b]; } }
            if(product != 1 && product != N) fac.push_back(product);
        }

    AnyConfiguration cfg;
    #define o(name, T,R,c,d,E,F,m,a,t,s, condition) for(std::size_t f=0; f<T.size(); ++f) \
    { \
        std::size_t mult=0, add=0, temp=0, r1 = T[f], r2 = R; \
        if(!(condition)) \
        { \
            /*if(!f) std::printf("For N=%zu, ignoring %s because !(%s) - smallest=%zu\n", N,#name, #condition, smallest);*/ \
            continue; \
        } \
        float straight = (s); \
        lk.unlock(); \
        if(c) { auto& S = configure_any(c); mult += (d)*S.mult; add += (d)*S.add; temp += (d)*S.temp; straight *= S.straight; } \
        if(E) { auto& S = configure_any(E); mult += (F)*S.mult; add += (F)*S.add; temp += (F)*S.temp; straight *= S.straight; } \
        lk.lock(); \
        if(auto i = data.find(N); i != data.end()) return i->second; \
        mult += (m); \
        add  += (a); \
        temp += (t); \
        /* Each complex multiplication consists of 4 multiplications and 2 additions */ \
        /* Each complex addition consists of two real additions */ \
        std::size_t real_mul = mult * 4, real_add = mult * 2 + add * 2; \
        float influence = 1/*0.75 + (1-straight)*0.25*/; \
        std::size_t slowness = (real_mul + real_add) * influence + temp*64; \
        /*std::printf("For N=%zu, %s (r1=%zu r2=%zu) gets %zu mult (%zu real, %zu add (%zu real), straight=%g, %zu temp\n", \
            N, #name, r1,r2, mult,real_mul, add, real_add, straight, temp);*/ \
        if(best == ~std::size_t{} || (slowness < best || (slowness == best && straight > cfg.straight))) \
        { \
            best = slowness; \
            cfg.method   = AnyMethod::name; \
            cfg.r1       = r1; \
            cfg.r2       = r2; \
            cfg.mult = mult; \
            cfg.add = add; \
            cfg.temp = temp; \
            cfg.straight = straight; \
        } \
    }
    LIST_METHODS(o)
    #undef o
    //cfg.buf.resize(cfg.r1 * cfg.r2);

    resolving.erase(N);
    /*std::size_t real_mul = chose_mult * 4, real_add = chose_mult * 3 + chose_add * 2;
    std::printf("N=%zu uses %s r1=%zu r2=%zu (%zu mult (%zu real), %zu add (%zu real), straight=%g, %zu temp)\n",
        N, chose, cfg.r1,cfg.r2,
        chose_mult,real_mul, chose_add,real_add, chose_straight, chose_temp);*/
    return data.emplace(N, std::move(cfg)).first->second;
}

void FFT_any_fast(
    const complex* input,
    complex* output,
    exp_vector exps,
    std::size_t N,
    std::size_t instride, std::size_t outstride)
{
    const auto& cfg = configure_any(N);
    std::size_t r1 = cfg.r1, r2 = cfg.r2;
    //auto& buf = cfg.buf;
    cvector buf(N);

    switch(cfg.method)
    {
        case AnyMethod::DFT:       DFT(input, output, exps, N, instride, outstride); break;
        case AnyMethod::Rader:     FFT_rader(input, output, N, instride, outstride); break;
        case AnyMethod::Bluestein: FFT_bluestein(input, output, N, instride, outstride, 0, r1); break;
        case AnyMethod::Radix2:
        {
            // Because of cache efficiency reasons, use a smaller exps table in deeper recursions.
            //if(N >= 256 && exps.data->size() > N) { exps = ExpCircle(N); }
            std::size_t half = N/2;
            complex* even = output, *odd = output + half*outstride;
            FFT_any_fast(input,          even, exps, half, instride*2, outstride);
            FFT_any_fast(input+instride, odd,  exps, half, instride*2, outstride);
            if(outstride == 1)
            {
                for(std::size_t a = 0; a < half; ++a) odd[a] *= exps(a, N);
                for(std::size_t a = 0; a < half; ++a)
                {
                    complex part1 = even[a] + odd[a];
                    complex part2 = even[a] - odd[a];
                    even[a] = part1;
                    odd[a]  = part2;
                }
            }
            else
            {
                for(std::size_t a = 0; a < half; ++a) odd[a*outstride] *= exps(a, N);
                for(std::size_t a = 0; a < half; ++a)
                {
                    complex part1 = even[a*outstride] + odd[a*outstride];
                    complex part2 = even[a*outstride] - odd[a*outstride];
                    even[a*outstride] = part1;
                    odd[a*outstride]  = part2;
                }
            }
            return;
        }
        case AnyMethod::Tukey_Rec1:
            // Because of cache efficiency reasons, use a smaller exps table in deeper recursions.
            //if(N >= 256 && exps.data->size() > N) { exps = ExpCircle(N); }

            //#pragma omp parallel for
            for(std::size_t k0 = 0; k0 < r2; ++k0)
                FFT_any_fast(input + k0*instride, &buf[k0], exps, r1, r2*instride, r2);

            //#pragma omp parallel for
            if(outstride == 1)
            {
                for(std::size_t p=0, n1 = 0; n1 < r2; ++n1)
                    for(std::size_t q = 0, n0 = 0; n0 < r1; ++n0, ++p)
                    {
                        complex value = 0;
                        for(std::size_t z = 0, k0 = 0; k0 < r2; ++k0, z += p, ++q)
                            value += buf[q] * exps(z, N);
                        output[p] = value;
                    }
            }
            else
            {
                for(std::size_t p=0, n1 = 0; n1 < r2; ++n1)
                    for(std::size_t q = 0, n0 = 0; n0 < r1; ++n0, ++p)
                    {
                        complex value = 0;
                        for(std::size_t z = 0, k0 = 0; k0 < r2; ++k0, z += p, ++q)
                            value += buf[q] * exps(z, N);
                        output[p*outstride] = value;
                    }
            }
            return;
        case AnyMethod::Tukey_Rec2:
            // Because of cache efficiency reasons, use a smaller exps table in deeper recursions.
            //if(N >= 256 && exps.data->size() > N) { exps = ExpCircle(N); }

            //#pragma omp parallel for
            for(std::size_t k0 = 0; k0 < r2; ++k0)
                FFT_any_fast(input + k0*instride, &buf[k0], exps, r1, r2*instride, r2);

            for(std::size_t q = 0, n0=0; n0<r1; ++n0)
                for(std::size_t z = 0, k0 = 0; k0 < r2; ++k0, z += n0, ++q)
                    buf[q] *= exps(z, N);

            //#pragma omp parallel for
            for(std::size_t q = 0, n0 = 0; n0 < r1; ++n0, q += r2)
                FFT_any_fast(&buf[q], &output[n0*outstride], exps, r2, 1,r1*outstride);
            return;
    }
}

cvector FFT_any_fast(const cvector& input)
{
    const auto N = input.size();
    exp_vector exps = ExpCircle(N);
    cvector result(N);
    FFT_any_fast(&input[0], &result[0], exps, N, 1, 1);
    return result;
}
