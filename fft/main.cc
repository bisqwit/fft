#include "dft.hh"
#include "fft_fftw.hh"
#include "fft_radix2.hh"
#include "fft_tukey.hh"
#include "fft_rader.hh"
#include "fft_bluestein.hh"
#include "fft_any.hh"
#include "factor.hh"
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <type_traits>
#include <mutex>

#define VERIFY 1

static double avg_speed = 0; unsigned avg_mul = 0;

#include <fftw3.h>
static std::mutex fftw_lock, io_lock;
extern double DFT_counter;
void RunTest(cvector (*func)(const cvector&),
             const cvector& input,
             const cvector& expected,
             const char* what)
{
    bool is_fftw = (func == (cvector(*)(const cvector&))FFT_fftw);
    cvector output(input.size());

    fftwf_plan plan;

    if(is_fftw)
    {
      std::lock_guard<std::mutex> lk(fftw_lock);
      FFT_fftw(input, output, plan, 1);
    }

    ExpCircle(input.size());
    DFT_counter = 0;
    auto begin = std::chrono::steady_clock::now();

    if(is_fftw)
    {
        FFT_fftw(input, output, plan, 2);
    }
    else
    {
        output = func(input);
    }
    auto end   = std::chrono::steady_clock::now();
    double duration = std::chrono::duration<double>(end-begin).count();
    duration += 3.5289973e-9 * DFT_counter; DFT_counter = 0;

    const double threshold = 0.1;
    if(duration < threshold)// || true)
    {
        std::size_t reps = std::size_t(std::max(1., threshold / duration));
        // Start clock anew just in case the testcase has expensive first init
        begin = std::chrono::steady_clock::now();
        for(std::size_t a=0; a<reps; ++a)
        {
            if(is_fftw)
            {
                FFT_fftw(input, output, plan, 2);
            }
            else
            {
                output = func(input);
            }
        }
        end   = std::chrono::steady_clock::now();
        duration = std::chrono::duration<double>(end-begin).count();
        duration += 3.5289973e-9 * DFT_counter; DFT_counter = 0;
        duration /= reps;
    }

    if(is_fftw)
    {
      std::lock_guard<std::mutex> lk(fftw_lock);
      FFT_fftw(input, output, plan, 3);
    }

    //cvector (*temp)(const cvector&) = DFT;
    //if(func == temp) duration /= 48.0; // compensation for threads
    double diff = 0;
    const char* warning = "";
    #if VERIFY
    /* Create a hash of the generated output for validation purposes */
    for(std::size_t a=0; a<output.size(); ++a)
    {
        //double absdiff = std::abs(std::complex<double>(output[a])) - std::abs(std::complex<double>(expected[a]));
        //double argdiff = std::arg(std::complex<double>(output[a])) - std::arg(std::complex<double>(expected[a]));
        //diff += absdiff*absdiff + argdiff*argdiff;
        auto d = std::abs(std::complex<double>(output[a]) - std::complex<double>(expected[a]));
        diff += d*d;
    }
    diff /= output.size();
    if(diff > 1e-3f) warning = "\33[31m - WARNING\33[m";
    #endif

    std::lock_guard<std::mutex> lk(io_lock);
    std::printf("  Measuring N = %zu...\n", input.size());
    std::printf("    %16.8f ms (%s) - Mean squared difference: %.3e - TS %.6f ns%s\n", duration*1e3, what, diff, duration*1e9/input.size(), warning);
    avg_speed += duration*1e9/input.size(); ++avg_mul;
    //std::printf("%g ", duration*1e9/input.size());
    std::fflush(stdout);
}

#define tukey(how) [[maybe_unused]] static cvector FFT_tukey_##how(const cvector&input) { return FFT_tukey_fast(input, Tukey_Split_##how); }
tukey(r1)
tukey(half)
tukey(r2)
tukey(r1_rec_always)
tukey(half_rec_always)
tukey(r2_rec_always)
tukey(r1_rec_small_r1)
tukey(half_rec_small_r1)
tukey(r2_rec_small_r1)
tukey(r1_rec_small_r2)
tukey(half_rec_small_r2)
tukey(r2_rec_small_r2)
tukey(r1_rec_smallm1_r1)
tukey(half_rec_smallm1_r1)
tukey(r2_rec_smallm1_r1)
tukey(r1_rec_smallm1_r2)
tukey(half_rec_smallm1_r2)
tukey(r2_rec_smallm1_r2)

static unsigned sel = 0;

void Measure(std::size_t N)
{
    cvector test(N);

    for(std::size_t a = 0; a < N; ++a)
        test[a] = std::polar(1.f, 20.f * std::rand());

    //std::printf("  Measuring N = %zu...\n", N);
    cvector expected;
#if VERIFY
    { std::lock_guard<std::mutex> lk(fftw_lock);
      expected = FFT_fftw(test);
    }
#endif

    RunTest(FFT_fftw,   test, expected, "FFTW (Reference FFT from FFTW3)");
    if(true || N < 100000 || N%7==0)
    {
        //RunTest(DFT,        test, expected, "DFT (singlethread)");
    }
#if 1
    if(N > 2 && (N & (N-1)) == 0) // If N is a power of two
    {
        RunTest(FFT_radix2, test, expected, "FFT (Cooley-Tukey Radix2 version)");
        RunTest(FFT_radix2_fast, test, expected, "FFT (Cooley-Tukey Radix2 version fast)");
    }
#endif
    if(smallest_factor(N) != N) //N < 200000 || small_factors_only(N))
    {
        RunTest(FFT_tukey, test, expected, "FFT (Cooley-Tukey generic version)");

        //RunTest(FFT_tukey_r1,   test, expected, "FFT (Cooley-Tukey R1)");
        //RunTest(FFT_tukey_half, test, expected, "FFT (Cooley-Tukey Half)");
        if(sel <= 3) RunTest(FFT_tukey_r2,   test, expected, "FFT (Cooley-Tukey R2)");
        if(sel <= 3) RunTest(FFT_tukey_r1_rec_always,   test, expected, "FFT (Cooley-Tukey R1 RecY)");
        /*RunTest(FFT_tukey_half_rec_always, test, expected, "FFT (Cooley-Tukey Half RecY)");
        RunTest(FFT_tukey_r2_rec_always,   test, expected, "FFT (Cooley-Tukey R2 RecY)");
        RunTest(FFT_tukey_r1_rec_small_r1,   test, expected, "FFT (Cooley-Tukey R1 RecR1)");
        RunTest(FFT_tukey_half_rec_small_r1, test, expected, "FFT (Cooley-Tukey Half RecR1)");
        RunTest(FFT_tukey_r2_rec_small_r1,   test, expected, "FFT (Cooley-Tukey R2 RecR1)");
        RunTest(FFT_tukey_r1_rec_small_r2,   test, expected, "FFT (Cooley-Tukey R1 RecR2)");
        RunTest(FFT_tukey_half_rec_small_r2, test, expected, "FFT (Cooley-Tukey Half RecR2)");
        RunTest(FFT_tukey_r2_rec_small_r2,   test, expected, "FFT (Cooley-Tukey R2 RecR2)");
        RunTest(FFT_tukey_r1_rec_smallm1_r1,   test, expected, "FFT (Cooley-Tukey R1 RecM1R1)");
        RunTest(FFT_tukey_half_rec_smallm1_r1, test, expected, "FFT (Cooley-Tukey Half RecM1R1)");
        RunTest(FFT_tukey_r2_rec_smallm1_r1,   test, expected, "FFT (Cooley-Tukey R2 RecM1R1)");
        RunTest(FFT_tukey_r1_rec_smallm1_r2,   test, expected, "FFT (Cooley-Tukey R1 RecM1R2)");
        RunTest(FFT_tukey_half_rec_smallm1_r2, test, expected, "FFT (Cooley-Tukey Half RecMR2)");
        RunTest(FFT_tukey_r2_rec_smallm1_r2,   test, expected, "FFT (Cooley-Tukey R2 RecM1R2)");*/
    }
#if 1
    if(N > 1 && smallest_factor(N) == N)
    {
        RunTest(FFT_rader, test, expected, "FFT (Rader)");
        /*
        RunTest(FFT_bluestein,      test, expected, "FFT (Bluestein - prime)");
        RunTest(FFT_bluestein_pow2, test, expected, "FFT (Bluestein_pow2 - prime)");
        RunTest(FFT_bluestein_fac2_3, test, expected, "FFT (Bluestein_fac2_3 - prime)");
        */
    }
    if(N > 0)
    {
        RunTest(FFT_bluestein,      test, expected, "FFT (Bluestein - kaikki)");
        RunTest(FFT_bluestein_pow2, test, expected, "FFT (Bluestein_pow2 - kaikki)");
        RunTest(FFT_bluestein_fac2_3, test, expected, "FFT (Bluestein_fac2_3 - kaikki)");
    }
#endif
    if(sel <= 3) RunTest(FFT_any_fast, test, expected, "FFT (any-method new fast)");
    //fftwf_export_wisdom_to_filename("fftw-wisdom.dat");
}

template<typename... I>
void RunTests(const char* what, I... n) requires (std::is_integral_v<I> && ...)
{
    std::printf("Testing with %s\n", what);
    (Measure(n), ...);
}
#include "exp.hh"
#include <set>
#include <random>
int main(int argc, char **argv)
{
    /*{ int r = fftwf_import_wisdom_from_filename("fftw-wisdom.dat");
    std::printf("FFTW wisdom loaded, status=%d\n", r); }
    fftwf_set_timelimit(1.0);*/

    #pragma omp parallel for num_threads(32)
    for(unsigned n=1; n<=(1<<16); ++n) ExpCircle(n);
    #pragma omp parallel for num_threads(32)
    for(unsigned n=0; n<=25; ++n) ExpCircle(1<<n);
    //std::printf("ready\n"); std::fflush(stdout);

#if 0
    for(;;)
      #pragma omp parallel for num_threads(4)
      for(unsigned m=0; m<65536; ++m)
        for(unsigned n=1; n<=4; ++n)
            Measure(n);
#endif

#if 1
  #if 0 // Testing just Radix2
    //std::printf("Speeds: "); std::fflush(stdout);
    //#pragma omp parallel for num_threads(44) schedule(dynamic,1)
    for(unsigned m=0; m<1; ++m)
        for(unsigned power=0; power<=27; ++power)
        {
          unsigned n = 1 << power;
          //#pragma omp parallel for if(power>3) num_threads(28-power)
          for(unsigned q=0; q<3; ++q)
            Measure(n);
        }
    //std::printf(" - avg: %.6g ns/sample (%u tests)\n", avg_speed/avg_mul, avg_mul); std::fflush(stdout);
    return 0;
  #endif
  #if 1 // Testing just FFTW
    {std::vector<unsigned> q;
    sel = std::atoi(argv[2]);
    for(unsigned m=std::atoi(argv[1]), n=m, max = std::min(1048576*200u, n+800000); n<=max; ++n)
    {
        /*//std::printf("test %u\n", n);
        unsigned k = n;
        while(k%2==0) k/=2;
        while(k%3==0) k/=3;
        while(k%5==0) k/=5;
        while(k%7==0) k/=7;
        if(k == 1)*/
        if(small_factors_only(n))
            q.push_back(n);
    }
    for(unsigned m=std::atoi(argv[1]), n=m, max = std::min(1048576*200u, n+800000); n<=max; )
        // Tukey R1,R2.. multiset done up to 230
    {
        //std::size_t fac[16], nfac = get_factors(n, fac, 16);
        // up to 4000, go by +1
        // up to 32000, go by 
        q.push_back(n);
        //q.push_back(n);
        /*
        //#pragma omp parallel for num_threads(4)
        //for(unsigned q=0; q<4; ++q)
        //    Measure(n);
        //if(n >= 4096) ++n;
        //if(n >= 8000) ++n;
        if(n >= 5000) n+=2;
        if(n >= 7500) n+=2;
        if(n >= 10000) n+=2;
        if(n >= 15000) n+=2;
        if(n >= 20000) n+=2;
        if(n >= 30000) n+=2;
        if(n >= 40000) n+=2;
        if(n >= 60000) n+=2;
        if(n >= 80000) n+=2;
        if(n >= 120000) n+=2;
        if(n >= 160000) n+=2;
        if(n >= 240000) n+=2;
        if(n >= 320000) n+=2;
        if(n >= 480000) n+=2;
        if(n >= 640000) n+=2;
        if(n >= 960000) n+=2;
        if(n >= 1280000) n+=2;
        */
        //n += 1 + std::rand() % std::max(1, int(std::round(std::pow(n, 1.1677) / 7658.)));
        n += 1 + std::rand() % std::max(1, int(std::round(std::pow(n, 0.8451) / 24.504)));
        /*
        solve( [3000^x/y = 1.5, 1048576^x/y = 1400]);  --- saadaan x=1.6777, y=7658
        solve( [100^x/y = 2, 1048576^x/y = 5000]);  --- saadaan x=0.8451, y=24.504
            B^x/y = K
            B^x   = y*K
            log(B^x) = log(y*K)
            x*log(B) = log(y)+log(K)
            x*log(B)-log(y) = log(K)
            x-log(y)/log(B) = log(K-B)
            x = log(K-B) + log(y)/log(B)
            y = B^x/K
            Sijoitetaan y:
            x1 = log(K1-B1) + log(B2^x2/K2)/log(B1)
               = log(K1-B1) + (x2*log(B2)-log(K2))/log(B1)
               = log(K1-B1) + x2*log(B2)/log(B1) - log(K2)/log(B1)
        */
    }
    std::sort(q.begin(), q.end()); q.erase(std::unique(q.begin(), q.end()), q.end());
    std::printf("%zu tests:", q.size());
    for(auto v: q) std::printf(" %zu", v);
    std::printf("\n");    
    std::fflush(stdout);
    #pragma omp parallel for num_threads(35) schedule(dynamic,1)
    for(unsigned a=0; a<q.size(); ++a)
        Measure(q[a]);
    return 0;}
  #endif
  #if 0 // Testing just small factors
    {std::vector<unsigned> q;
    for(unsigned n=1; n<=1048576*2; ++n)
        //if(small_factors_only(n))
    {
        q.push_back(n);
        q.push_back(n);
        q.push_back(n);
        q.push_back(n);
    }
    std::printf("%zu tests\n", q.size());
    //#pragma omp parallel for num_threads(4) schedule(dynamic,2)
    for(unsigned a=0; a<q.size(); ++a)
        Measure(q[a]);
    return 0;}
  #endif
  #if 0 // Testing just primes
    {std::vector<unsigned> q;
    for(unsigned m=std::atoi(argv[1]), n=m, max = std::min(1048576*8u, n+30000); n<max; )
    {
        while(smallest_factor(n)!=n) ++n;
        q.push_back(n);
        //q.push_back(n);
        n += 1;
        n += std::rand() % std::max(1, int(std::round(n / (6000. / std::log10(n)))));
    }
    std::printf("%zu tests\n", q.size());
    //#pragma omp parallel for num_threads(6) schedule(dynamic,2)
    for(unsigned a=0; a<q.size(); ++a)
        Measure(q[a]);
    return 0;}
  #endif
  #if 0 // Testing just DFT
    {std::vector<unsigned> t;
    //for(unsigned n=0; n<80; ++n) t.push_back(0);
    /*for(unsigned n=0; n<=26; ++n)
    for(unsigned m=0; m<400; ++m)
        t.push_back(1<<n);*/
    if(false)for(unsigned n=1; n<=1048576*64; ++n)
        if(small_factors_only(n))
            t.push_back(n);
    for(unsigned m=0; m<10; ++m)
    for(unsigned n=1; n<=1048576*64; )
    {
        t.push_back(n);
        if(n < 13000) ++n; else n += std::max(1u, std::rand() % unsigned(1+10*std::log(n)));
    }
    for(unsigned a=0; a<8; ++a)
    for(unsigned m=0; m<=26; ++m)
    {
    unsigned l = 1<<((26-m)/2);//6300 / (m+1);// 
    //unsigned l = t.size();
    /*
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(t.begin(), t.end(), g);
    */
    l = std::max(8u, std::min(32u, l));
    //#pragma omp parallel for num_threads(40) schedule(dynamic,32)
    for(unsigned n=0; n<l; ++n)
    {
        //Measure(t[n]);
        Measure(1<<m);
    }
    }
    return 0;}
#endif

    std::set<unsigned> t;
    std::size_t temp[16];
    for(unsigned n=1; n<=1048576; ++n)
        if(false
        || smallest_factor(n)==n
        || small_factors_only(n)
        || get_factors(n, temp, 8) <= 3
          )
            t.insert(n);
    std::vector<unsigned> tests(t.begin(), t.end());
    //std::printf("Making %zu exp circles\n", tests.size()); std::fflush(stdout);
    //#pragma omp parallel for num_threads(16)
    //for(unsigned n=0; n<tests.size(); ++n) ExpCircle(tests[n]);

    std::random_device rd;
    std::mt19937 g(rd());
    //std::shuffle(tests.begin(), tests.end(), g);

    //for(unsigned n=1; n<2048; ++n) Measure(n);
    std::printf("Measuring\n"); std::fflush(stdout);
    for(unsigned n=0; n<tests.size(); ++n) Measure(tests[n]);
#endif
    //Measure(41);return 0;
    /*
    RunTests("powers of 2", 16, 256, 4096, 16384, 65536, 262144);
    RunTests("powers of 3", 9, 81, 729, 6561, 59049, 177147);
    RunTests("powers of 5", 25, 625, 15625, 78125);
    RunTests("mixed factors", 30, 900, 18900, 147000);
    RunTests("primes", 3, 7, 17, 173, 971, 2113, 5393, 37813, 59359, 139901, 200183, 401987);
    RunTests("primes", 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271);
    RunTests("primes", 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601);
    RunTests("primes", 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953);
    RunTests("primes", 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279);
    RunTests("primes", 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597);
    RunTests("primes", 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, 1901, 1907, 1913, 1931, 1933);
    RunTests("primes", 1949, 1951, 1973, 1979, 1987, 1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, 2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129, 2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, 2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281);
    RunTests("primes", 2287, 2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, 2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, 2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, 2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, 2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, 2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741);
    */
}
