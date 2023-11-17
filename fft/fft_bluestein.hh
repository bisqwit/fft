#include "defs.hh"

// 2, 3, 5, 7,   11, 13, 17, 19
void FFT_bluestein(const complex* input, complex* output, std::size_t N,
                   std::size_t instride, std::size_t outstride,
                   unsigned char maxfac, std::size_t m = 0);

cvector FFT_bluestein(const cvector& input);
cvector FFT_bluestein_pow2(const cvector& input);
cvector FFT_bluestein_fac2_3(const cvector& input);
