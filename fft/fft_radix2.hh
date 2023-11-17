#include "defs.hh"

cvector FFT_radix2(const cvector& input);
cvector FFT_radix2_fast(const cvector& input);

void FFT_radix2_fast(
    const complex*    input,
    complex*          output,
    exp_vector        exps,
    std::size_t N, std::size_t instride, std::size_t outstride);
