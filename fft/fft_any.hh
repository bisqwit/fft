#include "defs.hh"
#include "exp.hh"

cvector FFT_any_fast(const cvector& input);
void FFT_any_fast(
    const complex*    input,
    complex*          output,
    exp_vector        exps,
    std::size_t N, std::size_t instride, std::size_t outstride);
