#include "defs.hh"

void FFT_rader(const complex* input, complex* output, std::size_t N,
               std::size_t instride, std::size_t outstride);

cvector FFT_rader(const cvector& input);
