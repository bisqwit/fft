#include "defs.hh"

cvector FFT_fftw(const cvector& input);
extern "C" {
    struct fftwf_plan_s;
}
void FFT_fftw(const cvector& input, cvector& output, struct fftwf_plan_s*& plan, unsigned phase);
