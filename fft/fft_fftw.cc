#include "fft_fftw.hh"
#include <fftw3.h>
#include <mutex>

static std::mutex mutex;
cvector FFT_fftw(const cvector& input)
{
    const std::size_t N = input.size();
    cvector output(N);
    fftwf_plan plan;

    {std::lock_guard<std::mutex> lk(mutex);
     plan = fftwf_plan_dft_1d(N,
        (fftwf_complex*)&input[0],
        (fftwf_complex*)&output[0],
        FFTW_FORWARD, FFTW_ESTIMATE);
    }
    fftwf_execute(plan);
    {std::lock_guard<std::mutex> lk(mutex);
     fftwf_destroy_plan(plan);
    }
    return output;
}

void FFT_fftw(const cvector& input, cvector& output, fftwf_plan& plan, unsigned phase)
{
    const std::size_t N = input.size();
    switch(phase)
    {
        case 1:
            plan = fftwf_plan_dft_1d(N,
                (fftwf_complex*)&input[0],
                (fftwf_complex*)&output[0],
                FFTW_FORWARD, FFTW_ESTIMATE);
            break;
        case 2:
            fftwf_execute(plan);
            break;
        case 3:
            fftwf_destroy_plan(plan);
    }
}
