#include "dft.hh"
#include "fft_radix2.hh"

cvector FFT_radix2(const cvector& input)
{
    const auto N = input.size();
    if(N <= 1) return input;

    std::size_t half = N/2;

    exp_vector exps = ExpCircle(N);
    cvector even(half), odd(half);
    for(std::size_t a = 0; a < half; ++a) even[a] = input[a*2];
    for(std::size_t a = 0; a < half; ++a) odd[a]  = input[a*2+1];
    even = FFT_radix2(even);
    odd  = FFT_radix2(odd);
    for(std::size_t a = 0; a < half; ++a) odd[a] *= exps(a,N);

    even.resize(N);
    for(std::size_t a = 0; a < half; ++a)
    {
        even[a+half] = even[a] - odd[a];
        even[a]      = even[a] + odd[a];
    }
    return even;
}

#include "fft_any.hh"

cvector FFT_radix2_fast(const cvector& input)
{
    const auto N = input.size();
    exp_vector exps = ExpCircle(N);
    cvector output(N);
    FFT_radix2_fast(&input[0], &output[0], exps, N, 1,1);
    return output;
}

void FFT_radix2_fast(
    const complex*    input,
    complex*          output,
    exp_vector        exps,
    std::size_t N, std::size_t instride, std::size_t outstride)
{
    if(N <= 4) { DFT(input,output,exps,N,instride,outstride); return; }

    // Because of cache efficiency reasons, use a smaller exps table in deeper recursions.
    //if(N >= 256 && exps.data->size() > N) exps = ExpCircle(N);

    if(N%2 != 0) [[unlikely]]
    {
        FFT_any_fast(input,output,exps, N, instride,outstride);
        return;
    }

    std::size_t half = N/2;
    complex* even = output, *odd = output + half*outstride;
    FFT_radix2_fast(input,          even, exps, half, instride*2, outstride);
    FFT_radix2_fast(input+instride, odd,  exps, half, instride*2, outstride);
    for(std::size_t a = 0; a < half; ++a) odd[a*outstride] *= exps(a, N);

    for(std::size_t a = 0; a < half; ++a)
    {
        complex part1 = even[a*outstride] + odd[a*outstride];
        complex part2 = even[a*outstride] - odd[a*outstride];
        even[a*outstride] = part1;
        odd[a*outstride]  = part2;
    }
}
