#include "dft.hh"
#include "fft_tukey.hh"
#include "fft_any.hh"
#include "factor.hh"

cvector FFT_tukey(const cvector& input, exp_vector exps)
{
    const auto N = input.size(), r2 = smallest_factor(N), r1 = N/r2;
    if(r1 <= 1) { return FFT_any_fast(input); }// Factorization failed, revert to DFT.

    cvector result(N), temp(N);
    //#pragma omp parallel for
    for(std::size_t k0 = 0; k0 < r2; ++k0)
    {
        cvector slice(r1);
        for(std::size_t k1 = 0; k1 < r1; ++k1) slice[k1] = input[k1*r2 + k0];
        auto subresult = FFT_tukey(slice);
        for(std::size_t k1 = 0; k1 < r1; ++k1) temp[k0*r1 + k1] = subresult[k1];
    }

    //#pragma omp parallel for collapse(2)
    for(std::size_t n1 = 0; n1 < r2; ++n1)
    for(std::size_t n0 = 0; n0 < r1; ++n0)
    {
        complex value = 0;
        for(std::size_t k0 = 0; k0 < r2; ++k0)
            value += temp[k0*r1+n0] * exps((n1*r1+n0)*k0, N);
        result[n1*r1 + n0] = value;
    }
    return result;
}

cvector FFT_tukey(const cvector& input)
{
    return FFT_tukey(input, ExpCircle(input.size()));
}

template<int splitdir, int use_rec_method>
void FFT_tukey_fast(
    const complex* input,
    complex* output,
    exp_vector exps,
    std::size_t N,
    std::size_t instride, std::size_t outstride)
{
    //if(N >= 256 && exps.data->size() > N) { exps = ExpCircle(N); }

    std::size_t r1,r2;
    switch(splitdir)
    {
        default:
        case 0: r1 = half_factors(N); r2 = N/r1; break;
        case -1: r1 = smallest_factor(N); r2 = N/r1; break;
        case 1: r2 = smallest_factor(N); r1 = N/r2; break;
    }
    if(r1 == 1 || r2 == 1)
    {
        // Factorization failed
        FFT_any_fast(input, output, exps, N, instride, outstride);
        return;
    }

    bool use_rec = false;
    switch(use_rec_method)
    {
        case 0: break;
        case 9: use_rec = true; break;
        case -1: use_rec = small_factors_only(r1); break;
        case  1: use_rec = small_factors_only(r2); break;
        case -2: use_rec = smallest_factor(r1)==r1 && small_factors_only(r1-1); break;
        case  2: use_rec = smallest_factor(r2)==r2 && small_factors_only(r2-1); break;
    }

    if(use_rec)
    {
        cvector temp(N);
        //#pragma omp parallel for
        for(std::size_t k0 = 0; k0 < r2; ++k0)
            FFT_tukey_fast<splitdir,use_rec_method>(input + k0*instride, &temp[k0], exps, r1, r2*instride, r2);
        for(std::size_t n0=0; n0<r1; ++n0)
            for(std::size_t k0 = 0; k0 < r2; ++k0)
                temp[k0 + n0*r2] *= exps(k0*n0, N);

        //#pragma omp parallel for
        for(std::size_t n0 = 0; n0 < r1; ++n0)
            FFT_tukey_fast<splitdir,use_rec_method>(&temp[n0*r2], &output[n0*outstride], exps, r2, 1,r1*outstride);
    }
    else
    {
        cvector temp(N);
        //#pragma omp parallel for
        for(std::size_t k0 = 0; k0 < r2; ++k0)
            FFT_tukey_fast<splitdir,use_rec_method>(input + k0*instride, &temp[k0], exps, r1, r2*instride, r2);

        //#pragma omp parallel for
        for(std::size_t n1 = 0; n1 < r2; ++n1)
            for(std::size_t n0 = 0; n0 < r1; ++n0)
            {
                complex value = 0;
                for(std::size_t k0 = 0; k0 < r2; ++k0)
                    value += temp[k0 + n0*r2] * exps((n1*r1+n0)*k0, N);
                output[(n1*r1 + n0)*outstride] = value;
            }
    }
}

cvector FFT_tukey_fast(const cvector& input, TukeySplit how)
{
    const auto N = input.size();
    exp_vector exps = ExpCircle(N);
    cvector result(N);
    switch(how)
    {
    case Tukey_Split_r1:                  FFT_tukey_fast<-1,0>(&input[0], &result[0], exps, N, 1, 1); break;
    case Tukey_Split_half:                FFT_tukey_fast< 0,0>(&input[0], &result[0], exps, N, 1, 1); break;
    case Tukey_Split_r2:                  FFT_tukey_fast< 1,0>(&input[0], &result[0], exps, N, 1, 1); break;
    case Tukey_Split_r1_rec_always:       FFT_tukey_fast<-1,9>(&input[0], &result[0], exps, N, 1, 1); break;
    case Tukey_Split_half_rec_always:     FFT_tukey_fast< 0,9>(&input[0], &result[0], exps, N, 1, 1); break;
    case Tukey_Split_r2_rec_always:       FFT_tukey_fast< 1,9>(&input[0], &result[0], exps, N, 1, 1); break;
    case Tukey_Split_r1_rec_small_r1:     FFT_tukey_fast<-1,-1>(&input[0], &result[0], exps, N, 1, 1); break;
    case Tukey_Split_half_rec_small_r1:   FFT_tukey_fast< 0,-1>(&input[0], &result[0], exps, N, 1, 1); break;
    case Tukey_Split_r2_rec_small_r1:     FFT_tukey_fast< 1,-1>(&input[0], &result[0], exps, N, 1, 1); break;
    case Tukey_Split_r1_rec_small_r2:     FFT_tukey_fast<-1,1>(&input[0], &result[0], exps, N, 1, 1); break;
    case Tukey_Split_half_rec_small_r2:   FFT_tukey_fast< 0,1>(&input[0], &result[0], exps, N, 1, 1); break;
    case Tukey_Split_r2_rec_small_r2:     FFT_tukey_fast< 1,1>(&input[0], &result[0], exps, N, 1, 1); break;
    case Tukey_Split_r1_rec_smallm1_r1:   FFT_tukey_fast<-1,-2>(&input[0], &result[0], exps, N, 1, 1); break;
    case Tukey_Split_half_rec_smallm1_r1: FFT_tukey_fast< 0,-2>(&input[0], &result[0], exps, N, 1, 1); break;
    case Tukey_Split_r2_rec_smallm1_r1:   FFT_tukey_fast< 1,-2>(&input[0], &result[0], exps, N, 1, 1); break;
    case Tukey_Split_r1_rec_smallm1_r2:   FFT_tukey_fast<-1,2>(&input[0], &result[0], exps, N, 1, 1); break;
    case Tukey_Split_half_rec_smallm1_r2: FFT_tukey_fast< 0,2>(&input[0], &result[0], exps, N, 1, 1); break;
    case Tukey_Split_r2_rec_smallm1_r2:   FFT_tukey_fast< 1,2>(&input[0], &result[0], exps, N, 1, 1); break;
    }
    return result;
}
