#include "dft.hh"

cvector DFT_simple(const cvector& input)
{
    const auto N = input.size();
    auto exps = ExpCircle(N);
    cvector result(N);
    //#pragma omp parallel for
    for(std::size_t a = 0; a < N; ++a)
    {
        complex value {};
        for(std::size_t b = 0; b < N; ++b)
            value += input[b] * exps(a*b, N);
        result[a] = value;
    }
    return result;
}

#include "dft.hh"
#include <numbers>
#include <tuple>
double DFT_counter = 0;

void DFT(const complex* in, complex* out, exp_vector exps,
         std::size_t N, std::size_t is, std::size_t os)
{
    complex a,b,c,d, p13{-.5f,std::numbers::sqrt3_v<float>/-2}, p23=std::conj(p13), i{0,1};
    switch(N)
    {
        case 0: return;
        case 1: *out = *in; break;
        case 2: std::tie(a,b) = std::tuple(in[0], in[is]);
                std::tie(out[0], out[os]) = std::tuple{a+b, a-b};
                break;
        case 3: std::tie(a,b,c) = std::tuple(in[0], in[is], in[is*2]);
                std::tie(out[0], out[os], out[os*2]) = std::tuple(a+b+c, a+b*p13+c*p23, a+b*p23+c*p13);
                break;
        case 4: std::tie(a,b,c,d) = std::tuple(in[0], in[is], in[is*2], in[is*3]);
                std::tie(out[0], out[os], out[os*2], out[os*3])
                    = std::tuple(a+b+c+d, a-b*i-c+d*i, a-b+c-d, a+b*i-c-d*i);
                break;
        default:
            /*if(N >= 16) -- Cheat to get results fast from Tukey general
            {
                DFT_counter += std::pow(N, 2.037517);
                return;
            }*/
            //#pragma omp parallel for
            if(is == 1)
                if(os == 1)
                    for(std::size_t a = 0; a < N; ++a)
                    {
                        complex value {};
                        for(std::size_t z = 0, b = 0; b < N; ++b, z+=a)
                            value += in[b] * exps(z, N);
                        out[a] = value;
                    }
                else
                    for(std::size_t a = 0; a < N; ++a)
                    {
                        complex value {};
                        for(std::size_t z = 0, b = 0; b < N; ++b, z+=a)
                            value += in[b] * exps(z, N);
                        out[a*os] = value;
                    }
            else if(os == 1)
                for(std::size_t a = 0; a < N; ++a)
                {
                    complex value {};
                    for(std::size_t z = 0, b = 0; b < N; ++b, z+=a)
                        value += in[b*is] * exps(z, N);
                    out[a] = value;
                }
            else
                for(std::size_t a = 0; a < N; ++a)
                {
                    complex value {};
                    for(std::size_t z = 0, b = 0; b < N; ++b, z+=a)
                        value += in[b*is] * exps(z, N);
                    out[a*os] = value;
                }
    }
}

cvector DFT(const cvector& input)
{
    const auto N = input.size();
    if(N <= 1) return input;
    auto exps = ExpCircle(N);
    cvector result(N);
    DFT(&input[0], &result[0], exps, N, 1, 1);
    return result;
}
