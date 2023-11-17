#include <algorithm>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdio>

#include "trend.hh"
#include <set>

#include "fft_any.hh"
#include "fft_fftw.hh"
#include <gd.h>

static double Burg(const float* data, std::size_t N, float* coeff, std::size_t M)
{
    for(std::size_t a=0; a<M; ++a) coeff[a] = 0;
    std::vector<float> b1(N), b2(N), aa(M);
    double p = 0;
    for(std::size_t a=0; a<N; ++a) p += data[a] * data[a];
    double xms = p / N;
    b1[0]   = data[0];
    b2[N-2] = data[N-1];
    for(std::size_t a=1; a<N-1; ++a) b1[a] = b2[a-1] = data[a];
    for(std::size_t a=0; a<M; ++a)
    {
        double num = 0, denum = 0;
        for(std::size_t b=0; b<N-a; ++b)
            num   += b1[b]*b2[b],
            denum += b1[b]*b1[b] + b2[b]*b2[b];
        coeff[a] = 2*num/denum;
        xms *= 1 - coeff[a]*coeff[a];
        for(std::size_t b=1; b<=a; ++b)
            coeff[b-1] = aa[b-1] - coeff[a] * aa[a-b];
        if(a < M-1)
        {
            for(std::size_t b=0; b<=a; ++b) aa[b] = coeff[b];
            for(std::size_t b=1; b<=N-(a+1)-1; ++b)
                b1[b-1] -= aa[a] * b2[b-1],
                b2[b-1] = b2[b-1+1] - aa[a]*b1[b-1+1];
        }
    }
    return xms;
}

static void Analyze(const float* data, std::size_t size, std::size_t rate)
{
    double length    = size / double(rate); /* Length in seconds */
    double increment = length/4000, window = 0.013;
    std::size_t maxfreq = 44100, interesting[2] = {0,5000}, x = 0; // Max 8000 Hz
    std::size_t columns = std::ceil(length/increment)+1;
    std::vector<double> result(maxfreq * columns);
    for(double start = 0; start < length; start += increment, ++x)
    {
        /* Generate complex sample */
        cvector sample(maxfreq);
        for(std::size_t p=start*rate, a=0; a<maxfreq; ++a,++p) { sample[a] = p < size ? data[p] : 0; }
        /* Apply Gaussian window */
        std::size_t wnd = window*rate;
        double wndsum = 0;
        for(std::size_t a=0; a<maxfreq; ++a)
        {
            double edge  = std::exp(-12.0), p = double(a) / wnd - 0.5;
            double value = (a<wnd ? (std::exp(-48.*p*p) - edge) / (1-edge) : 0);
            sample[a] *= value; wndsum += value;
        }
        // Normalize it
        for(std::size_t a=0; a<maxfreq; ++a) sample[a] /= wndsum;
        /* Convert */
        sample = FFT_fftw(sample);
        /* Convert complex FFT into real */
        std::vector<float> real(maxfreq), realpre(maxfreq);
        for(std::size_t a=0; a<maxfreq; ++a) { real[a] = realpre[a] = std::abs(sample[a]/* * (a?sample[maxfreq-a]:complex(1))*/); }
        /* Convert absolute signal to dB */
        for(std::size_t a=0; a<maxfreq/2; ++a) real[a] = 10 * std::log10((real[a] + 1e-30) / 4e-10) + 6*std::log2((a+!a)/1000.0);
        /*
        float min=4e40,max=-4e40, dbrange = 10;
        for(std::size_t a=1; a<maxfreq/2; ++a) min = std::min(min, real[a]);
        for(std::size_t a=1; a<maxfreq/2; ++a) max = std::max(max, real[a]);
        //min = std::max(min, max - dbrange);
        //for(std::size_t a=0; a<maxfreq; ++a) min = std::min(min, real[a]);
        printf("min=%g, max=%g\n", min,max);
        */
        float min = 55, max = 90;
        for(std::size_t a=0; a<maxfreq; ++a) real[a] = std::clamp((real[a] - min) / (max-min), 0.f, 1.f);
        //for(std::size_t a=0; a<maxfreq; ++a) real[a] = std::clamp(real[a], 0.f, 70.f) / 70.;
        /* Render into result */
        for(std::size_t a=0; a<maxfreq/2; ++a)
        {
            result[x*maxfreq + a]             = real[a];
            result[x*maxfreq + a + maxfreq/2] = realpre[a];
        }
    }
    /* Convert the result into a picture */
    gdImagePtr im = gdImageCreateTrueColor(columns, interesting[1]-interesting[0]);
    auto color = [](double f)
    {
        f = std::pow(1-f, 1.2f); int b = f*768, r = b-256, g = b-256;
        //int r = (1-f)*255, g=r,b=g;
        return (std::clamp(r,0,255)<<16)+(std::clamp(g,0,255)<<8)+(std::clamp(b,0,255)<<0);
    };
    for(std::size_t x=0; x<columns; ++x)
    {
        std::vector<double> freq, data;
        for(std::size_t a=interesting[0]; a<=interesting[1]; ++a)
        {
            freq.push_back( a );
            data.push_back( result[x*maxfreq + a] );
        }
        TrendPoly<12> est;
        est.estimate(freq, data);
        /* For zero points for the first derivative */
        auto diff = est.differentiate();
        auto diff2 = diff.differentiate();
        std::set<int> cand;
        printf("%zu: f = %s, f' = %s, f'' = %s\n",
            freq.size(),
            est.translate().c_str(), diff.translate().c_str(), diff2.translate().c_str());
        for(auto f: freq)
        {
            double f2 = diff.solve_newton(f, 1e-12, 200);
            //printf("\tguess for %g=%g\n", f, f2);
            if(f2 >= 50)
            {
                if(diff2.predict(f2) < 0) { cand.insert(std::round(f2)); }
            }
        }
        printf("guess @x=%zu:", x); for(auto c: cand) printf(" %d", c); printf("\n");
        for(std::size_t y=0; y <= (interesting[1]-interesting[0]); ++y)
        {
            int ry = interesting[1]-y;
            gdImageSetPixel(im, x,y, color(result[x*maxfreq + ry]*.7));
            if(cand.find(ry) != cand.end())
            {
                gdImageFilledEllipse(im, x,y, 16,16, 0x000000);
                //gdImageSetPixel(im, x,y, 0x000000);
            }
            if(y%200 == 0) gdImageSetPixel(im, x,y, 0x5555FF);
            if(y%1000 == 0) gdImageSetPixel(im, x,y, 0xFF5555);
        }
    }
    FILE* fp = std::fopen("test.png", "wb");
    gdImagePng(im, fp);
    std::fclose(fp);
}
int main()
{
    std::ifstream f("test.raw", std::ios::binary);
    f.seekg(0, std::ios::end);
    std::vector<char> data(f.tellg());
    f.seekg(0, std::ios::beg);
    f.read(&data[0], data.size());
    Analyze(reinterpret_cast<float*>(&data[0]), data.size() / sizeof(float), 44100);
}
