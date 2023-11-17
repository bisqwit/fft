#include <algorithm>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <gd.h>
#include "fft_fftw.hh"

static const struct Formant
{
    char symbol[8];
    int f1, f2;
} formants[] =
{
    { "i", 240, 2400 },
    { "y", 235, 2100 },
    { "e", 390, 2300 },
    { "ɛ", 610, 1900 },
    { "ø", 370, 1900 },
    { "œ", 585, 1710 },
    { "a", 850, 1610 },
    { "ɶ", 820, 1530 },
    { "ɑ", 750, 940 },
    { "ɒ", 700, 760 },
    { "ʌ", 600, 1170 },
    { "ɔ", 500, 700 },
    { "ɤ", 460, 1310 },
    { "o", 360, 640 },
    { "ɯ", 300, 1390 },
    { "u", 250, 595 },
};
int main()
{
    /* Read input file */
    std::ifstream f("test.raw", std::ios::binary);
    f.seekg(0, std::ios::end);
    std::vector<char> buf(f.tellg());
    f.seekg(0, std::ios::beg);
    f.read(&buf[0], buf.size());
    /* Analyze the data */
    const float* data = reinterpret_cast<float*>(&buf[0]);
    std::size_t size = buf.size() / sizeof(float), rate = 44100;
    double length    = size / double(rate); /* Length in seconds */
    double increment = length/4000, window = 0.013;
    std::size_t maxfreq = 44100, interesting[2] = {0,5000}, x = 0;
    std::size_t columns = std::ceil(length/increment)+1;
    std::vector<double> result(maxfreq * columns);
    for(double start = 0; start < length; start += increment, ++x)
    {
        /* Generate complex sample */
        cvector sample(maxfreq);
        for(std::size_t p=start*rate, a=0; a<maxfreq; ++a,++p) { sample[a] = p < size ? data[p] : 0; }
        /* Apply Gaussian window */
        std::size_t wnd = window*rate;
        double wndsum = 0, C = 48, edge = std::exp(-C/4);
        for(std::size_t a=0; a<maxfreq; ++a)
        {
            double value = (std::exp(-C*std::pow(double(a)/wnd-0.5, 2.)) - edge) / (1-edge);
            sample[a] *= value; wndsum += value;
        }
        // Normalize it
        for(std::size_t a=0; a<maxfreq; ++a) sample[a] /= wndsum;
        /* Transform */
        sample = FFT_fftw(sample);
        /* Convert complex FFT into real power spectrum */
        std::vector<float> real(maxfreq/2);
        for(std::size_t a=0; a<maxfreq/2; ++a) { real[a] = std::abs(sample[a]); }
        /* Convert absolute signal to dB */
        for(std::size_t a=0; a<maxfreq/2; ++a) { real[a] = 10 * std::log10((real[a] + 1e-30) / 4e-10) + 6*std::log2((a+!a)/1000.0); }
        /* Clamp to 60 ... 90 dB range */
        float min = 60, max = 90;
        for(std::size_t a=0; a<maxfreq/2; ++a) { real[a] = std::clamp((real[a] - min) / (max-min), 0.f, 1.f); }
        /* Render into result */
        for(std::size_t a=0; a<maxfreq/2; ++a) { result[x*maxfreq + a] = real[a]; }
    }
    /* Convert the result into a picture */
    gdImagePtr im = gdImageCreateTrueColor(columns, interesting[1]-interesting[0]);
    struct totuus
    {
        unsigned first,last;
        char totuus[8];
        unsigned hits[16];
    } chart[] =
    {
    #define X(n) unsigned(n*columns/250)
        { X(0), X(35),  "-", {} },
        { X(35), X(42), "ʋ", {} },
        { X(42), X(49), "ɑ", {} },
        { X(49), X(57), "e", {} },
        { X(57), X(70), "-", {} },
        { X(70), X(76), "k", {} },
        { X(76), X(82), "ɑ", {} },
        { X(82), X(99), "m", {} },
        { X(99), X(127), "ɑ", {} },
        { X(127), X(143), "p", {} },
        { X(143), X(149), "ɑ", {} },
        { X(149), X(162), "l", {} },
        { X(162), X(172), "o", {} },
        { X(172), X(186), "n", {} },
        { X(186), X(192), "r", {} },
        { X(192), X(213), "ɑ", {} },
        { X(213), X(220), "k", {} },
        { X(220), X(231), "e", {} },
        { X(231), X(244), "n", {} },
        { X(244), X(250), "e", {} },
    #undef X
    };
    for(std::size_t x=0; x<columns; ++x)
    {
        double best = 0; unsigned i = 0, besti = 0;
        const char* name = nullptr;
        for(const auto& f: formants)
        {
            double power1 = result[x*maxfreq + f.f1];
            double power2 = result[x*maxfreq + f.f2];
            double score = power1*power1 + power2*power2;
            if(score > best) { best = score; name = f.symbol; besti = i; }
            ++i;
        }
        unsigned p=0;
        while(chart[p].last <= x) ++p;
        ++chart[p].hits[besti];

        //std::printf("x=%zu: symbol=%s\n", x, name);
        for(std::size_t y=0; y < (interesting[1]-interesting[0]); ++y)
        {
            float f = std::pow(1 - result[x*maxfreq + (interesting[1]-y)], 1.2f);
            //if(x == 800) std::printf("\t%zu\t%.8g\n", interesting[1]-y, f);
            int b = f*768, g = b-256, r = g-256;
            unsigned color = (std::clamp(r,0,255)<<16)+(std::clamp(g,0,255)<<8)+(std::clamp(b,0,255)<<0);
            gdImageSetPixel(im, x,y, color);
            if(y%200 == 0) gdImageSetPixel(im, x,y, 0x5555FF);
            if(y%1000 == 0) gdImageSetPixel(im, x,y, 0xFF5555);
        }
    }
    for(auto& d: chart)
    {
        std::string vowel = d.totuus;
        if(vowel == "ɑ") vowel = "A";
        if(vowel == "ʋ") vowel = "V";
        std::printf("$%.2f$--$%.2f$ & \\textipa{%s}",
            d.first*length/columns, d.last*length/columns,
            vowel.c_str());
        for(unsigned n=0; n<16; ++n)
            printf(" & $%.0f$\\,\\%%", d.hits[n] * 100. / (d.last-d.first));
        std::printf(" \\\\\n");
    }
    /* Write out the picture */
    FILE* fp = std::fopen("test.png", "wb");
    gdImagePng(im, fp);
    std::fclose(fp);
}
/*
Ground truth:
250 frames
    0-34 = silence
    35-42 = v
    42-47 = a
    48-52 = e
    53-57 = y
    70-75 = k
    76-82 = a
    83-98 = m
    99-126 = a
    127-142 = p
    143-148 = a
    149-161 = l
    162-171 = o
    172-185 = n
    186-191 = r
    192-201 = a
    213-219 = k
    220-230 = e
    231-243 = n
    244-250 = e
*/