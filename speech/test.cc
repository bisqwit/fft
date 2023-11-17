#include <algorithm>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdio>
#include "fft_any.hh"
#include "fft_fftw.hh"

static void Analyze(float* data, std::size_t size, std::size_t rate)
{
    cvector sample(size);
    for(std::size_t p=0; p<size; ++p) sample[p] = data[p];
    sample = FFT_fftw(sample);
    for(std::size_t p=0; p<size; ++p) data[p] = std::abs(sample[p]) * 75 / size;
}
int main()
{
    std::ifstream f("test.raw", std::ios::binary);
    f.seekg(0, std::ios::end);
    std::vector<char> data(f.tellg());
    f.seekg(0, std::ios::beg);
    f.read(&data[0], data.size());
    Analyze(reinterpret_cast<float*>(&data[0]), data.size() / sizeof(float), 44100);

    std::ofstream f2("test2.raw", std::ios::binary);
    f2.write(&data[0], data.size());
}
