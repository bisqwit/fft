#pragma once
#include <vector>
#include <complex>

using complex = std::complex<float>;
using cvector = std::vector<complex>;

inline complex swap(complex b) { return complex(b.imag(), b.real()); }
