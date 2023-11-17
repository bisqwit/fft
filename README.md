# Fast Fourier Transform - implementation in C++

It’s a collection of Fast Fourier Transform algorithms,
featured in Joel Yliluoma’s master’s thesis (“pro gradu”)
in University of Helsinki.

## Algorithms

Supported algorithms:

* Discrete Fourier Transform (naive DFT)
* Cooley-Tukey FFT (radix-2)
* Cooley-Tukey FFT (generic)
* Rader’s FFT
* Bluestein’s FFT
* Combined (heuristic) algorithm
* FFTW (for reference)

## Structure

* `defs.hh`: Definitions for complex types and complex vectors.
* `exp.hh` and `exp.cc`: A module for handling exp-vectors containing $e^{-i2\pi x/N}$ values.
* `factor.hh` and `factor.cc`: A module for integer factorization utilities.
* `dft.hh` and `dft.cc`: The naive DFT. Also works as a fallback for any algorithm. Contains optimized FFTs for sizes ≤ 4.
* `fft_fftw.hh` and `fft_fftw.cc`: FFT through libFFTW3.
* `fft_radix2.hh` and `fft_radix2.cc`: Fast FFT for powers-of-two lengths (Cooley-Tukey radix-2).
* `fft_tukey.hh` and `fft_tukey.cc`: Generic Cooley-Tukey algorithm. Requires that the input length is a composite number.
* `fft_rader.hh` and `fft_rader.cc`: Rader’s FFT. Requires that the input length is a prime number. Performs sub-FFTs using FFT_any for length $N-1$.
* `fft_bluestein.hh` and `fft_bluestein.cc`: Bluestein’s FFT. Performs sub-FFTs using FFT_any for length $m$. Provided are three versions:
* * `FFT_bluestein` selects smallest convolution length $m ≥ 2N-1$ such that $m$ consists of only these factors: 2, 3, 5 or 7.
* * `FFT_bluestein_pow2` selects smallest convolution length $m ≥ 2N-1$ such that $m$ consists of only these factors: 2.
* * `FFT_bluestein_fac2_3` selects smallest convolution length $m ≥ 2N-1$ such that $m$ consists of only these factors: 2 or 3.
* `fft_any.hh` and `fft_any.cc`: Combined FFT. Uses heuristics to decide the fastest method for any input size (excluding FFTW).

Also

* `render.cc`: See below
* `planner.cc`: A tool for generating FFTW wisdom files

## Graph rendering tool

The program `render.cc` is used to render
speed comparison graphs for the thesis.
This tool was created because LibreOffice Calc
was too sluggish to operate with CSV files
containing tens of thousands of rows and multiple columns.

## Speech recognition test

The program `main.cc` in the `speech/` branch
is used for analyzing vowel formants in speech by FFT.
