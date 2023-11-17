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

## Graph rendering tool

The program `render.cc` is used to render
speed comparison graphs for the thesis.
This tool was created because LibreOffice Calc
was too sluggish to operate with CSV files
containing tens of thousands of rows and multiple columns.

## Speech recognition test

The program `main.cc` in the `speech/` branch
is used for analyzing vowel formants in speech by FFT.
