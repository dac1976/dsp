# DSP (Digital Signal Processing) #
## Introduction ##
As part of my day job I frequently have to make use of various digital signal processing techniques. At work we use Intel's IPP libraries but I decided that for small projects at home or for platforms and operating systems where these libraries aren't available or practical then having a lightweight header only C++14 library would be quite useful.

It is licensed under the terms of LGPL 3.0 and the relevant documentation for this can be found at the top of each source file and in the LICENSE text file.

The code is the work of me (Duncan Crutchley) (<dac1976github@outlook.com>).

Copyright (C) 2018 onwards Duncan Crutchley.

## Requirements ##
Any compatible C++14 compiler can be used. The code is header only and has no external dependencies. It uses just pure C++. It has been tested with Clang, GCC and MSVC2015 and 2017. It has been tested on Linux and Windows 10 64bit OSes. I've also run it successfully on my Samsung Galaxy S9 phone using the C4Droid (based on GCC) compiler available in the Google Play Store.

## Contents ##
The library contains classes and helper functions for the following.

* Mathematical Utilties: discrete convolution, Bessel function, Sinc functions, Sinusoidal function, GCD (greatest common divisor), power of 2 checker.
* Pi: Multi-precision static template functions to generate variants of Pi, e.g. Pi, 2Pi, Pi/2 etc.
* Roots: Multi-precision static template functions to generate variants of Sqrt(2), e.g. Sqrt(2), 2Sqrt, 1/Sqrt(2) etc. 
* Signal Generators: template functions to generate single- and multi-tone sinusoidal signals.
* Window Functions: template functions and functors to evaluate window functions and generate coefficients, e.g. flat-top, Hann, Hamming, rectangle, sinc, Bartlett, Blackman, Kaiser and Lanczos.
* Filters: template functions and functors to compute FIR filter coeffcients, e.g. low-, high-, band-pass and notch filters.
* Resampling: template functors to perform filter-based signal data resampling, including computing up- and down-sample factors.
* Fast Fourier Transforms: template functors to perform various forward and inverse FFT related activities. This also includes fast FFT-based convolution.

## Usage ##
Simply include dsp.hpp and start uing the code. Everything is in the "dsp" namespace.

The code is commented using Doxygen style comments so check the documentation in the /docs/html/ folder and open index.html in your browser of choice.

For usage examples I recommend looking at the unit test code in the file "test_main.cpp".

## Comments ##
I'll fix bugs and improve the code when necessary but make no guarantees on how often this happens. I provide no warranty or support for any issues that are encountered while using it. Although if you are really stuck email me at the provided address and if I have the time I will try to help or fix the issue if it's within my power to do so
