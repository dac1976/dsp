// This file is part of DSP library containing useful reusable
// signal processing utility classes.
//
// Copyright (C) 2018 Duncan Crutchley
// Contact <dac1976github@outlook.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License and GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU General Public License
// and GNU Lesser General Public License along with this program. If
// not, see <http://www.gnu.org/licenses/>.
#include <iostream>
#include <fstream>
#include <numeric>
#include "dsp.hpp"
#include "test_timer.hpp"

/*! \brief dsp namespace */
namespace dsp
{

/*! \brief Struct to hold stats data. */
struct Statistics
{
	/*! \brief Variable to hold mean. */
    double mean{0.};
    /*! \brief Variable to hold max. */
    double max{0.};
    /*! \brief Variable to hold min. */
    double min{0.};
    /*! \brief Variable to hold stdandard deviation. */
    double stdDev{0.};
};

/*!
 * \brief Function to compute statistics.
 * \return A pair reporting the number of tests failed and the number of tests run.
 */
template <typename Iter> Statistics ComputeStats(Iter first, Iter last)
{
    auto count = std::distance(first, last);

    if (count <= 0)
    {
        return Statistics{};
    }

    Statistics stats;
    auto       minMaxIters = std::minmax_element(first, last);
    stats.min              = *minMaxIters.first;
    stats.max              = *minMaxIters.second;
    stats.mean             = std::accumulate(first, last, 0.0);
    stats.mean /= static_cast<double>(count);

    for (auto i = first; i != last; ++i)
    {
        auto arg = *i - stats.mean;
        stats.stdDev += arg * arg;
    }

    stats.stdDev /= static_cast<double>(count - 1);
    stats.stdDev = sqrt(stats.stdDev);

    return stats;
}

/*!
 * \brief Unit test function for Window functions.
 * \return A pair reporting the number of tests failed and the number of tests run.
 */
std::pair<int, int> TestWindowFunction()
{
    std::pair<int, int> result{0, 17};
    std::cout << "Testing dsp::WindowFunction class..." << std::endl;
    std::cout << "[Test 1 - FlatTop1Generator]" << std::endl;

    test::Timer       t;
    FlatTop1Generator ft1;
    window_fn_d       wf1(ft1, 1025, true);
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    if ((std::abs(wf1.CoherentGain() - 1.0) < 1.e-2) && (std::abs(wf1.PowerGain() - 3.77) < 1.e-2))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        std::cout << "Test passed? = false" << std::endl;
        ++result.first;
    }

    std::cout << "[Test 2 - FlatTop2Generator]" << std::endl;

    t.Reset();
    FlatTop2Generator ft2;
    window_fn_d       wf2(ft2, 1025, true);
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    if ((std::abs(wf2.CoherentGain() - 0.28) < 1.e-2) &&
        (std::abs(wf2.PowerGain() - 0.234) < 1.e-2))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        std::cout << "Test passed? = false" << std::endl;
        ++result.first;
    }

    std::cout << "[Test 3 - FlatTop3Generator]" << std::endl;

    t.Reset();
    FlatTop3Generator ft3;
    window_fn_d       wf3(ft3, 1025, true);
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    if ((std::abs(wf3.CoherentGain() - 0.22) < 1.e-2) &&
        (std::abs(wf3.PowerGain() - 0.175) < 1.e-2))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        std::cout << "Test passed? = false" << std::endl;
        ++result.first;
    }

    std::cout << "[Test 4 - FlatTop4Generator]" << std::endl;

    t.Reset();
    FlatTop4Generator ft4;
    window_fn_d       wf4(ft4, 1025, true);
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    if ((std::abs(wf4.CoherentGain() - 0.999) < 1.e-2) &&
        (std::abs(wf4.PowerGain() - 3.42) < 1.e-2))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        std::cout << "Test passed? = false" << std::endl;
        ++result.first;
    }

    std::cout << "[Test 5 - FlatTop5Generator]" << std::endl;

    t.Reset();
    FlatTop5Generator ft5;
    window_fn_d       wf5(ft5, 1025, true);
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    if ((std::abs(wf5.CoherentGain() - 1.0) < 1.e-2) && (std::abs(wf5.PowerGain() - 3.46) < 1.e-2))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        std::cout << "Test passed? = false" << std::endl;
        ++result.first;
    }

    std::cout << "[Test 6 - FlatTop6Generator]" << std::endl;

    t.Reset();
    FlatTop6Generator ft6;
    window_fn_d       wf6(ft6, 1025, true);
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    if ((std::abs(wf6.CoherentGain() - 1.0) < 1.e-2) && (std::abs(wf6.PowerGain() - 3.85) < 1.e-2))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        std::cout << "Test passed? = false" << std::endl;
        ++result.first;
    }

    std::cout << "[Test 7 - FlatTop7Generator]" << std::endl;

    t.Reset();
    FlatTop7Generator ft7;
    window_fn_d       wf7(ft7, 1025, true);
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    if ((std::abs(wf7.CoherentGain() - 0.19) < 1.e-2) &&
        (std::abs(wf7.PowerGain() - 0.154) < 1.e-2))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        std::cout << "Test passed? = false" << std::endl;
        ++result.first;
    }

    std::cout << "[Test 8 - HannGenerator]" << std::endl;

    t.Reset();
    HannGenerator hann;
    window_fn_d   wf8(hann, 1025, true);
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    if ((std::abs(wf8.CoherentGain() - 0.5) < 1.e-2) && (std::abs(wf8.PowerGain() - 0.375) < 1.e-2))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        std::cout << "Test passed? = false" << std::endl;
        ++result.first;
    }

    std::cout << "[Test 9 - HammingGenerator]" << std::endl;

    t.Reset();
    HammingGenerator ham;
    window_fn_d      wf9(ham, 1025, true);
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    if ((std::abs(wf9.CoherentGain() - 0.54) < 1.e-2) &&
        (std::abs(wf9.PowerGain() - 0.397) < 1.e-2))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        std::cout << "Test passed? = false" << std::endl;
        ++result.first;
    }

    std::cout << "[Test 10 - RectangleGenerator]" << std::endl;

    t.Reset();
    RectangleGenerator rect;
    window_fn_d        wf10(rect, 1025, true);
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    if ((std::abs(wf10.CoherentGain() - 1.0) < 1.e-2) && (std::abs(wf10.PowerGain() - 1.0) < 1.e-2))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        std::cout << "Test passed? = false" << std::endl;
        ++result.first;
    }

    std::cout << "[Test 11 - BartlettGenerator]" << std::endl;

    t.Reset();
    BartlettGenerator bart;
    window_fn_d       wf11(bart, 1025, true);
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    if ((std::abs(wf11.CoherentGain() - 0.5) < 1.e-2) &&
        (std::abs(wf11.PowerGain() - 0.333) < 1.e-2))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        std::cout << "Test passed? = false" << std::endl;
        ++result.first;
    }

    std::cout << "[Test 12 - ExactBlackmanGenerator]" << std::endl;

    t.Reset();
    ExactBlackmanGenerator exBlack;
    window_fn_d            wf12(exBlack, 1025, true);
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    if ((std::abs(wf12.CoherentGain() - 0.42) < 1.e-2) &&
        (std::abs(wf12.PowerGain() - 0.309) < 1.e-2))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        std::cout << "Test passed? = false" << std::endl;
        ++result.first;
    }

    std::cout << "[Test 13 - BlackmanGenerator]" << std::endl;

    t.Reset();
    BlackmanGenerator black;
    window_fn_d       wf13(black, 1025, true);
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    if ((std::abs(wf13.CoherentGain() - 0.42) < 1.e-2) &&
        (std::abs(wf13.PowerGain() - 0.305) < 1.e-2))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        std::cout << "Test passed? = false" << std::endl;
        ++result.first;
    }

    std::cout << "[Test 14 - KaiserGenerator beta = 2*Pi]" << std::endl;

    t.Reset();
    KaiserGenerator kaiser1(6.283185307);
    window_fn_d     wf14(kaiser1, 1025, true);
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    if ((std::abs(wf14.CoherentGain() - 0.49) < 1.e-2) &&
        (std::abs(wf14.PowerGain() - 0.359) < 1.e-2))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        std::cout << "Test passed? = false" << std::endl;
        ++result.first;
    }

    std::cout << "[Test 15 - KaiserGenerator beta = 3*Pi]" << std::endl;

    t.Reset();
    KaiserGenerator kaiser2(9.424777961);
    window_fn_d     wf15(kaiser2, 1025, true);
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    if ((std::abs(wf15.CoherentGain() - 0.404) < 1.e-2) &&
        (std::abs(wf15.PowerGain() - 0.292) < 1.e-2))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        std::cout << "Test passed? = false" << std::endl;
        ++result.first;
    }

    std::cout << "[Test 16 - LanczosGenerator]" << std::endl;

    t.Reset();
    LanczosGenerator lanczos;
    window_fn_d      wf16(lanczos, 1025, true);
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    if ((std::abs(wf16.CoherentGain() - 0.59) < 1.e-2) &&
        (std::abs(wf16.PowerGain() - 0.452) < 1.e-2))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        std::cout << "Test passed? = false" << std::endl;
        ++result.first;
    }

    std::cout << "[Test 17 - Apply Hann Window]" << std::endl;
    std::vector<double> data(1024, 1.);

    t.Reset();
    wf8(data.begin(), data.end(), data.begin());
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    if (data == wf8.Coefficients())
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        std::cout << "Test passed? = false" << std::endl;
        ++result.first;
    }

    return result;
}

/*!
 * \brief Unit test function for Convolve function.
 * \return A pair reporting the number of tests failed and the number of tests run.
 */
std::pair<int, int> TestConvolve()
{
    std::pair<int, int> result{0, 4};
    std::cout << "Testing dsp::Convolve function..." << std::endl;
    std::cout << "[Test 1]" << std::endl;
    std::cout << "input 1 = {1, 1, 1, 1, 1, 1}" << std::endl;
    std::cout << "input 2 = {1, 1, 1, 1, 1, 1}" << std::endl;
    std::cout << "expected output = {1, 2, 3, 4, 5, 6, 5, 4, 3, 2, 1}" << std::endl;

    std::vector<double> v1{1., 1., 1., 1., 1., 1.};
    std::vector<double> v2{1., 1., 1., 1., 1., 1.};
    std::vector<double> v3(v1.size() + v2.size() - 1);
    std::vector<double> output{1., 2., 3., 4., 5., 6., 5., 4., 3., 2., 1.};

    test::Timer t;
    Convolve(v1.begin(), v1.end(), v2.begin(), v2.end(), v3.begin());
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    int numCommas = 10;
    std::cout << "actual output = {";

    for (auto v : v3)
    {
        std::cout << v;

        if (numCommas > 0)
        {
            std::cout << ", ";
        }

        --numCommas;
    }

    std::cout << "}" << std::endl;
    auto pass = output == v3;
    std::cout << "Test passed? = " << std::boolalpha << pass << std::endl;

    if (!pass)
    {
        ++result.first;
    }

    std::cout << "[Test 2]" << std::endl;
    std::cout << "input 1 = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}" << std::endl;
    std::cout << "input 2 = {1, 1, 1, 1, 1, 1}" << std::endl;
    std::cout << "expected output = {1, 3, 6, 9, 12, 15, 18, 21, 24, 27, 19, 10}" << std::endl;

    v1 = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10.};
    v2 = {1., 1., 1.};
    v3.resize(v1.size() + v2.size() - 1);
    output = {1., 3., 6., 9., 12., 15., 18., 21., 24., 27., 19., 10.};

    t.Reset();
    Convolve(v1.begin(), v1.end(), v2.begin(), v2.end(), v3.begin());
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    numCommas = 11;
    std::cout << "actual output = {";

    for (auto v : v3)
    {
        std::cout << v;

        if (numCommas > 0)
        {
            std::cout << ", ";
        }

        --numCommas;
    }

    std::cout << "}" << std::endl;
    pass = output == v3;
    std::cout << "Test passed? = " << std::boolalpha << pass << std::endl;

    if (!pass)
    {
        ++result.first;
    }

    std::cout << "[Test 3]" << std::endl;
    std::cout << "input 1 = {1, 1, 1, 1, 1, 1}" << std::endl;
    std::cout << "input 2 = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}" << std::endl;
    std::cout << "expected output = {1, 3, 6, 9, 12, 15, 18, 21, 24, 27, 19, 10}" << std::endl;

    t.Reset();
    Convolve(v2.begin(), v2.end(), v1.begin(), v1.end(), v3.begin());
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    numCommas = 11;
    std::cout << "actual output = {";

    for (auto v : v3)
    {
        std::cout << v;

        if (numCommas > 0)
        {
            std::cout << ", ";
        }

        --numCommas;
    }

    std::cout << "}" << std::endl;
    pass = output == v3;
    std::cout << "Test passed? = " << std::boolalpha << pass << std::endl;

    if (!pass)
    {
        ++result.first;
    }

    std::cout << "[Test 4 - Large vectors]" << std::endl;
    std::vector<double> big1(1001);
    std::vector<double> big2(46500);
    std::vector<double> bigResult(big1.size() + big2.size() - 1);
    std::vector<double> bigResult2(big1.size() + big2.size() - 1);

    std::generate(big1.begin(), big1.end(), []() {
        static int i = 1;
        return i++;
    });

    std::generate(big2.begin(), big2.end(), []() {
        static int i = 1;
        return i++;
    });

    double elapsed = 0.0;
    int    numRuns = 10;
    std::cout << "[Perform discrete convolution, iterations = " << numRuns << "]" << std::endl;

    for (int i = 0; i < numRuns; ++i)
    {
        t.Reset();
        Convolve(big2.begin(), big2.end(), big1.begin(), big1.end(), bigResult.begin());
        elapsed += t.Elapsed();
    }

    std::cout << "\tDuration (mean) " << elapsed / numRuns
              << "s , number of iterations: " << numRuns << std::endl;

    std::cout << "[Create FFT convolution object]" << std::endl;
    t.Reset();
    fft_convolve_d fftConvolve(big1.size(), big2.size());
    std::cout << "\tDuration" << t.Elapsed() << "s" << std::endl;

    elapsed = 0.0;
    numRuns = 10;
    std::cout << "[Perform FFT convolution, iterations = " << numRuns << "]" << std::endl;

    for (int i = 0; i < numRuns; ++i)
    {
        t.Reset();
        fftConvolve(big2.begin(), big2.end(), big1.begin(), big1.end(), bigResult2.begin());
        elapsed += t.Elapsed();
    }

    std::cout << "\tDuration (mean) " << elapsed / numRuns
              << "s , number of iterations: " << numRuns << std::endl;

    auto stats1 = ComputeStats(bigResult.begin(), bigResult.end());
    auto stats2 = ComputeStats(bigResult2.begin(), bigResult2.end());

    if ((std::abs(stats1.min - stats2.min) < 1.e-1) &&
        (std::abs(stats1.max - stats2.max) < 1.e-1) &&
        (std::abs(stats1.mean - stats2.mean) < 1.e-1) &&
        (std::abs(stats1.stdDev - stats2.stdDev) < 1.e-1))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }

    return result;
}

/*!
 * \brief Unit test function for Bessel function.
 * \return A pair reporting the number of tests failed and the number of tests run.
 */
std::pair<int, int> TestBessel()
{
    std::pair<int, int> result{0, 3};
    std::cout << "Testing dsp::Bessel function..." << std::endl;
    std::cout << "[Test 1]" << std::endl;
    std::cout << "input = 0" << std::endl;
    std::cout << "expected output = 1" << std::endl;

    test::Timer t;
    auto        value = Bessel(0.0);
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    std::cout << "actual output = " << value << std::endl;

    if (std::abs(value - 1) > 1.e-9)
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }
    else
    {
        std::cout << "Test passed? = true" << std::endl;
    }

    std::cout << "[Test 2]" << std::endl;
    std::cout << "input = 2" << std::endl;

    double expected = 2.2795853023359909;
    std::cout << "expected output = " << expected << std::endl;

    t.Reset();
    value = Bessel(2.0);
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    std::cout << "actual output = " << value << std::endl;

    if (std::abs(value - expected) > 1.e-9)
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }
    else
    {
        std::cout << "Test passed? = true" << std::endl;
    }

    std::cout << "[Test 3]" << std::endl;
    std::cout << "input = 3" << std::endl;

    expected = 4.8807925856077325;
    std::cout << "expected output = " << expected << std::endl;

    t.Reset();
    value = Bessel(3.0);
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    std::cout << "actual output = " << value << std::endl;

    if (std::abs(value - expected) > 1.e-9)
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }
    else
    {
        std::cout << "Test passed? = true" << std::endl;
    }

    return result;
}

/*!
 * \brief Unit test function for Sinc and SincNorm functions.
 * \return A pair reporting the number of tests failed and the number of tests run.
 */
std::pair<int, int> TestSinc()
{
    std::pair<int, int> result{0, 4};
    std::cout << "Testing dsp::Sinc function..." << std::endl;
    std::cout << "[Test 1]" << std::endl;
    std::cout << "input = 0" << std::endl;
    std::cout << "expected output = 1" << std::endl;

    test::Timer t;
    auto        value = Sinc(0.0);
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    std::cout << "actual output = " << value << std::endl;

    if (std::abs(value - 1) > 1.e-9)
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }
    else
    {
        std::cout << "Test passed? = true" << std::endl;
    }

    std::cout << "[Test 2]" << std::endl;
    std::cout << "input = 1" << std::endl;

    auto expected = sin(1);
    std::cout << "expected output = " << expected << std::endl;

    t.Reset();
    value = Sinc(1.);
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    std::cout << "actual output = " << value << std::endl;

    if (std::abs(value - expected) > 1.e-9)
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }
    else
    {
        std::cout << "Test passed? = true" << std::endl;
    }

    std::cout << std::endl << "Testing dsp::SincNorm function..." << std::endl;
    std::cout << "[Test 1]" << std::endl;
    std::cout << "input = 0" << std::endl;
    std::cout << "expected output = 1" << std::endl;

    t.Reset();
    value = SincNorm(0.0);
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    std::cout << "actual output = " << value << std::endl;

    if (std::abs(value - 1) > 1.e-9)
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }
    else
    {
        std::cout << "Test passed? = true" << std::endl;
    }

    std::cout << "[Test 2]" << std::endl;
    std::cout << "input = 1" << std::endl;

    expected = sin(Pi<double>()) / Pi<double>();
    std::cout << "expected output = " << expected << std::endl;

    t.Reset();
    value = SincNorm(1.);
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    std::cout << "actual output = " << value << std::endl;

    if (std::abs(value - expected) > 1.e-9)
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }
    else
    {
        std::cout << "Test passed? = true" << std::endl;
    }

    return result;
}

/*!
 * \brief Unit test function for Sine function.
 * \return A pair reporting the number of tests failed and the number of tests run.
 */
std::pair<int, int> TestSine()
{
    std::pair<int, int> result(0, 5);
    std::cout << "Testing dsp::Sine function..." << std::endl;
    std::cout << "[Test 1]" << std::endl;
    std::cout << "input = {5, 0, 1, 0, 0}" << std::endl;
    std::cout << "expected output = 0" << std::endl;

    test::Timer t;
    auto        value = dsp::Sine(5., 0., 1., 0., 0.);
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    std::cout << "actual output = " << value << std::endl;

    if (std::abs(value) > 1.e-9)
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }
    else
    {
        std::cout << "Test passed? = true" << std::endl;
    }

    std::cout << "[Test 2]" << std::endl;
    std::cout << "input = {5, 0.25, 1, 0, 0}" << std::endl;
    std::cout << "expected output = 5" << std::endl;

    t.Reset();
    value = dsp::Sine(5., 0.25, 1., 0., 0.);
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    std::cout << "actual output = " << value << std::endl;

    if (std::abs(value - 5.) > 1.e-9)
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }
    else
    {
        std::cout << "Test passed? = true" << std::endl;
    }

    std::cout << "[Test 3]" << std::endl;
    std::cout << "input = {5, 0.25, 1, 0, 0}" << std::endl;
    std::cout << "expected output = -5" << std::endl;

    t.Reset();
    value = dsp::Sine(5., 0.75, 1., 0., 0.);
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    std::cout << "actual output = " << value << std::endl;

    if (std::abs(value - -5.) > 1.e-9)
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }
    else
    {
        std::cout << "Test passed? = true" << std::endl;
    }

    std::cout << "[Test 4]" << std::endl;
    std::cout << "input = {5, 0.25, 1, 0, 5}" << std::endl;
    std::cout << "expected output = 0" << std::endl;

    t.Reset();
    value = dsp::Sine(5., 0.75, 1., 0., 5.);
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    std::cout << "actual output = " << value << std::endl;

    if (std::abs(value) > 1.e-9)
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }
    else
    {
        std::cout << "Test passed? = true" << std::endl;
    }

    std::cout << "[Test 5]" << std::endl;
    std::cout << "input = {5, 0, 1, Pi/2, 0}" << std::endl;
    std::cout << "expected output = 5" << std::endl;

    t.Reset();
    value = dsp::Sine(5., 0., 1., dsp::HalfPi<double>(), 0.);
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    std::cout << "actual output = " << value << std::endl;

    if (std::abs(value - 5) > 1.e-9)
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }
    else
    {
        std::cout << "Test passed? = true" << std::endl;
    }

    return result;
}

/*!
 * \brief Unit test function for filter functions using fast convolution.
 * \return A pair reporting the number of tests failed and the number of tests run.
 */
std::pair<int, int> TestFilters1(bool logToFile)
{
    std::pair<int, int> result(0, 4);

    std::cout << "Testing dsp::FilterHolder class (fast convolution)..." << std::endl;

    ToneParams t1{10., 50., 0., 0.};
    ToneParams t2{5., 150., 0., 0.};
    ToneParams t3{2., 500., 0., 0.};

    auto multiToneSignal1 = MultiTone<double>({t1, t2, t3}, 2000., 2000);

    if (logToFile)
    {
        std::ofstream ofs("multiToneSignal1.csv");

        for (auto s : multiToneSignal1)
        {
            ofs << s << std::endl;
        }
    }

    auto multiToneSignal2 = MultiTone<double>({t2, t3}, 2000., 2000);

    if (logToFile)
    {
        std::ofstream ofs("multiToneSignal2.csv");

        for (auto s : multiToneSignal2)
        {
            ofs << s << std::endl;
        }
    }

    auto toneSignal1 = Tone<double>(t1, 2000., 2000);

    if (logToFile)
    {
        std::ofstream ofs("toneSignal1.csv");

        for (auto s : toneSignal1)
        {
            ofs << s << std::endl;
        }
    }

    auto toneSignal2 = Tone<double>(t2, 2000., 2000);

    if (logToFile)
    {
        std::ofstream ofs("toneSignal2.csv");

        for (auto s : toneSignal2)
        {
            ofs << s << std::endl;
        }
    }

    auto toneSignal3 = Tone<double>(t3, 2000., 2000);

    if (logToFile)
    {
        std::ofstream ofs("toneSignal3.csv");

        for (auto s : toneSignal3)
        {
            ofs << s << std::endl;
        }
    }

    std::cout << "[Test 1 - FirLowPassFilter]" << std::endl;

    KaiserGenerator kaiser(10.);

    test::Timer t;
    auto        filterCoeffs = FirLowPassFilter(451, 100., 2000., kaiser);
    std::cout << "\t(FirLowPassFilter) Duration " << t.Elapsed() << "s" << std::endl;

    if (logToFile)
    {
        std::ofstream ofs("lowPass.csv");

        for (auto fc : filterCoeffs)
        {
            ofs << fc << std::endl;
        }
    }

    t.Reset();
    filter_hldr_d lpFilter(2000, filterCoeffs, true);
    std::cout << "\t(filter_hldr_d) Duration " << t.Elapsed() << "s" << std::endl;

    std::vector<double> filteredSignalLp(multiToneSignal1.size());

    t.Reset();
    lpFilter(multiToneSignal1.begin(), multiToneSignal1.end(), filteredSignalLp.begin(), true);
    std::cout << "\t(Apply filter) Duration " << t.Elapsed() << "s" << std::endl;

    auto stats     = ComputeStats(std::next(filteredSignalLp.begin(), 100),
                              std::next(filteredSignalLp.begin(), 1100));
    auto statsComp = ComputeStats(toneSignal1.begin(), toneSignal1.end());

    if ((std::abs(stats.min - statsComp.min) < 1.e-1) &&
        (std::abs(stats.max - statsComp.max) < 1.e-1) &&
        (std::abs(stats.mean - statsComp.mean) < 1.e-2) &&
        (std::abs(stats.stdDev - statsComp.stdDev) < 1.e-2))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }

    if (logToFile)
    {
        std::ofstream ofs("filteredLowPassSignal.csv");

        for (auto fs : filteredSignalLp)
        {
            ofs << fs << std::endl;
        }
    }

    std::cout << "[Test 2 - FirHighPassFilter]" << std::endl;

    t.Reset();
    filterCoeffs = FirHighPassFilter(451, 400., 2000., kaiser);
    std::cout << "\t(FirHighPassFilter) Duration " << t.Elapsed() << "s" << std::endl;

    if (logToFile)
    {
        std::ofstream ofs("highPass.csv");

        for (auto fc : filterCoeffs)
        {
            ofs << fc << std::endl;
        }
    }

    t.Reset();
    filter_hldr_d hpFilter(2000, filterCoeffs, true);
    std::cout << "\t(filter_hldr_d) Duration " << t.Elapsed() << "s" << std::endl;

    std::vector<double> filteredSignalHp(multiToneSignal1.size());

    t.Reset();
    hpFilter(multiToneSignal1.begin(), multiToneSignal1.end(), filteredSignalHp.begin(), true);
    std::cout << "\t(Apply filter) Duration " << t.Elapsed() << "s" << std::endl;

    stats     = ComputeStats(std::next(filteredSignalHp.begin(), 100),
                         std::next(filteredSignalHp.begin(), 1100));
    statsComp = ComputeStats(toneSignal3.begin(), toneSignal3.end());

    if ((std::abs(stats.min - statsComp.min) < 1.e-1) &&
        (std::abs(stats.max - statsComp.max) < 1.e-1) &&
        (std::abs(stats.mean - statsComp.mean) < 1.e-2) &&
        (std::abs(stats.stdDev - statsComp.stdDev) < 1.e-2))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }

    if (logToFile)
    {
        std::ofstream ofs("filteredHighPassSignal.csv");

        for (auto fs : filteredSignalHp)
        {
            ofs << fs << std::endl;
        }
    }

    std::cout << "[Test 3 - FirBandPassFilter]" << std::endl;

    t.Reset();
    filterCoeffs = FirBandPassFilter(451, 150., 100., 2000., kaiser);
    std::cout << "\t(FirBandPassFilter) Duration " << t.Elapsed() << "s" << std::endl;

    if (logToFile)
    {
        std::ofstream ofs("bandPass.csv");

        for (auto fc : filterCoeffs)
        {
            ofs << fc << std::endl;
        }
    }

    t.Reset();
    filter_hldr_d bpFilter(2000, filterCoeffs, true);
    std::cout << "\t(filter_hldr_d) Duration " << t.Elapsed() << "s" << std::endl;

    std::vector<double> filteredSignalBp(multiToneSignal1.size());

    t.Reset();
    bpFilter(multiToneSignal1.begin(), multiToneSignal1.end(), filteredSignalBp.begin(), true);
    std::cout << "\t(Apply filter) Duration " << t.Elapsed() << "s" << std::endl;

    stats     = ComputeStats(std::next(filteredSignalBp.begin(), 100),
                         std::next(filteredSignalBp.begin(), 1100));
    statsComp = ComputeStats(toneSignal2.begin(), toneSignal2.end());

    if ((std::abs(stats.min - statsComp.min) < 1.e-1) &&
        (std::abs(stats.max - statsComp.max) < 1.e-1) &&
        (std::abs(stats.mean - statsComp.mean) < 1.e-2) &&
        (std::abs(stats.stdDev - statsComp.stdDev) < 1.e-2))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }

    if (logToFile)
    {
        std::ofstream ofs("filtereBandPassSignal.csv");

        for (auto fs : filteredSignalHp)
        {
            ofs << fs << std::endl;
        }
    }

    std::cout << "[Test 4 - FirNotchFilter]" << std::endl;

    t.Reset();
    filterCoeffs = FirNotchFilter(451, 150., 10., 2000., kaiser);
    std::cout << "\t(FirNotchFilter) Duration " << t.Elapsed() << "s" << std::endl;

    if (logToFile)
    {
        std::ofstream ofs("notch.csv");

        for (auto fc : filterCoeffs)
        {
            ofs << fc << std::endl;
        }
    }

    t.Reset();
    filter_hldr_d nFilter(2000, filterCoeffs, true);
    std::cout << "\t(filter_hldr_d) Duration " << t.Elapsed() << "s" << std::endl;

    std::vector<double> filteredSignalN(multiToneSignal1.size());

    t.Reset();
    nFilter(multiToneSignal2.begin(), multiToneSignal2.end(), filteredSignalN.begin(), true);
    std::cout << "\t(Apply filter) Duration " << t.Elapsed() << "s" << std::endl;

    stats     = ComputeStats(std::next(filteredSignalN.begin(), 100),
                         std::next(filteredSignalN.begin(), 1100));
    statsComp = ComputeStats(toneSignal3.begin(), toneSignal3.end());

    if ((std::abs(stats.min - statsComp.min) < 1.e-1) &&
        (std::abs(stats.max - statsComp.max) < 1.e-1) &&
        (std::abs(stats.mean - statsComp.mean) < 1.e-2) &&
        (std::abs(stats.stdDev - statsComp.stdDev) < 1.e-2))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }

    if (logToFile)
    {
        std::ofstream ofs("filtereNotchSignal.csv");

        for (auto fs : filteredSignalN)
        {
            ofs << fs << std::endl;
        }
    }

    return result;
}

/*!
 * \brief Unit test function for filter functions using slow convolution.
 * \return A pair reporting the number of tests failed and the number of tests run.
 */
std::pair<int, int> TestFilters2(bool logToFile)
{
    std::pair<int, int> result(0, 4);

    std::cout << "Testing dsp::FilterHolder class (slow convolution)..." << std::endl;

    ToneParams t1{10., 50., 0., 0.};
    ToneParams t2{5., 150., 0., 0.};
    ToneParams t3{2., 500., 0., 0.};

    auto multiToneSignal1 = MultiTone<double>({t1, t2, t3}, 2000., 2000);

    if (logToFile)
    {
        std::ofstream ofs("multiToneSignal1.csv");

        for (auto s : multiToneSignal1)
        {
            ofs << s << std::endl;
        }
    }

    auto multiToneSignal2 = MultiTone<double>({t2, t3}, 2000., 2000);

    if (logToFile)
    {
        std::ofstream ofs("multiToneSignal2.csv");

        for (auto s : multiToneSignal2)
        {
            ofs << s << std::endl;
        }
    }

    auto toneSignal1 = Tone<double>(t1, 2000., 2000);

    if (logToFile)
    {
        std::ofstream ofs("toneSignal1.csv");

        for (auto s : toneSignal1)
        {
            ofs << s << std::endl;
        }
    }

    auto toneSignal2 = Tone<double>(t2, 2000., 2000);

    if (logToFile)
    {
        std::ofstream ofs("toneSignal2.csv");

        for (auto s : toneSignal2)
        {
            ofs << s << std::endl;
        }
    }

    auto toneSignal3 = Tone<double>(t3, 2000., 2000);

    if (logToFile)
    {
        std::ofstream ofs("toneSignal3.csv");

        for (auto s : toneSignal3)
        {
            ofs << s << std::endl;
        }
    }

    std::cout << "[Test 1 - FirLowPassFilter]" << std::endl;

    KaiserGenerator kaiser(10.);

    test::Timer t;
    auto        filterCoeffs = FirLowPassFilter(451, 100., 2000., kaiser);
    std::cout << "\t(FirLowPassFilter) Duration " << t.Elapsed() << "s" << std::endl;

    if (logToFile)
    {
        std::ofstream ofs("lowPass.csv");

        for (auto fc : filterCoeffs)
        {
            ofs << fc << std::endl;
        }
    }

    t.Reset();
    filter_hldr_d lpFilter(2000, filterCoeffs, false);
    std::cout << "\t(filter_hldr_d) Duration " << t.Elapsed() << "s" << std::endl;

    std::vector<double> filteredSignalLp(multiToneSignal1.size());

    t.Reset();
    lpFilter(multiToneSignal1.begin(), multiToneSignal1.end(), filteredSignalLp.begin(), true);
    std::cout << "\t(Apply filter) Duration " << t.Elapsed() << "s" << std::endl;

    auto stats     = ComputeStats(std::next(filteredSignalLp.begin(), 100),
                              std::next(filteredSignalLp.begin(), 1100));
    auto statsComp = ComputeStats(toneSignal1.begin(), toneSignal1.end());

    if ((std::abs(stats.min - statsComp.min) < 1.e-1) &&
        (std::abs(stats.max - statsComp.max) < 1.e-1) &&
        (std::abs(stats.mean - statsComp.mean) < 1.e-2) &&
        (std::abs(stats.stdDev - statsComp.stdDev) < 1.e-2))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }

    if (logToFile)
    {
        std::ofstream ofs("filteredLowPassSignal.csv");

        for (auto fs : filteredSignalLp)
        {
            ofs << fs << std::endl;
        }
    }

    std::cout << "[Test 2 - FirHighPassFilter]" << std::endl;

    t.Reset();
    filterCoeffs = FirHighPassFilter(451, 400., 2000., kaiser);
    std::cout << "\t(FirHighPassFilter) Duration " << t.Elapsed() << "s" << std::endl;

    if (logToFile)
    {
        std::ofstream ofs("highPass.csv");

        for (auto fc : filterCoeffs)
        {
            ofs << fc << std::endl;
        }
    }

    t.Reset();
    filter_hldr_d hpFilter(2000, filterCoeffs, false);
    std::cout << "\t(filter_hldr_d) Duration " << t.Elapsed() << "s" << std::endl;

    std::vector<double> filteredSignalHp(multiToneSignal1.size());

    t.Reset();
    hpFilter(multiToneSignal1.begin(), multiToneSignal1.end(), filteredSignalHp.begin(), true);
    std::cout << "\t(Apply filter) Duration " << t.Elapsed() << "s" << std::endl;

    stats     = ComputeStats(std::next(filteredSignalHp.begin(), 100),
                         std::next(filteredSignalHp.begin(), 1100));
    statsComp = ComputeStats(toneSignal3.begin(), toneSignal3.end());

    if ((std::abs(stats.min - statsComp.min) < 1.e-1) &&
        (std::abs(stats.max - statsComp.max) < 1.e-1) &&
        (std::abs(stats.mean - statsComp.mean) < 1.e-2) &&
        (std::abs(stats.stdDev - statsComp.stdDev) < 1.e-2))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }

    if (logToFile)
    {
        std::ofstream ofs("filteredHighPassSignal.csv");

        for (auto fs : filteredSignalHp)
        {
            ofs << fs << std::endl;
        }
    }

    std::cout << "[Test 3 - FirBandPassFilter]" << std::endl;

    t.Reset();
    filterCoeffs = FirBandPassFilter(451, 150., 100., 2000., kaiser);
    std::cout << "\t(FirBandPassFilter) Duration " << t.Elapsed() << "s" << std::endl;

    if (logToFile)
    {
        std::ofstream ofs("bandPass.csv");

        for (auto fc : filterCoeffs)
        {
            ofs << fc << std::endl;
        }
    }

    t.Reset();
    filter_hldr_d bpFilter(2000, filterCoeffs, false);
    std::cout << "\t(filter_hldr_d) Duration " << t.Elapsed() << "s" << std::endl;

    std::vector<double> filteredSignalBp(multiToneSignal1.size());

    t.Reset();
    bpFilter(multiToneSignal1.begin(), multiToneSignal1.end(), filteredSignalBp.begin(), true);
    std::cout << "\t(Apply filter) Duration " << t.Elapsed() << "s" << std::endl;

    stats     = ComputeStats(std::next(filteredSignalBp.begin(), 100),
                         std::next(filteredSignalBp.begin(), 1100));
    statsComp = ComputeStats(toneSignal2.begin(), toneSignal2.end());

    if ((std::abs(stats.min - statsComp.min) < 1.e-1) &&
        (std::abs(stats.max - statsComp.max) < 1.e-1) &&
        (std::abs(stats.mean - statsComp.mean) < 1.e-2) &&
        (std::abs(stats.stdDev - statsComp.stdDev) < 1.e-2))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }

    if (logToFile)
    {
        std::ofstream ofs("filtereBandPassSignal.csv");

        for (auto fs : filteredSignalHp)
        {
            ofs << fs << std::endl;
        }
    }

    std::cout << "[Test 4 - FirNotchFilter]" << std::endl;

    t.Reset();
    filterCoeffs = FirNotchFilter(451, 150., 10., 2000., kaiser);
    std::cout << "\t(FirNotchFilter) Duration " << t.Elapsed() << "s" << std::endl;

    if (logToFile)
    {
        std::ofstream ofs("notch.csv");

        for (auto fc : filterCoeffs)
        {
            ofs << fc << std::endl;
        }
    }

    t.Reset();
    filter_hldr_d nFilter(2000, filterCoeffs, false);
    std::cout << "\t(filter_hldr_d) Duration " << t.Elapsed() << "s" << std::endl;

    std::vector<double> filteredSignalN(multiToneSignal1.size());

    t.Reset();
    nFilter(multiToneSignal2.begin(), multiToneSignal2.end(), filteredSignalN.begin(), true);
    std::cout << "\t(Apply filter) Duration " << t.Elapsed() << "s" << std::endl;

    stats     = ComputeStats(std::next(filteredSignalN.begin(), 100),
                         std::next(filteredSignalN.begin(), 1100));
    statsComp = ComputeStats(toneSignal3.begin(), toneSignal3.end());

    if ((std::abs(stats.min - statsComp.min) < 1.e-1) &&
        (std::abs(stats.max - statsComp.max) < 1.e-1) &&
        (std::abs(stats.mean - statsComp.mean) < 1.e-2) &&
        (std::abs(stats.stdDev - statsComp.stdDev) < 1.e-2))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }

    if (logToFile)
    {
        std::ofstream ofs("filtereNotchSignal.csv");

        for (auto fs : filteredSignalN)
        {
            ofs << fs << std::endl;
        }
    }

    return result;
}

/*!
 * \brief Unit test function for GCD.
 * \return A pair reporting the number of tests failed and the number of tests run.
 */
std::pair<int, int> TestGcd()
{
    std::pair<int, int> result(0, 2);

    std::cout << "Testing dsp::Gcd function..." << std::endl;
    std::cout << "[Test 1 - Gcd(48, 36)]" << std::endl;

    size_t a = 48;
    size_t b = 36;

    test::Timer t;
    auto        g = Gcd(a, b);
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    if (12 == g)
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }

    std::cout << "[Test 2 - Gcd(2680, 496)]" << std::endl;

    a = 2680;
    b = 496;

    t.Reset();
    g = Gcd(a, b);
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    if (8 == g)
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }

    return result;
}

/*!
 * \brief Unit test function for resample function using fast FFT convolution.
 * \return A pair reporting the number of tests failed and the number of tests run.
 */
std::pair<int, int> TestResampling1(bool logToFile)
{
    std::pair<int, int> result(0, 6);

    std::cout << "Testing dsp::Resample class (fast convolution)..." << std::endl;
    std::cout << "[Test 1 - Resample by 93/13]" << std::endl;

    ToneParams t1{10., 1., 0., 0.};

    auto originalSignal = Tone<double>(t1, 100., 500);

    if (logToFile)
    {
        std::ofstream ofs("originalSignal1.csv");

        for (auto s : originalSignal)
        {
            ofs << s << std::endl;
        }
    }

    test::Timer      t;
    Resample<double> resampler1(originalSignal.size(), 93, 13, 100., 50., 1001, 10., true);
    std::cout << "\t(Resample<double>) Duration " << t.Elapsed() << "s" << std::endl;

    std::vector<double> resampledSignal1(resampler1.ResampledSize());

    t.Reset();
    resampler1(originalSignal.begin(), originalSignal.end(), resampledSignal1.begin());
    std::cout << "\tPerform resampling, duration " << t.Elapsed() << "s" << std::endl;

    if (logToFile)
    {
        std::ofstream ofs("resampledSignal1a.csv");

        for (auto fs : resampledSignal1)
        {
            ofs << fs << std::endl;
        }
    }

    auto stats     = ComputeStats(resampledSignal1.begin(), resampledSignal1.end());
    auto statsComp = ComputeStats(originalSignal.begin(), originalSignal.end());

    if ((std::abs(stats.min - statsComp.min) < 1.e-1) &&
        (std::abs(stats.max - statsComp.max) < 1.e-1) &&
        (std::abs(stats.mean - statsComp.mean) < 1.e-2) &&
        (std::abs(stats.stdDev - statsComp.stdDev) < 1.e-2))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }

    std::cout << "[Test 2 - Resample by 100/1]" << std::endl;

    t.Reset();
    Resample<double> resampler2(originalSignal.size(), 100, 1, 100., 50., 1001, 10., true);
    std::cout << "\t(Resample<double>) Duration " << t.Elapsed() << "s" << std::endl;

    std::vector<double> resampledSignal2(resampler2.ResampledSize());

    t.Reset();
    resampler2(originalSignal.begin(), originalSignal.end(), resampledSignal2.begin());
    std::cout << "\tPerform resampling, duration " << t.Elapsed() << "s" << std::endl;

    if (logToFile)
    {
        std::ofstream ofs("resampledSignal1b.csv");

        for (auto fs : resampledSignal2)
        {
            ofs << fs << std::endl;
        }
    }

    stats     = ComputeStats(resampledSignal2.begin(), resampledSignal2.end());
    statsComp = ComputeStats(originalSignal.begin(), originalSignal.end());

    if ((std::abs(stats.min - statsComp.min) < 1.e-1) &&
        (std::abs(stats.max - statsComp.max) < 1.e-1) &&
        (std::abs(stats.mean - statsComp.mean) < 1.e-2) &&
        (std::abs(stats.stdDev - statsComp.stdDev) < 1.e-2))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }

    std::cout << "[Test 3 - Resample by 1/5]" << std::endl;
    ToneParams t2{10., 100., 0., 0.};

    originalSignal = Tone<double>(t2, 10000., 5000);

    if (logToFile)
    {
        std::ofstream ofs("originalSignal2.csv");

        for (auto s : originalSignal)
        {
            ofs << s << std::endl;
        }
    }

    t.Reset();
    Resample<double> resampler3(originalSignal.size(), 1, 5, 10000., 1000., 1001, 10., true);
    std::cout << "\t(Resample<double>) Duration " << t.Elapsed() << "s" << std::endl;

    std::vector<double> resampledSignal3(resampler3.ResampledSize());

    t.Reset();
    resampler3(originalSignal.begin(), originalSignal.end(), resampledSignal3.begin());
    std::cout << "\tPerform resampling, duration " << t.Elapsed() << "s" << std::endl;

    if (logToFile)
    {
        std::ofstream ofs("resampledSignal2a.csv");

        for (auto fs : resampledSignal3)
        {
            ofs << fs << std::endl;
        }
    }

    stats     = ComputeStats(resampledSignal3.begin(), resampledSignal3.end());
    statsComp = ComputeStats(originalSignal.begin(), originalSignal.end());

    if ((std::abs(stats.min - statsComp.min) < 1.e-1) &&
        (std::abs(stats.max - statsComp.max) < 1.e-1) &&
        (std::abs(stats.mean - statsComp.mean) < 1.e-2) &&
        (std::abs(stats.stdDev - statsComp.stdDev) < 1.e-2))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }

    std::cout << "[Test 4 - Resample by 2/3]" << std::endl;

    t.Reset();
    Resample<double> resampler4(originalSignal.size(), 2, 3, 10000., 3333.333333, 1001, 10., true);
    std::cout << "\t(Resample<double>) Duration " << t.Elapsed() << "s" << std::endl;

    std::vector<double> resampledSignal4(resampler4.ResampledSize());

    t.Reset();
    resampler4(originalSignal.begin(), originalSignal.end(), resampledSignal4.begin());
    std::cout << "\tPerform resampling, duration " << t.Elapsed() << "s" << std::endl;

    if (logToFile)
    {
        std::ofstream ofs("resampledSignal2b.csv");

        for (auto fs : resampledSignal4)
        {
            ofs << fs << std::endl;
        }
    }

    stats     = ComputeStats(resampledSignal4.begin(), resampledSignal4.end());
    statsComp = ComputeStats(originalSignal.begin(), originalSignal.end());

    if ((std::abs(stats.min - statsComp.min) < 1.e-1) &&
        (std::abs(stats.max - statsComp.max) < 1.e-1) &&
        (std::abs(stats.mean - statsComp.mean) < 1.e-2) &&
        (std::abs(stats.stdDev - statsComp.stdDev) < 1.e-2))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }

    std::cout << "[Test 5 - Compute resample factors (27.65421)]" << std::endl;

    t.Reset();
    auto factors = dsp::ComputeResampleFactors(27.65421);
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    auto absDiff = std::abs(
        (static_cast<double>(factors.first) / static_cast<double>(factors.second)) - 27.65421);

    if (absDiff < 5.e-2)
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }

    std::cout << "[Test 6 - Compute resample factors (0.8659)]" << std::endl;

    t.Reset();
    factors = dsp::ComputeResampleFactors(0.8659);
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    absDiff = std::abs((static_cast<double>(factors.first) / static_cast<double>(factors.second)) -
                       0.8659);

    if (absDiff < 5.e-2)
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }

    return result;
}

/*!
 * \brief Unit test function for resample function using slow convolution.
 * \return A pair reporting the number of tests failed and the number of tests run.
 */
std::pair<int, int> TestResampling2(bool logToFile)
{
    std::pair<int, int> result(0, 6);

    std::cout << "Testing dsp::Resample class (slow convolution)..." << std::endl;
    std::cout << "[Test 1 - Resample by 93/13]" << std::endl;

    ToneParams t1{10., 1., 0., 0.};

    auto originalSignal = Tone<double>(t1, 100., 500);

    if (logToFile)
    {
        std::ofstream ofs("originalSignal1.csv");

        for (auto s : originalSignal)
        {
            ofs << s << std::endl;
        }
    }

    test::Timer      t;
    Resample<double> resampler1(originalSignal.size(), 93, 13, 100., 50., 1001, 10., false);
    std::cout << "\t(Resample<double>) Duration " << t.Elapsed() << "s" << std::endl;

    std::vector<double> resampledSignal1(resampler1.ResampledSize());

    t.Reset();
    resampler1(originalSignal.begin(), originalSignal.end(), resampledSignal1.begin());
    std::cout << "\tPerform resampling, duration " << t.Elapsed() << "s" << std::endl;

    if (logToFile)
    {
        std::ofstream ofs("resampledSignal1a.csv");

        for (auto fs : resampledSignal1)
        {
            ofs << fs << std::endl;
        }
    }

    auto stats     = ComputeStats(resampledSignal1.begin(), resampledSignal1.end());
    auto statsComp = ComputeStats(originalSignal.begin(), originalSignal.end());

    if ((std::abs(stats.min - statsComp.min) < 1.e-1) &&
        (std::abs(stats.max - statsComp.max) < 1.e-1) &&
        (std::abs(stats.mean - statsComp.mean) < 1.e-2) &&
        (std::abs(stats.stdDev - statsComp.stdDev) < 1.e-2))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }

    std::cout << "[Test 2 - Resample by 100/1]" << std::endl;

    t.Reset();
    Resample<double> resampler2(originalSignal.size(), 100, 1, 100., 50., 1001, 10., false);
    std::cout << "\t(Resample<double>) Duration " << t.Elapsed() << "s" << std::endl;

    std::vector<double> resampledSignal2(resampler2.ResampledSize());

    t.Reset();
    resampler2(originalSignal.begin(), originalSignal.end(), resampledSignal2.begin());
    std::cout << "\tPerform resampling, duration " << t.Elapsed() << "s" << std::endl;

    if (logToFile)
    {
        std::ofstream ofs("resampledSignal1b.csv");

        for (auto fs : resampledSignal2)
        {
            ofs << fs << std::endl;
        }
    }

    stats     = ComputeStats(resampledSignal2.begin(), resampledSignal2.end());
    statsComp = ComputeStats(originalSignal.begin(), originalSignal.end());

    if ((std::abs(stats.min - statsComp.min) < 1.e-1) &&
        (std::abs(stats.max - statsComp.max) < 1.e-1) &&
        (std::abs(stats.mean - statsComp.mean) < 1.e-2) &&
        (std::abs(stats.stdDev - statsComp.stdDev) < 1.e-2))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }

    std::cout << "[Test 3 - Resample by 1/5]" << std::endl;
    ToneParams t2{10., 100., 0., 0.};

    originalSignal = Tone<double>(t2, 10000., 5000);

    if (logToFile)
    {
        std::ofstream ofs("originalSignal2.csv");

        for (auto s : originalSignal)
        {
            ofs << s << std::endl;
        }
    }

    t.Reset();
    Resample<double> resampler3(originalSignal.size(), 1, 5, 10000., 1000., 1001, 10., false);
    std::cout << "\t(Resample<double>) Duration " << t.Elapsed() << "s" << std::endl;

    std::vector<double> resampledSignal3(resampler3.ResampledSize());

    t.Reset();
    resampler3(originalSignal.begin(), originalSignal.end(), resampledSignal3.begin());
    std::cout << "\tPerform resampling, duration " << t.Elapsed() << "s" << std::endl;

    if (logToFile)
    {
        std::ofstream ofs("resampledSignal2a.csv");

        for (auto fs : resampledSignal3)
        {
            ofs << fs << std::endl;
        }
    }

    stats     = ComputeStats(resampledSignal3.begin(), resampledSignal3.end());
    statsComp = ComputeStats(originalSignal.begin(), originalSignal.end());

    if ((std::abs(stats.min - statsComp.min) < 1.e-1) &&
        (std::abs(stats.max - statsComp.max) < 1.e-1) &&
        (std::abs(stats.mean - statsComp.mean) < 1.e-2) &&
        (std::abs(stats.stdDev - statsComp.stdDev) < 1.e-2))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }

    std::cout << "[Test 4 - Resample by 2/3]" << std::endl;

    t.Reset();
    Resample<double> resampler4(
        originalSignal.size(), 2, 3, 10000., 3333.3333333, 1001, 10., false);
    std::cout << "\t(Resample<double>) Duration " << t.Elapsed() << "s" << std::endl;

    std::vector<double> resampledSignal4(resampler4.ResampledSize());

    t.Reset();
    resampler4(originalSignal.begin(), originalSignal.end(), resampledSignal4.begin());
    std::cout << "\tPerform resampling, duration " << t.Elapsed() << "s" << std::endl;

    if (logToFile)
    {
        std::ofstream ofs("resampledSignal2b.csv");

        for (auto fs : resampledSignal4)
        {
            ofs << fs << std::endl;
        }
    }

    stats     = ComputeStats(resampledSignal4.begin(), resampledSignal4.end());
    statsComp = ComputeStats(originalSignal.begin(), originalSignal.end());

    if ((std::abs(stats.min - statsComp.min) < 1.e-1) &&
        (std::abs(stats.max - statsComp.max) < 1.e-1) &&
        (std::abs(stats.mean - statsComp.mean) < 1.e-2) &&
        (std::abs(stats.stdDev - statsComp.stdDev) < 1.e-2))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }

    std::cout << "[Test 5 - Compute resample factors (27.65421)]" << std::endl;

    t.Reset();
    auto factors = dsp::ComputeResampleFactors(27.65421);
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    auto absDiff = std::abs(
        (static_cast<double>(factors.first) / static_cast<double>(factors.second)) - 27.65421);

    if (absDiff < 5.e-2)
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }

    std::cout << "[Test 6 - Compute resample factors (0.8659)]" << std::endl;

    t.Reset();
    factors = dsp::ComputeResampleFactors(0.8659);
    std::cout << "\tDuration " << t.Elapsed() << "s" << std::endl;

    absDiff = std::abs((static_cast<double>(factors.first) / static_cast<double>(factors.second)) -
                       0.8659);

    if (absDiff < 5.e-2)
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }

    return result;
}

/*!
 * \brief Unit test function for FFT magnitude functions.
 * \return A pair reporting the number of tests failed and the number of tests run.
 */
std::pair<int, int> TestComplexFftToMagnitude(bool logToFile)
{
    std::pair<int, int> result(0, 2);
    double              totalTime = 0.0;

    std::cout << "Testing dsp::ComplexFft class (to magnitude FFT)..." << std::endl;

    // Create signal.
    ToneParams t1{10., 3000., 0., 0.};
    ToneParams t2{5., 6000., 0., 0.};
    ToneParams t3{2., 12000., 0., 0.};

    auto multiToneSignal = MultiTone<double>({t1, t2, t3}, 256000., 1024);

    // Create window function.
    std::cout << "[Compute Window Function]" << std::endl;
    test::Timer t;
    window_fn_d window(HannGenerator{}, multiToneSignal.size() + 1, true);
    auto        elapsed = t.Elapsed();
    std::cout << "\tDuration " << elapsed << "s" << std::endl;
    totalTime += elapsed;

    // Apply window function.
    std::cout << "[Apply Window Function to data]" << std::endl;
    t.Reset();
    window(multiToneSignal.begin(), multiToneSignal.end(), multiToneSignal.begin());
    elapsed = t.Elapsed();
    std::cout << "\tDuration " << elapsed << "s" << std::endl;
    totalTime += elapsed;

    // Generate forward complex FFT result.
    std::cout << "[Compute forward complex FFT]" << std::endl;
    t.Reset();
    auto complexData = complex_fft_d::Forward(multiToneSignal.begin(), multiToneSignal.end());
    elapsed          = t.Elapsed();
    std::cout << "\tDuration " << elapsed << "s" << std::endl;
    totalTime += elapsed;

    // Normalise complex FFT data.
    std::cout << "[Normalise complex FFT]" << std::endl;
    t.Reset();
    complex_fft_d::Normalise(complexData);
    elapsed = t.Elapsed();
    std::cout << "\tDuration " << elapsed << "s" << std::endl;
    totalTime += elapsed;

    // Zero bottom 3 BINs.
    complexData[0] = {0, 0};
    complexData[1] = {0, 0};
    complexData[2] = {0, 0};

    // Correct window coherent gain.
    std::cout << "[Gain correct complex FFT]" << std::endl;
    t.Reset();
    window_fn_d::ApplyGainCorrection(
        complexData.begin(), complexData.end(), complexData.begin(), window.CoherentGain());
    elapsed = t.Elapsed();
    std::cout << "\tDuration " << elapsed << "s" << std::endl;
    totalTime += elapsed;

    // To magnitude real valued result.
    std::cout << "[Compute magnitude (real) FFT]" << std::endl;
    t.Reset();
    complex_fft_d::ToMagnitude(complexData);
    elapsed = t.Elapsed();
    std::cout << "\tDuration " << elapsed << "s" << std::endl;
    totalTime += elapsed;

    if (logToFile)
    {
        std::ofstream ofs("MagnitudeFft.csv");
        auto          halfSize = complexData.size() / 2;

        for (size_t i = 0; i < halfSize; ++i)
        {
            ofs << complexData[i].real() << std::endl;
        }
    }

    auto bin13DIff = std::abs(complexData[12].real() - 10.);
    auto bin25Diff = std::abs(complexData[24].real() - 5.);
    auto bin49Diff = std::abs(complexData[48].real() - 2.);

    std::cout << "\tTotal FFT Processing Duration " << totalTime << "s" << std::endl;

    if ((bin13DIff < 0.1) && (bin25Diff < 0.1) && (bin49Diff < 0.1))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }

    std::cout << "[Test 2 - Magnitude Fft functor]" << std::endl;
    multiToneSignal = MultiTone<double>({t1, t2, t3}, 256000., 1024);

    std::cout << "[Create functor]" << std::endl;
    t.Reset();
    magnitude_fft_d magnitudeFft(HannGenerator{}, multiToneSignal.size());
    std::cout << "\tDuration " << elapsed << "s" << std::endl;
    totalTime += elapsed;

    std::vector<double> resultSpectrum(multiToneSignal.size());

    std::cout << "[Execute functor]" << std::endl;
    t.Reset();
    magnitudeFft(multiToneSignal.begin(), multiToneSignal.end(), resultSpectrum);
    std::cout << "\tDuration " << elapsed << "s" << std::endl;
    totalTime += elapsed;

    bin13DIff = std::abs(resultSpectrum[12] - 10.);
    bin25Diff = std::abs(resultSpectrum[24] - 5.);
    bin49Diff = std::abs(resultSpectrum[48] - 2.);

    if ((bin13DIff < 0.1) && (bin25Diff < 0.1) && (bin49Diff < 0.1))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }

    return result;
}

/*!
 * \brief Unit test function for FFT 3 BIN Sum functions.
 * \return A pair reporting the number of tests failed and the number of tests run.
 */
std::pair<int, int> TestComplexFftTo3BinSum(bool logToFile)
{
    std::pair<int, int> result(0, 2);
    double totalTime = 0.0;

    std::cout << "Testing dsp::ComplexFft class (to 3-BIN sum FFT)..." << std::endl;

    // Create signal.
    ToneParams t1{10., 3000., 0., 0.};
    ToneParams t2{5., 6000., 0., 0.};
    ToneParams t3{2., 12000., 0., 0.};

    auto multiToneSignal = MultiTone<double>({t1, t2, t3}, 256000., 1024);

    std::cout << "[Test 1 - 3Bin sum FFT]" << std::endl;
    std::cout << "[Compute Window Function]" << std::endl;
    test::Timer t;
    window_fn_d window(HannGenerator{}, multiToneSignal.size() + 1, true);
    auto        elapsed = t.Elapsed();
    std::cout << "\tDuration " << elapsed << "s" << std::endl;
    totalTime += elapsed;

    // Apply window function.
    std::cout << "[Apply Window Function to data]" << std::endl;
    t.Reset();
    window(multiToneSignal.begin(), multiToneSignal.end(), multiToneSignal.begin());
    elapsed = t.Elapsed();
    std::cout << "\tDuration " << elapsed << "s" << std::endl;
    totalTime += elapsed;

    // Generate forward complex FFT result.
    std::cout << "[Compute forward complex FFT]" << std::endl;
    t.Reset();
    auto complexData = complex_fft_d::Forward(multiToneSignal.begin(), multiToneSignal.end());
    elapsed          = t.Elapsed();
    std::cout << "\tDuration " << elapsed << "s" << std::endl;
    totalTime += elapsed;

    // Normalise complex FFT data.
    std::cout << "[Normalise complex FFT]" << std::endl;
    t.Reset();
    complex_fft_d::Normalise(complexData);
    elapsed = t.Elapsed();
    std::cout << "\tDuration " << elapsed << "s" << std::endl;
    totalTime += elapsed;

    // Zero bottom 3 BINs.
    complexData[0] = {0, 0};
    complexData[1] = {0, 0};
    complexData[2] = {0, 0};

    // Convert to power spectrum.
    std::cout << "[Compute power spectrum of complex FFT]" << std::endl;
    t.Reset();
    complex_fft_d::ToPower(complexData);
    elapsed = t.Elapsed();
    std::cout << "\tDuration " << elapsed << "s" << std::endl;
    totalTime += elapsed;

    if (logToFile)
    {
        std::ofstream ofs("PowerSpectrum.csv");
        auto          halfSize = complexData.size() / 2;

        for (size_t i = 0; i < halfSize; ++i)
        {
            ofs << complexData[i].real() << std::endl;
        }
    }

    // Correct gain of power spectrum
    std::cout << "[Gain correct power spectrum (real) FFT]" << std::endl;
    t.Reset();
    window_fn_d::ApplyGainCorrection(
        complexData.begin(), complexData.end(), complexData.begin(), window.CombinedGain());
    elapsed = t.Elapsed();
    std::cout << "\tDuration " << elapsed << "s" << std::endl;
    totalTime += elapsed;

    if (logToFile)
    {
        std::ofstream ofs("CorrectedPowerSpectrum.csv");
        auto          halfSize = complexData.size() / 2;

        for (size_t i = 0; i < halfSize; ++i)
        {
            ofs << complexData[i].real() << std::endl;
        }
    }

    // Convert to magnitude using 3-BIN summation.
    std::cout << "[Compute 3-BIN summed (real) FFT]" << std::endl;
    t.Reset();
    complex_fft_d::To3BinSum(complexData);
    elapsed = t.Elapsed();
    std::cout << "\tDuration " << elapsed << "s" << std::endl;
    totalTime += elapsed;

    auto halfSize = complexData.size() / 2;

    if (logToFile)
    {
        std::ofstream ofs("3BinSumFft.csv");

        for (size_t i = 0; i < halfSize; ++i)
        {
            ofs << complexData[i].real() << std::endl;
        }
    }

    std::cout << "\tTotal FFT Processing Duration " << totalTime << "s" << std::endl;

    auto bin13DIff = std::abs(complexData[12].real() - 10.);
    auto bin25Diff = std::abs(complexData[24].real() - 5.);
    auto bin49Diff = std::abs(complexData[48].real() - 2.);

    if ((bin13DIff < 0.1) && (bin25Diff < 0.1) && (bin49Diff < 0.1))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }

    std::cout << "[Test 2 - 3Bin sum FFT functor]" << std::endl;
    multiToneSignal = MultiTone<double>({t1, t2, t3}, 256000., 1024);

    std::cout << "[Create functor]" << std::endl;
    t.Reset();
    three_bin_sum_fft_d threeBinSumFft(HannGenerator{}, multiToneSignal.size());
    std::cout << "\tDuration " << elapsed << "s" << std::endl;
    totalTime += elapsed;

    std::vector<double> resultSpectrum(multiToneSignal.size());

    std::cout << "[Execute functor]" << std::endl;
    t.Reset();
    threeBinSumFft(multiToneSignal.begin(), multiToneSignal.end(), resultSpectrum);
    std::cout << "\tDuration " << elapsed << "s" << std::endl;
    totalTime += elapsed;

    bin13DIff = std::abs(resultSpectrum[12] - 10.);
    bin25Diff = std::abs(resultSpectrum[24] - 5.);
    bin49Diff = std::abs(resultSpectrum[48] - 2.);

    if ((bin13DIff < 0.1) && (bin25Diff < 0.1) && (bin49Diff < 0.1))
    {
        std::cout << "Test passed? = true" << std::endl;
    }
    else
    {
        ++result.first;
        std::cout << "Test passed? = false" << std::endl;
    }

    return result;
}

} // namespace dsp

/*! \brief Test harness main function. */
int main()
{
    bool                logToFile = true;
    std::pair<int, int> totals(0, 0);

    auto result = dsp::TestConvolve();
    totals.first += result.first;
    totals.second += result.second;

    std::cout << "Tests failed = " << result.first << ", Tests run = " << result.second << std::endl
              << std::endl;

    result = dsp::TestBessel();
    totals.first += result.first;
    totals.second += result.second;

    std::cout << "Tests failed = " << result.first << ", Tests run = " << result.second << std::endl
              << std::endl;

    result = dsp::TestSinc();
    totals.first += result.first;
    totals.second += result.second;

    std::cout << "Tests failed = " << result.first << ", Tests run = " << result.second << std::endl
              << std::endl;

    result = dsp::TestSine();
    totals.first += result.first;
    totals.second += result.second;

    std::cout << "Tests failed = " << result.first << ", Tests run = " << result.second << std::endl
              << std::endl;

    result = dsp::TestWindowFunction();
    totals.first += result.first;
    totals.second += result.second;

    std::cout << "Tests failed = " << result.first << ", Tests run = " << result.second << std::endl
              << std::endl;

    result = dsp::TestFilters1(logToFile);
    totals.first += result.first;
    totals.second += result.second;

    result = dsp::TestFilters2(logToFile);
    totals.first += result.first;
    totals.second += result.second;

    std::cout << "Tests failed = " << result.first << ", Tests run = " << result.second << std::endl
              << std::endl;

    result = dsp::TestGcd();
    totals.first += result.first;
    totals.second += result.second;

    std::cout << "Tests failed = " << result.first << ", Tests run = " << result.second << std::endl
              << std::endl;

    result = dsp::TestComplexFftToMagnitude(logToFile);
    totals.first += result.first;
    totals.second += result.second;

    std::cout << "Tests failed = " << result.first << ", Tests run = " << result.second << std::endl
              << std::endl;

    result = dsp::TestComplexFftTo3BinSum(logToFile);
    totals.first += result.first;
    totals.second += result.second;

    std::cout << "Tests failed = " << result.first << ", Tests run = " << result.second << std::endl
              << std::endl;

    result = dsp::TestResampling1(logToFile);
    totals.first += result.first;
    totals.second += result.second;

    std::cout << "Tests failed = " << result.first << ", Tests run = " << result.second << std::endl
              << std::endl;

    result = dsp::TestResampling2(logToFile);
    totals.first += result.first;
    totals.second += result.second;

    std::cout << "Tests failed = " << result.first << ", Tests run = " << result.second << std::endl
              << std::endl;

    std::cout << "Total tests failed = " << totals.first << ", Total tests run = " << totals.second
              << std::endl
              << std::endl;

    std::cout << "Press any key followed by enter to quit." << std::endl;

    char temp;
    std::cin >> temp;

    return 0;
}