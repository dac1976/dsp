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

/*!
 * \file dsp_math.hpp
 * \brief File containing generic math algorithms and functions.
 */

#ifndef DSP_MATH_HPP
#define DSP_MATH_HPP

#include <iterator>
#include <vector>
#include "dsp_errors.hpp"
#include "dsp_pi.hpp"

/*! \brief dsp namespace */
namespace dsp
{

/*!
 * \brief Perform discrete convolution on 2 data ranges.
 * \param[in] first1 - Start of the first range.
 * \param[in] last1 - End of the first range.
 * \param[in] first2 - Start of the second range.
 * \param[in] last2 - End of the second range.
 * \param[out] result - Start of the result range.
 *
 * Iterator result must point to a range with space for
 * M + N - 1 elements, where range 1 contains M elements
 * and range 2 contains N elements. None of the input or
 * output ranges should overlap otherwise the behaviour
 * is undefined.
 *
 * Formula for the discrete convolution, y(n), of finite
 * sequences x(n) and h(n), where x(n) has M values and
 * h(n) has N values, is:
 *
 * y(n) = SUM[k=0..n]{x(k)*h(N-1-k)}
 * for n=0..M+N-2
 */
template <typename InIter1, typename InIter2, typename OutIter>
void Convolve(InIter1 first1, InIter1 last1, InIter2 first2, InIter2 last2, OutIter result)
{
    using value_t = typename std::iterator_traits<OutIter>::value_type;

    auto M = std::distance(first1, last1);
    DSP_ASSERT_THROW(M > 0, "range 1 invalid");

    auto N = std::distance(first2, last2);
    DSP_ASSERT_THROW(N > 0, "range 2 invalid");

    auto K = M + N - 1;

    for (auto k = 0; k < K; ++k)
    {
        *result   = value_t(0);
        auto kMin = k >= N - 1 ? k - (N - 1) : 0;
        auto kMax = k < M - 1 ? k : M - 1;

        for (auto kSub = kMin; kSub <= kMax; ++kSub)
        {
            *result +=
                static_cast<value_t>(*std::next(first1, kSub) * *std::next(first2, k - kSub));
        }

        std::advance(result, 1);
    }
}

/*!
 * \brief Zeroth-order modified Bessel function of the first kind
 * \param[in] x - Value to compute bessel function for.
 * \return Bessel value of x.
 *
 * We'll use this to create a Kaiser window.
 */
template <typename FloatType> FloatType Bessel(FloatType x)
{
    static_assert(std::is_floating_point<FloatType>::value, "invalid floating point type");

    FloatType sum{0};

    for (int i = 1; i < 10; ++i)
    {
        auto xToIpower = std::pow(x / 2, static_cast<FloatType>(i));
        int  factorial = 1;

        for (int j = 1; j <= i; ++j)
        {
            factorial *= j;
        }

        sum += std::pow(xToIpower / factorial, static_cast<FloatType>(2));
    }

    return sum + static_cast<FloatType>(1);
}

/*!
 * \brief Sinc function.
 * \param[in] x - Value to compute sinc function for.
 * \param[in] limitThreshold - at x = 0 sinc(x) is defined as 1 use this limit to control when to
 * set sinc = 1.
 * \return Sinc value.
 *
 * This is the unnormalised "classic" version of sinc.
 */
template <typename FloatType>
FloatType Sinc(FloatType x, FloatType limitThreshold = FloatType(1.e-9))
{
    static_assert(std::is_floating_point<FloatType>::value, "invalid floating point type");

    if (std::abs(x) < limitThreshold)
    {
        return static_cast<FloatType>(1);
    }

    return sin(x) / x;
}

/*!
 * \brief Normalised sinc function.
 * \param[in] x - Value to compute sinc function for.
 * \param[in] limitThreshold - at x = 0 sinc(x) is defined as 1 use this limit to control when to
 * set sinc = 1.
 * \return Normalised sinc value.
 *
 * This is the normalised version of sinc.
 */
template <typename FloatType>
FloatType SincNorm(FloatType x, FloatType limitThreshold = FloatType(1.e-9))
{
    static_assert(std::is_floating_point<FloatType>::value, "invalid floating point type");

    if (std::abs(x) < limitThreshold)
    {
        return FloatType(1);
    }

    auto X = Pi<FloatType>() * x;
    return sin(X) / X;
}

/*!
 * \brief Sinusoidal equation.
 * \param[in] amplitude - Peak amplitude of sine wave.
 * \param[in] time - Time point in seconds.
 * \param[in] frequency - Signal frequency in Hz.
 * \param[in] phase - Phase offset in radians.
 * \param[in] offset - Amplitude offset.
 * \return Point on sine wave for given time.
 *
 * y(t) = A.sin(2.Pi.f.t + p) + o
 */
template <typename FloatType>
FloatType Sine(FloatType amplitude, FloatType time, FloatType frequency, FloatType phase,
               FloatType offset)
{
    static_assert(std::is_floating_point<FloatType>::value, "invalid floating point type");
    return (amplitude * sin((TwoPi<FloatType>() * frequency * time) + phase)) + offset;
}

/*!
 * \brief Computes GCD of two unsigned values using Binary GCD algorithm.
 * \param[in] a - First value.
 * \param[in] b - Second value.
 * \return The GCD of the two inputs.
 */
template <typename UintType> UintType Gcd(UintType a, UintType b)
{
    static_assert(std::is_unsigned<UintType>::value, "invalid unsigned type");

    if (a == b)
    {
        return a;
    }

    if (0 == a)
    {
        return b;
    }

    if (0 == b)
    {
        return a;
    }

    if (~a & 1)
    {
        if (b & 1)
        {
            return Gcd(a >> 1, b);
        }
        else
        {
            return Gcd(a >> 1, b >> 1) << 1;
        }
    }

    if (~b & 1)
    {
        return Gcd(a, b >> 1);
    }

    if (a > b)
    {
        return Gcd((a - b) >> 1, b);
    }

    return Gcd((b - a) >> 1, a);
}

/*!
 * \brief Check if a positive integer is a power of 2.
 * \param[in] n - the number to test.
 * \return True if a power of 2 false otherwise.
 */
template <typename IntType> bool IsPowerOf2(IntType n)
{
    static_assert(std::is_integral<IntType>::value, "invalid integer type");

    return (n > 0) && ((n & (n - 1)) == 0);
}

} // namespace dsp

#endif // DSP_MATH_HPP