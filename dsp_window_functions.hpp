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
 * \file dsp_window_functions.hpp
 * \brief File containing generic implementations of window functions.
 */
#ifndef DSP_WINDOW_FUNCTIONS_HPP
#define DSP_WINDOW_FUNCTIONS_HPP

#include <array>
#include <functional>
#include <algorithm>
#include "dsp_math.hpp"

/*! \brief dsp namespace */
namespace dsp
{

/*!
 * \brief Generate a vector of window coefficients.
 * \param[in,out] windowCoeffs - Pre-sized vector to receive window coefficients.
 * \param[in] evalCoeff - function or functor to compute a single window coefficient.
 *
 * Uses symmetry to optimise computation time needed to generate the
 * window coefficients.
 */
template <typename FloatType, typename EvalFunc>
void WindowGenerator(std::vector<FloatType>& windowCoeffs, EvalFunc evalCoeff)
{
    auto size = windowCoeffs.size();
    DSP_ASSERT_THROW(size > 1, "invalid window size");
    auto sizeMinusOne = static_cast<FloatType>(size - 1);
    auto halfSize     = size >> 1;

    for (size_t n = 0, nRev = size - 1; n < halfSize; ++n, --nRev)
    {
        windowCoeffs[n]    = evalCoeff(static_cast<FloatType>(n), sizeMinusOne);
        windowCoeffs[nRev] = windowCoeffs[n];
    }

    if (1 == size % 2)
    {
        windowCoeffs[halfSize] = evalCoeff(static_cast<FloatType>(halfSize), sizeMinusOne);
    }
}

/*!
 * \brief Flat top equation coefficient evaluator
 * \param[in] n - The current window coefficient index as a floating point type T.
 * \param[in] sizeMinusOne - The window length minus 1 as a floating point type T.
 * \param[in] eqCoeffsFirst - Iterator to first flat-top equation coefficient.
 * \param[in] eqCoeffsLast - Iterator to one past the final flat-top equation coefficient.
 * \return The nth flat-top window coefficient as a floating point type T.
 *
 * Formula:
 *
 * w(n) = a0 - a1.cos(2.pi.n / N-1)
 *           + a2.cos(4.pi.n / N-1)
 *           - a3.cos(6.pi.n / N-1)
 *           + a4.cos(8.pi.n / N-1)
 *           - ...
 */
template <typename FloatType, typename Iter>
FloatType EvaluateFlatTopCoefficient(FloatType n, FloatType sizeMinusOne, Iter eqCoeffsFirst,
                                     Iter eqCoeffsLast)
{
    static const auto TWO_PI      = TwoPi<FloatType>();
    auto              numEqCoeffs = std::distance(eqCoeffsFirst, eqCoeffsLast);
    DSP_ASSERT_THROW(numEqCoeffs > 1, "invalid number of equation coefficients");
    auto twoPiN      = TWO_PI * n;
    auto sign        = static_cast<FloatType>(-1);
    auto windowCoeff = *eqCoeffsFirst;
    int  i           = 1;

    for (auto itr = std::next(eqCoeffsFirst); itr != eqCoeffsLast;
         std::advance(itr, 1), sign *= static_cast<FloatType>(-1), ++i)
    {
        windowCoeff += sign * *itr * std::cos((static_cast<FloatType>(i) * twoPiN) / sizeMinusOne);
    }

    return windowCoeff;
}

/*! \brief Flat-top generator: ISO 18431-1. */
class FlatTop1Generator final
{
public:
    /*!
     * \brief Generate a vector of window coefficients.
     * \param[in,out] windowCoeffs - Pre-sized vector to receive window coefficients.
     */
    template <typename FloatType> void operator()(std::vector<FloatType>& windowCoeffs) const
    {
        WindowGenerator(windowCoeffs,
                        std::bind(&FlatTop1Generator::EvaluateCoefficient<FloatType>,
                                  std::placeholders::_1,
                                  std::placeholders::_2));
    }

private:
    /*!
     * \brief Evaluate nth window coefficient.
     * \param[in] n - The current window coefficient index as a floating point type T.
     * \param[in] sizeMinusOne - The window length minus 1 as a floating point type T.
     * \return The nth window coefficient as a floating type T.
     */
    template <typename FloatType>
    static FloatType EvaluateCoefficient(FloatType n, FloatType sizeMinusOne)
    {
        static const std::array<FloatType, 5> eqCoeffs = {static_cast<FloatType>(1),
                                                          static_cast<FloatType>(1.933),
                                                          static_cast<FloatType>(1.286),
                                                          static_cast<FloatType>(0.388),
                                                          static_cast<FloatType>(0.0322)};
        return EvaluateFlatTopCoefficient(n, sizeMinusOne, eqCoeffs.begin(), eqCoeffs.end());
    }
};

/*! \brief Flat-top generator: 2 point. */
class FlatTop2Generator final
{
public:
    /*!
     * \brief Generate a vector of window coefficients.
     * \param[in,out] windowCoeffs - Pre-sized vector to receive window coefficients.
     */
    template <typename FloatType> void operator()(std::vector<FloatType>& windowCoeffs) const
    {
        WindowGenerator(windowCoeffs,
                        std::bind(&FlatTop2Generator::EvaluateCoefficient<FloatType>,
                                  std::placeholders::_1,
                                  std::placeholders::_2));
    }

private:
    /*!
     * \brief Evaluate nth window coefficient.
     * \param[in] n - The current window coefficient index as a floating point type T.
     * \param[in] sizeMinusOne - The window length minus 1 as a floating point type T.
     * \return The nth window coefficient as a floating type T.
     */
    template <typename FloatType>
    static FloatType EvaluateCoefficient(FloatType n, FloatType sizeMinusOne)
    {
        static const std::array<FloatType, 3> eqCoeffs = {static_cast<FloatType>(0.2810639),
                                                          static_cast<FloatType>(0.5208972),
                                                          static_cast<FloatType>(0.1980399)};
        return EvaluateFlatTopCoefficient(n, sizeMinusOne, eqCoeffs.begin(), eqCoeffs.end());
    }
};

/*! \brief Flat-top generator: alternate 4 point. */
class FlatTop3Generator final
{
public:
    /*!
     * \brief Generate a vector of window coefficients.
     * \param[in,out] windowCoeffs - Pre-sized vector to receive window coefficients.
     */
    template <typename FloatType> void operator()(std::vector<FloatType>& windowCoeffs) const
    {
        WindowGenerator(windowCoeffs,
                        std::bind(&FlatTop3Generator::EvaluateCoefficient<FloatType>,
                                  std::placeholders::_1,
                                  std::placeholders::_2));
    }

private:
    /*!
     * \brief Evaluate nth window coefficient.
     * \param[in] n - The current window coefficient index as a floating point type T.
     * \param[in] sizeMinusOne - The window length minus 1 as a floating point type T.
     * \return The nth window coefficient as a floating type T.
     */
    template <typename FloatType>
    static FloatType EvaluateCoefficient(FloatType n, FloatType sizeMinusOne)
    {
        static const std::array<FloatType, 5> eqCoeffs = {static_cast<FloatType>(0.21557895),
                                                          static_cast<FloatType>(0.41663158),
                                                          static_cast<FloatType>(0.277263158),
                                                          static_cast<FloatType>(0.083578947),
                                                          static_cast<FloatType>(0.006947368)};
        return EvaluateFlatTopCoefficient(n, sizeMinusOne, eqCoeffs.begin(), eqCoeffs.end());
    }
};

/*! \brief Flat-top generator: 3 point HP P301. */
class FlatTop4Generator final
{
public:
    /*!
     * \brief Generate a vector of window coefficients.
     * \param[in,out] windowCoeffs - Pre-sized vector to receive window coefficients.
     */
    template <typename FloatType> void operator()(std::vector<FloatType>& windowCoeffs) const
    {
        WindowGenerator(windowCoeffs,
                        std::bind(&FlatTop4Generator::EvaluateCoefficient<FloatType>,
                                  std::placeholders::_1,
                                  std::placeholders::_2));
    }

private:
    /*!
     * \brief Evaluate nth window coefficient.
     * \param[in] n - The current window coefficient index as a floating point type T.
     * \param[in] sizeMinusOne - The window length minus 1 as a floating point type T.
     * \return The nth window coefficient as a floating type T.
     */
    template <typename FloatType>
    static FloatType EvaluateCoefficient(FloatType n, FloatType sizeMinusOne)
    {
        static const std::array<FloatType, 4> eqCoeffs = {static_cast<FloatType>(0.9994484),
                                                          static_cast<FloatType>(1.911456),
                                                          static_cast<FloatType>(1.076578),
                                                          static_cast<FloatType>(0.183162)};
        return EvaluateFlatTopCoefficient(n, sizeMinusOne, eqCoeffs.begin(), eqCoeffs.end());
    }
};

/*! \brief Flat-top generator: HP 4 point. */
class FlatTop5Generator final
{
public:
    /*!
     * \brief Generate a vector of window coefficients.
     * \param[in,out] windowCoeffs - Pre-sized vector to receive window coefficients.
     */
    template <typename FloatType> void operator()(std::vector<FloatType>& windowCoeffs) const
    {
        WindowGenerator(windowCoeffs,
                        std::bind(&FlatTop5Generator::EvaluateCoefficient<FloatType>,
                                  std::placeholders::_1,
                                  std::placeholders::_2));
    }

private:
    /*!
     * \brief Evaluate nth window coefficient.
     * \param[in] n - The current window coefficient index as a floating point type T.
     * \param[in] sizeMinusOne - The window length minus 1 as a floating point type T.
     * \return The nth window coefficient as a floating type T.
     */
    template <typename FloatType>
    static FloatType EvaluateCoefficient(FloatType n, FloatType sizeMinusOne)
    {
        static const std::array<FloatType, 5> eqCoeffs = {static_cast<FloatType>(1),
                                                          static_cast<FloatType>(1.869032),
                                                          static_cast<FloatType>(1.195972),
                                                          static_cast<FloatType>(0.035928),
                                                          static_cast<FloatType>(0.030916)};
        return EvaluateFlatTopCoefficient(n, sizeMinusOne, eqCoeffs.begin(), eqCoeffs.end());
    }
};

/*! \brief Flat-top generator: Modified HP P401 5 point. */
class FlatTop6Generator final
{
public:
    /*!
     * \brief Generate a vector of window coefficients.
     * \param[in,out] windowCoeffs - Pre-sized vector to receive window coefficients.
     */
    template <typename FloatType> void operator()(std::vector<FloatType>& windowCoeffs) const
    {
        WindowGenerator(windowCoeffs,
                        std::bind(&FlatTop6Generator::EvaluateCoefficient<FloatType>,
                                  std::placeholders::_1,
                                  std::placeholders::_2));
    }

private:
    /*!
     * \brief Evaluate nth window coefficient.
     * \param[in] n - The current window coefficient index as a floating point type T.
     * \param[in] sizeMinusOne - The window length minus 1 as a floating point type T.
     * \return The nth window coefficient as a floating type T.
     */
    template <typename FloatType>
    static FloatType EvaluateCoefficient(FloatType n, FloatType sizeMinusOne)
    {
        static const std::array<FloatType, 6> eqCoeffs = {static_cast<FloatType>(1),
                                                          static_cast<FloatType>(1.93774046310203),
                                                          static_cast<FloatType>(1.32530734987255),
                                                          static_cast<FloatType>(0.43206975880342),
                                                          static_cast<FloatType>(0.04359135851569),
                                                          static_cast<FloatType>(0.00015175580171)};
        return EvaluateFlatTopCoefficient(n, sizeMinusOne, eqCoeffs.begin(), eqCoeffs.end());
    }
};

/*! \brief Flat-top generator: Rohde & Schwartz 4 point. */
class FlatTop7Generator final
{
public:
    /*!
     * \brief Generate a vector of window coefficients.
     * \param[in,out] windowCoeffs - Pre-sized vector to receive window coefficients.
     */
    template <typename FloatType> void operator()(std::vector<FloatType>& windowCoeffs) const
    {
        WindowGenerator(windowCoeffs,
                        std::bind(&FlatTop7Generator::EvaluateCoefficient<FloatType>,
                                  std::placeholders::_1,
                                  std::placeholders::_2));
    }

private:
    /*!
     * \brief Evaluate nth window coefficient.
     * \param[in] n - The current window coefficient index as a floating point type T.
     * \param[in] sizeMinusOne - The window length minus 1 as a floating point type T.
     * \return The nth window coefficient as a floating type T.
     */
    template <typename FloatType>
    static FloatType EvaluateCoefficient(FloatType n, FloatType sizeMinusOne)
    {
        static const std::array<FloatType, 5> eqCoeffs = {static_cast<FloatType>(0.1881999),
                                                          static_cast<FloatType>(0.36923),
                                                          static_cast<FloatType>(0.28702),
                                                          static_cast<FloatType>(0.13077),
                                                          static_cast<FloatType>(0.02488)};
        return EvaluateFlatTopCoefficient(n, sizeMinusOne, eqCoeffs.begin(), eqCoeffs.end());
    }
};

/*! \brief Hann window generator. */
class HannGenerator final
{
public:
    /*!
     * \brief Generate a vector of window coefficients.
     * \param[in,out] windowCoeffs - Pre-sized vector to receive window coefficients.
     */
    template <typename FloatType> void operator()(std::vector<FloatType>& windowCoeffs) const
    {
        WindowGenerator(windowCoeffs,
                        std::bind(&HannGenerator::EvaluateCoefficient<FloatType>,
                                  std::placeholders::_1,
                                  std::placeholders::_2));
    }

private:
    /*!
     * \brief Evaluate nth window coefficient.
     * \param[in] n - The current window coefficient index as a floating point type T.
     * \param[in] sizeMinusOne - The window length minus 1 as a floating point type T.
     * \return The nth window coefficient as a floating type T.
     */
    template <typename FloatType>
    static FloatType EvaluateCoefficient(FloatType n, FloatType sizeMinusOne)
    {
        static const auto TWO_PI = TwoPi<FloatType>();
        auto              twoPiN = TWO_PI * n;
        return static_cast<FloatType>(0.5) *
               (static_cast<FloatType>(1) - std::cos(twoPiN / sizeMinusOne));
    }
};

/*! \brief Hamming window generator. */
class HammingGenerator final
{
public:
    /*!
     * \brief Generate a vector of window coefficients.
     * \param[in,out] windowCoeffs - Pre-sized vector to receive window coefficients.
     */
    template <typename FloatType> void operator()(std::vector<FloatType>& windowCoeffs) const
    {
        WindowGenerator(windowCoeffs,
                        std::bind(&HammingGenerator::EvaluateCoefficient<FloatType>,
                                  std::placeholders::_1,
                                  std::placeholders::_2));
    }

private:
    /*!
     * \brief Evaluate nth window coefficient.
     * \param[in] n - The current window coefficient index as a floating point type T.
     * \param[in] sizeMinusOne - The window length minus 1 as a floating point type T.
     * \return The nth window coefficient as a floating type T.
     */
    template <typename FloatType>
    static FloatType EvaluateCoefficient(FloatType n, FloatType sizeMinusOne)
    {
        static const auto a0     = static_cast<FloatType>(0.53836);
        static const auto a1     = static_cast<FloatType>(0.46164);
        static const auto TWO_PI = TwoPi<FloatType>();
        auto              twoPiN = TWO_PI * n;
        return a0 - (a1 * std::cos(twoPiN / sizeMinusOne));
    }
};

/*! \brief Rectangle window generator. */
struct RectangleGenerator final
{
    /*!
     * \brief Generate a vector of window coefficients.
     * \param[in,out] windowCoeffs - Pre-sized vector to receive window coefficients.
     */
    template <typename FloatType> void operator()(std::vector<FloatType>& windowCoeffs) const
    {
        auto size = windowCoeffs.size();
        DSP_ASSERT_THROW(size > 1, "invalid window size");
        std::fill(windowCoeffs.begin(), windowCoeffs.end(), static_cast<FloatType>(1));
    }
};

/*! \brief Bartlett window generator. */
class BartlettGenerator final
{
public:
    /*!
     * \brief Generate a vector of window coefficients.
     * \param[in,out] windowCoeffs - Pre-sized vector to receive window coefficients.
     */
    template <typename FloatType> void operator()(std::vector<FloatType>& windowCoeffs) const
    {
        WindowGenerator(windowCoeffs,
                        std::bind(&BartlettGenerator::EvaluateCoefficient<FloatType>,
                                  std::placeholders::_1,
                                  std::placeholders::_2));
    }

private:
/*!
     * \brief Evaluate nth window coefficient.
     * \param[in] n - The current window coefficient index as a floating point type T.
     * \param[in] sizeMinusOne - The window length minus 1 as a floating point type T.
     * \return The nth window coefficient as a floating type T.
     */
    template <typename FloatType>
    static FloatType EvaluateCoefficient(FloatType n, FloatType sizeMinusOne)
    {
        auto commonTerm = sizeMinusOne / 2;
        auto numerator  = n - commonTerm;
        return static_cast<FloatType>(1) - std::abs(numerator / commonTerm);
    }
};

/*! \brief Exact Blackman window generator. */
class ExactBlackmanGenerator final
{
public:
    /*!
     * \brief Generate a vector of window coefficients.
     * \param[in,out] windowCoeffs - Pre-sized vector to receive window coefficients.
     */
    template <typename FloatType> void operator()(std::vector<FloatType>& windowCoeffs) const
    {
        WindowGenerator(windowCoeffs,
                        std::bind(&ExactBlackmanGenerator::EvaluateCoefficient<FloatType>,
                                  std::placeholders::_1,
                                  std::placeholders::_2));
    }

private:
    /*!
     * \brief Evaluate nth window coefficient.
     * \param[in] n - The current window coefficient index as a floating point type T.
     * \param[in] sizeMinusOne - The window length minus 1 as a floating point type T.
     * \return The nth window coefficient as a floating type T.
     */
    template <typename FloatType>
    static FloatType EvaluateCoefficient(FloatType n, FloatType sizeMinusOne)
    {
        static const auto a0      = static_cast<FloatType>(7938. / 18608.);
        static const auto a1      = static_cast<FloatType>(9240. / 18608.);
        static const auto a2      = static_cast<FloatType>(1430. / 18608.);
        static const auto TWO_PI  = TwoPi<FloatType>();
        static const auto FOUR_PI = static_cast<FloatType>(2) * TwoPi<FloatType>();
        auto              twoPiN  = TWO_PI * n;
        auto              fourPiN = FOUR_PI * n;
        return a0 - (a1 * std::cos(twoPiN / sizeMinusOne)) +
               (a2 * std::cos(fourPiN / sizeMinusOne));
    }
};

/*! \brief Blackman window generator. */
class BlackmanGenerator final
{
public:
    /*!
     * \brief Generate a vector of window coefficients.
     * \param[in,out] windowCoeffs - Pre-sized vector to receive window coefficients.
     */
    template <typename FloatType> void operator()(std::vector<FloatType>& windowCoeffs) const
    {
        WindowGenerator(windowCoeffs,
                        std::bind(&BlackmanGenerator::EvaluateCoefficient<FloatType>,
                                  std::placeholders::_1,
                                  std::placeholders::_2));
    }

private:
    /*!
     * \brief Evaluate nth window coefficient.
     * \param[in] n - The current window coefficient index as a floating point type T.
     * \param[in] sizeMinusOne - The window length minus 1 as a floating point type T.
     * \return The nth window coefficient as a floating type T.
     */
    template <typename FloatType>
    static FloatType EvaluateCoefficient(FloatType n, FloatType sizeMinusOne)
    {
        static const auto a0      = static_cast<FloatType>(0.42);
        static const auto a1      = static_cast<FloatType>(0.5);
        static const auto a2      = static_cast<FloatType>(0.08);
        static const auto TWO_PI  = TwoPi<FloatType>();
        static const auto FOUR_PI = static_cast<FloatType>(2) * TwoPi<FloatType>();
        auto              twoPiN  = TWO_PI * n;
        auto              fourPiN = FOUR_PI * n;
        return a0 - (a1 * std::cos(twoPiN / sizeMinusOne)) +
               (a2 * std::cos(fourPiN / sizeMinusOne));
    }
};

/*! \brief Kaiser window generator. */
class KaiserGenerator final
{
public:
    /*!
     * \brief Generate Kaiser window coefficients.
     * \param[in] beta - Controls side-lobe roll-off, where beta == Pi*alpha.
     */
    explicit KaiserGenerator(double beta)
        : m_beta(beta)
    {
        DSP_ASSERT_THROW(m_beta > 0, "beta <= 0");
    }
    /*!
     * \brief Generate a vector of window coefficients.
     * \param[in,out] windowCoeffs - Pre-sized vector to receive window coefficients.
     */
    template <typename FloatType> void operator()(std::vector<FloatType>& windowCoeffs) const
    {
        WindowGenerator(windowCoeffs,
                        std::bind(&KaiserGenerator::EvaluateCoefficient<FloatType>,
                                  this,
                                  std::placeholders::_1,
                                  std::placeholders::_2));
    }

private:
    /*!
     * \brief Evaluate nth window coefficient.
     * \param[in] n - The current window coefficient index as a floating point type T.
     * \param[in] sizeMinusOne - The window length minus 1 as a floating point type T.
     * \return The nth window coefficient as a floating type T.
     */
    template <typename FloatType>
    FloatType EvaluateCoefficient(FloatType n, FloatType sizeMinusOne) const
    {
        auto term = ((static_cast<FloatType>(2) * n) / sizeMinusOne) - static_cast<FloatType>(1);
        auto arg = static_cast<FloatType>(m_beta) * sqrt(static_cast<FloatType>(1) - (term * term));
        return Bessel(arg) / Bessel(static_cast<FloatType>(m_beta));
    }

private:
    /*! \brief beta controls side-lobe roll-off, where beta == Pi*alpha. */
    double m_beta{1.};
};

/*! \brief Lanczos (sinc) window generator. */
class LanczosGenerator final
{
public:
    /*!
     * \brief Generate a vector of window coefficients.
     * \param[in,out] windowCoeffs - Pre-sized vector to receive window coefficients.
     */
    template <typename FloatType> void operator()(std::vector<FloatType>& windowCoeffs) const
    {
        WindowGenerator(windowCoeffs,
                        std::bind(&LanczosGenerator::EvaluateCoefficient<FloatType>,
                                  std::placeholders::_1,
                                  std::placeholders::_2));
    }

private:
    /*!
     * \brief Evaluate nth window coefficient.
     * \param[in] n - The current window coefficient index as a floating point type T.
     * \param[in] sizeMinusOne - The window length minus 1 as a floating point type T.
     * \return The nth window coefficient as a floating type T.
     */
    template <typename FloatType>
    static FloatType EvaluateCoefficient(FloatType n, FloatType sizeMinusOne)
    {
        auto arg = ((static_cast<FloatType>(2) * n) / sizeMinusOne) - static_cast<FloatType>(1);
        return SincNorm(arg);
    }
};

/*!
 * \brief Top-level "Window Function" class use this class instead of the above helper classes.
 *
 * Template typename T defines the window coefficient data type.
 */
template <typename FloatType> class WindowFunction final
{
public:
    /*! \brief Default constructor. */
    WindowFunction() = default;
    /*! \brief Destructor. */
    ~WindowFunction() = default;
    /*! \brief Default copy constructor. */
    WindowFunction(WindowFunction const&) = default;
    /*! \brief Default move constructor. */
    WindowFunction(WindowFunction&&) = default;
    /*! \brief Default copy assignment operator. */
    WindowFunction& operator=(WindowFunction const&) = default;
    /*! \brief Default move assignment operator. */
    WindowFunction& operator=(WindowFunction&&) = default;

    /*!
     * \brief Initialisation constructor.
     * \param[in] generator - Window function generator.
     * \param[in] size - Size required for the window.
     * \param[in] ignoreLastValue - If size is odd and window is to be applied to data for FFT
     * processing set this to true, else set to false.
     */
    template <typename Generator>
    WindowFunction(Generator const& generator, size_t size, bool ignoreLastValue)
        : m_windowCoefficients(size, static_cast<FloatType>(1))
        , m_ignoreLastValue(ignoreLastValue && (size % 2 == 1))
        , m_effectiveSize(ignoreLastValue ? m_windowCoefficients.size() - 1
                                          : m_windowCoefficients.size())
    {
        static_assert(std::is_floating_point<FloatType>::value, "invalid floating point type");
        generator(m_windowCoefficients);
        ComputeGains();
    }

    /*!
     * \brief Initialisation method.
     * \param[in] generator - Window function generator.
     * \param[in] size - Size required for the window.
     * \param[in] ignoreLastValue - If size is odd and window is to be applied to data for FFT
     * processing set this to true, else set to false.
     */
    template <typename Generator>
    void Initialise(Generator const& generator, size_t size, bool ignoreLastValue)
    {
        *this = std::move(WindowFunction(generator, size, ignoreLastValue));
    }

    /*!
     * \brief Get the coherent gain of the window coefficients.
     * \return The value as a floating point type T.
     */
    FloatType CoherentGain() const
    {
        return m_coherentGain;
    }

    /*!
     * \brief Get the power gain of the window coefficients.
     * \return The value as a floating point type T.
     */
    FloatType PowerGain() const
    {
        return m_powerGain;
    }

    /*!
     * \brief Get the combined gain of the window coefficients.
     * \return The value as a floating point type T.
     */
    FloatType CombinedGain() const
    {
        return m_coherentGain * m_powerGain;
    }

    /*!
     * \brief Apply gain correction to data range.
     * \param[in] dataFirst - Iterator to the first data sample to scale by the window.
     * \param[in] dataLast - Iterator to the one past the final data sample to scale by the window.
     * \param[in,out] resultFirst - Iterator the start of the output container to receive the scaled
     * data.
     * \param[in] gain - Gain to correct data for.
     */
    template <typename InIter, typename OutIter>
    static void ApplyGainCorrection(InIter dataFirst, InIter dataLast, OutIter resultFirst,
                                    FloatType gain)
    {
        for (auto itr = dataFirst; itr != dataLast;
             std::advance(itr, 1), std::advance(resultFirst, 1))
        {
            *resultFirst = *itr / gain;
        }
    }

    /*!
     * \brief Get the effective noise bandwidth of the window coefficients.
     * \return The value as a floating point type T.
     */
    FloatType EffectiveNoiseBandwidth() const
    {
        return m_enbw;
    }

    /*!
     * \brief Get the actual length the window.
     * \return The number of window coefficients.
     */
    size_t ActualSize() const
    {
        return m_windowCoefficients.size();
    }

    /*!
     * \brief Get the effective length the window.
     * \return The number of window coefficients to be used.
     *
     * When the window has sodd length and is to be used for FFT processing
     * this will be 1 less than the actual size.
     */
    size_t EffectiveSize() const
    {
        return m_effectiveSize;
    }

    /*!
     * \brief Get a copy of the window coefficients.
     * \return The number of window coefficients.
     */
    std::vector<FloatType> Coefficients() const
    {
        if (m_windowCoefficients.size() == m_effectiveSize)
        {
            return m_windowCoefficients;
        }
        else
        {
            return std::vector<FloatType>(
                m_windowCoefficients.begin(),
                std::next(m_windowCoefficients.begin(), static_cast<int>(m_effectiveSize)));
        }
    }

    /*!
     * \brief Apply the window coefficients to a block of data.
     * \param[in] dataFirst - Iterator to the first data sample to scale by the window.
     * \param[in] dataLast - Iterator to the one past the final data sample to scale by the window.
     * \param[in,out] resultFirst - Iterator the start of the output container to receive the scaled
     * data.
     *
     * This function can be used safely for in-place trasnformations.
     */
    template <typename InIter, typename OutIter>
    void operator()(InIter dataFirst, InIter dataLast, OutIter resultFirst) const
    {
        auto inputLength = std::distance(dataFirst, dataLast);
        DSP_ASSERT_THROW(inputLength == static_cast<decltype(inputLength)>(m_effectiveSize),
                         "invalid data size");

        // Generic lambda function to allow multiplying mixed types.
        // This is because we may want to window complex data not just
        // real valued data.
        auto multiply = [](auto const& x, auto const& y) { return x * y; };

        // Transform the data by applying the window coefficients using
        // piecewise multiplication.
        std::transform(dataFirst, dataLast, m_windowCoefficients.begin(), resultFirst, multiply);
    }

private:
    /*! \brief Compute the window gains. */
    void ComputeGains()
    {
        auto size = EffectiveSize();

        for (size_t i = 0; i < size; ++i)
        {
            m_coherentGain += m_windowCoefficients[i];
            m_enbw += m_windowCoefficients[i] * m_windowCoefficients[i];
        }

        auto enbwDivisor = m_coherentGain * m_coherentGain;

        if (std::abs(enbwDivisor) > static_cast<FloatType>(1.e-9))
        {
            m_enbw = (static_cast<FloatType>(size) * m_enbw) / enbwDivisor;
        }

        m_coherentGain /= static_cast<FloatType>(size);
        m_powerGain = m_coherentGain * m_coherentGain * m_enbw;
    }

private:
    /*! \brief Vector of window coefficients. */
    std::vector<FloatType> m_windowCoefficients{};
    /*! \brief Flag to define if we need to ignore the last window coefficient. */
    bool m_ignoreLastValue{false};
    /*! \brief The effective size of the window. */
    size_t m_effectiveSize{0};
    /*! \brief The coeherent gain of the window. */
    FloatType m_coherentGain{0.};
    /*! \brief The power gain of the window. */
    FloatType m_powerGain{0.};
    /*! \brief The effective noise bandwidth of the window. */
    FloatType m_enbw{0.};
};

/*! \brief Convenience typedef to WindowFunction<float>. */
using window_fn_f = WindowFunction<float>;
/*! \brief Convenience typedef to WindowFunction<double>. */
using window_fn_d = WindowFunction<double>;
/*! \brief Convenience typedef to WindowFunction<long double>. */
using window_fn_ld = WindowFunction<long double>;

} // namespace dsp

#endif // DSP_WINDOW_FUNCTIONS_HPP
