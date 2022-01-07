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
 * \file dsp_fft.hpp
 * \brief File containing generic FFT algorithms.
 */

#ifndef DSP_FFT_HPP
#define DSP_FFT_HPP

#include <complex>
#include <algorithm>
#include <functional>
#include <cstdint>
#include "dsp_roots.hpp"
#include "dsp_window_functions.hpp"

/*! \brief dsp namespace */
namespace dsp
{

/*! \brief Class to perform FFTs on complex data */
template <typename FloatType> class ComplexFFT final
{
public:
    /*! \brief Typedef to complex floating point value */
    using complex_type = std::complex<FloatType>;
    /*! \brief Typedef to complex floating point vector */
    using complex_vector = std::vector<complex_type>;
    /*! \brief Typedef to real floating point vector */
    using real_vector = std::vector<FloatType>;

public:
    /*! \brief Default constructor. */
    ComplexFFT()
    {
        static_assert(std::is_floating_point<FloatType>::value,
                      "FloatType must be either float, double or long double");
    }

    /*!
     * \brief Perform a forward FFT on the input range
     * \param[in] first - iterator to start of data samples.
     * \param[in] last - iterator to one past the end of the data samples.
     *
     * Input data range must be implicitly convertible to std::complex<FloatType>.
     */
    template <typename Iter> static complex_vector Forward(Iter first, Iter last)
    {
        auto N = static_cast<size_t>(std::distance(first, last));

        DSP_ASSERT_THROW(IsPowerOf2(N), "FFT size not a power of 2");

        complex_vector data(first, last);

        // Forward FFT of data.
        CooleyTukeyFFT(data);

        return data;
    }

    /*!
     * \brief Perform a forward FFT on complex data vector in-place.
     * \param[in] data - vector of complex data.
     */
    static void Forward(complex_vector& data)
    {
        DSP_ASSERT_THROW(IsPowerOf2(data.size()), "FFT size not a power of 2");

        // Forward FFT of data.
        CooleyTukeyFFT(data);
    }

    /*!
     * \brief Normalise complex FFT spectrum in-place.
     * \param[in] cplxFft - vector of complex FFT data.
     *
     * Normalise a complex FFT result by dividing
     * real and complex parts by the FFT size.
     */
    static void Normalise(complex_vector& cplxFft)
    {
        auto N = static_cast<FloatType>(cplxFft.size());

        for (auto& z : cplxFft)
        {
            z /= N;
        }
    }

    /*!
     * \brief Denormalise complex FFT spectrum in-place.
     * \param[in] cplxFft - vector of normalised complex FFT data.
     *
     * Denormalise a normalised complex FFT result by multiplying
     * real and complex parts by the FFT size.
     */
    static void Denormalise(complex_vector& cplxFft)
    {
        auto N = static_cast<FloatType>(cplxFft.size());

        for (auto& z : cplxFft)
        {
            z *= N;
        }
    }

    /*!
     * \brief Perform an inverse of complex FFT data in-place.
     * \param[in] cplxFft - complex data vector containing denormalised complex FFT.
     */
    static void Inverse(complex_vector& cplxFft)
    {
        DSP_ASSERT_THROW(IsPowerOf2(cplxFft.size()), "FFT size not a power of 2");

        // Conjugation lambda designed for speed.
        auto conjugate = [](complex_vector& cplxVector) {
            static constexpr auto minusOne = static_cast<FloatType>(-1);

            for (auto& z : cplxVector)
            {
                z.imag(z.imag() * minusOne);
            }
        };

        // Conjugate the FFT BIN values.
        conjugate(cplxFft);

        // Forward FFT of conjugates.
        CooleyTukeyFFT(cplxFft);

        // Conjugate the values again.
        conjugate(cplxFft);

        // Scale the numbers
        Normalise(cplxFft);
    }

    /*!
     * \brief Convert complex FFT data to a magnitude only spectrum in-place.
     * \param[in] cplxFft - complex data vector containing normalised complex FFT.
     * \param[in] zeroUnused - zero unused BINs in complex vector.
     * \param[in] fullSpectrum - keep the full spectrum.
     *
     * After calling this only the real part of each complex value
     * is valid, the imaginary part will be 0. Also because the power
     * is spread over a positive and negative half spectra we multiply
     * the value in each BIN by 2 to get the correct amplitude. The
     * negative half of the FFT is not needed unless we want the full
     * spectrum so if we only want the half spectrum we use so only the
     * first N/2 values.
     */
    static void ToMagnitude(complex_vector& cplxFft, bool zeroUnused = false,
                            bool fullSpectrum = false)
    {
        static constexpr auto two      = static_cast<FloatType>(2);
        static constexpr auto zero     = static_cast<FloatType>(0);
        auto                  halfSize = fullSpectrum ? cplxFft.size() : cplxFft.size() / 2;
        auto                  endItr   = std::next(cplxFft.begin(), static_cast<int32_t>(halfSize));

        for (auto zItr = cplxFft.begin(); zItr != endItr; ++zItr)
        {
            if (zItr != cplxFft.begin())
            {
                // This BIN does not contain the DC so
                // we multiply this by scalar.
                *zItr *= two;
            }

            auto magnitude = std::abs(*zItr);
            zItr->real(magnitude);
            zItr->imag(zero);
        }

        if (zeroUnused)
        {
            for (auto zItr = endItr; zItr != cplxFft.end(); ++zItr)
            {
                *zItr = {zero, zero};
            }
        }
    }

    /*!
     * \brief Convert complex FFT data to a magnitude only spectrum.
     * \param[in] cplxFft - complex data vector containing normalised complex FFT.
     * \param[out] magFft - The magnitude only spectrum as a real valued vector.
     * \param[in] fullSpectrum - keep the full spectrum.
     */
    static void ToMagnitude(complex_vector const& cplxFft, real_vector& magFft,
                            bool fullSpectrum = false)
    {
        static constexpr auto two      = static_cast<FloatType>(2);
        auto                  halfSize = fullSpectrum ? cplxFft.size() : cplxFft.size() / 2;
        auto                  endItr   = std::next(cplxFft.begin(), static_cast<int32_t>(halfSize));
        magFft.resize(halfSize);
        size_t index = 0;

        for (auto zItr = cplxFft.begin(); zItr != endItr; ++zItr, ++index)
        {
            if (zItr != cplxFft.begin())
            {
                // This BIN does not contain the DC so
                // we multiply this by scalar.
                magFft[index] = std::abs(*zItr * two);
            }
            else
            {
                magFft[index] = std::abs(*zItr);
            }
        }
    }

    /*!
     * \brief Convert complex FFT data to power spectrum in-place.
     * \param[in] cplxFft - complex data vector containing normalised complex FFT.
     * \param[in] zeroUnused - zero unused BINs in complex vector.
     * \param[in] fullSpectrum - keep the full spectrum.
     *
     * The power spectrum is computed as re[n]^2 + im[n]^2. Only the real
     * part of each complex value is valid, the imaginary part will be 0.
     * The negative half of the FFT is not needed unless we want the full
     * spectrum so if we only want the half spectrum we use so only the
     * first N/2 values.
     */
    static void ToPower(complex_vector& cplxFft, bool zeroUnused = false, bool fullSpectrum = false)
    {
        static constexpr auto zero     = static_cast<FloatType>(0);
        auto                  halfSize = fullSpectrum ? cplxFft.size() : cplxFft.size() / 2;
        auto                  endItr   = std::next(cplxFft.begin(), static_cast<int32_t>(halfSize));

        for (auto zItr = cplxFft.begin(); zItr != endItr; ++zItr)
        {
            auto power = std::norm(*zItr);
            zItr->real(power);
            zItr->imag(zero);
        }

        if (zeroUnused)
        {
            for (auto zItr = endItr; zItr != cplxFft.end(); ++zItr)
            {
                *zItr = {zero, zero};
            }
        }
    }

    /*!
     * \brief Convert complex FFT data to power spectrum.
     * \param[in] cplxFft - complex data vector containing normalised complex FFT.
     * \param[out] powerSpectrum - The power spectrum as a real valued vector.
     * \param[in] fullSpectrum - keep the full spectrum.
     */
    static void ToPower(complex_vector const& cplxFft, real_vector& powerSpectrum,
                        bool fullSpectrum = false)
    {
        auto halfSize = fullSpectrum ? cplxFft.size() : cplxFft.size() / 2;
        auto endItr   = std::next(cplxFft.begin(), static_cast<int32_t>(halfSize));
        powerSpectrum.resize(halfSize);
        size_t index = 0;

        for (auto zItr = cplxFft.begin(); zItr != endItr; ++zItr, ++index)
        {
            powerSpectrum[index] = std::norm(*zItr);
        }
    }

    /*!
     * \brief Convert power spectrum to PSD (power spectral density) in-place.
     * \param[in] powerSpectrum - complex data vector containing the power spectrum.
     * \param[in] binWidthHz - the bin width in Hz of the spectral data.
     * \param[in] zeroUnused - zero unused BINs in complex vector.
     * \param[in] fullSpectrum - keep the full spectrum.
     *
     * The real valued power spectrum stored in a complex vector
     * is converted to PSD values by dividing by BIN width. After
     * calling this function the half powerSpectrum contains the PSD in
     * its real components. The negative half of the FFT is not needed
     * unless we want the full spectrum so if we only want the half spectrum
     * we use so only the first N/2 values.
     */
    static void ToPsd(complex_vector& powerSpectrum, FloatType binWidthHz, bool zeroUnused = false,
                      bool fullSpectrum = false)
    {
        static constexpr auto zero = static_cast<FloatType>(0);
        auto halfSize              = fullSpectrum ? powerSpectrum.size() : powerSpectrum.size() / 2;
        auto endItr = std::next(powerSpectrum.begin(), static_cast<int32_t>(halfSize));

        for (auto psIter = powerSpectrum.begin(); psIter != endItr; std::advance(psIter, 1))
        {
            auto psd = psIter->real() / binWidthHz;
            psIter->real(psd);
        }

        if (zeroUnused)
        {
            for (auto psIter = endItr; psIter != powerSpectrum.end(); ++psIter)
            {
                *psIter = {zero, zero};
            }
        }
    }

    /*!
     * \brief Convert power spectrum to PSD (power spectral density) in-place.
     * \param[in] powerSpectrum - real valued data vector containing the power spectrum.
     * \param[in] binWidthHz - the bin width in Hz of the spectral data.
     *
     * After this function is called the powerSpectrum will contain the PSD values.
     */
    static void ToPsd(real_vector& powerSpectrum, FloatType binWidthHz)
    {
        static constexpr auto one = static_cast<FloatType>(1);

        for (auto& v : powerSpectrum)
        {
            v /= binWidthHz;
        }
    }

    /*!
     * \brief Convert power spectrum to PSD (power spectral density).
     * \param[in] powerSpectrum - complex data vector containing the power spectrum.
     * \param[in] binWidthHz - the bin width in Hz of the spectral data.
     * \param[out] psd - the real valued vector containing the PSD values.
     * \param[in] fullSpectrum - keep the full spectrum.
     *
     * The real valued power spectrum stored in a complex vector is converted to
     * PSD values by dividing by BIN width. The negative half of the FFT is not
     * needed  unless we want the full spectrum so if we only want the half spectrum
     * we use so only the first N/2 values.
     */
    static void ToPsd(complex_vector const& powerSpectrum, FloatType binWidthHz, real_vector& psd,
                      bool fullSpectrum = false)
    {
        auto halfSize = fullSpectrum ? powerSpectrum.size() : powerSpectrum.size() / 2;
        auto endItr   = std::next(powerSpectrum.begin(), static_cast<int32_t>(halfSize));
        psd.resize(halfSize);
        size_t index = 0;

        for (auto psIter = powerSpectrum.begin(); psIter != endItr;
             std::advance(psIter, 1), ++index)
        {
            psd[index] = psIter->real() / binWidthHz;
        }
    }

    /*!
     * \brief Convert power spectrum to PSD (power spectral density).
     * \param[in] powerSpectrum - real valued data vector containing te power spectrum.
     * \param[in] binWidthHz - the bin width in Hz of the spectral data.
     * \param[out] psd - the real valued vector containing the PSD values.
     *
     * The real valued power spectrum stored in a real valued vector
     * is converted to PSD values by dividing by BIN width.
     */
    static void ToPsd(const real_vector& powerSpectrum, FloatType binWidthHz, real_vector& psd)
    {
        psd.resize(powerSpectrum.size());
        size_t index = 0;

        for (auto psIter = powerSpectrum.begin(); psIter != powerSpectrum.end();
             std::advance(psIter, 1), ++index)
        {
            psd[index] = *psIter / binWidthHz;
        }
    }

    /*!
     * \brief Convert power spectrum to 3 BIN summed magnitude spectrum in-place.
     * \param[in] powerSpectrum - complex data vector containing the power spectrum.
     * \param[in] zeroUnused - zero unused BINs in complex vector.
     * \param[in] fullSpectrum - keep the full spectrum.
     *
     * The real valued power spectrum stored in a complex vector
     * is converted to 3-BIN summed magnitude spectrum. After
     * calling this function the powerSpectrum contains the result in
     * its real components. The negative half of the FFT is not needed
     * unless we want the full spectrum so if we only want the half spectrum
     * we use so only the first N/2 values.
     */
    static void To3BinSum(complex_vector& powerSpectrum, bool zeroUnused = false,
                          bool fullSpectrum = false)
    {
        static constexpr auto       one     = static_cast<FloatType>(1);
        static constexpr auto       zero    = static_cast<FloatType>(0);
        static const auto           sqrt2   = dsp::SqrtTwo<FloatType>();
        static const complex_vector convVec = {{one, zero}, {one, zero}, {one, zero}};
        auto           halfSize = fullSpectrum ? powerSpectrum.size() : powerSpectrum.size() / 2;
        auto           endItr   = std::next(powerSpectrum.begin(), static_cast<int32_t>(halfSize));
        complex_vector convRes(halfSize + convVec.size() - 1);

        Convolve(powerSpectrum.begin(), endItr, convVec.begin(), convVec.end(), convRes.begin());

        auto resItr = std::next(convRes.begin());
        auto resEnd = std::next(resItr, static_cast<int32_t>(halfSize));
        auto outItr = powerSpectrum.begin();

        while (resItr != resEnd)
        {
            // Convert from RMS power to peak magnitude.
            outItr->real(std::sqrt(resItr->real()) * sqrt2);
            std::advance(outItr, 1);
            std::advance(resItr, 1);
        }

        if (zeroUnused)
        {
            for (auto psIter = endItr; psIter != powerSpectrum.end(); ++psIter)
            {
                *psIter = {zero, zero};
            }
        }
    }

    /*!
     * \brief Convert power spectrum to 3 BIN summed magnitude spectrum in-place.
     * \param[in] powerSpectrum - real valued data vector containing the power spectrum.
     *
     * After this function is called the powerSpectrum will contain the 3 BIN summed
     * magnitude values.
     */
    static void To3BinSum(real_vector& powerSpectrum)
    {
        static constexpr auto    one     = static_cast<FloatType>(1);
        static const real_vector convVec = {one, one, one};
        auto                     sqrt2   = dsp::SqrtTwo<FloatType>();
        auto                     size    = powerSpectrum.size();
        real_vector              convRes(size + convVec.size() - 1);

        Convolve(powerSpectrum.begin(),
                 powerSpectrum.end(),
                 convVec.begin(),
                 convVec.end(),
                 convRes.begin());

        auto resBegin = std::next(convRes.begin());
        auto resEnd   = std::next(resBegin, static_cast<int32_t>(size));

        // Convert from RMS power to peak magnitude.
        auto transformer = [sqrt2](auto x) { return std::sqrt(x) * sqrt2; };

        std::transform(resBegin, resEnd, powerSpectrum.begin(), transformer);
    }

    /*!
     * \brief Convert power spectrum to 3 BIN summed magnitude spectrum.
     * \param[in] powerSpectrum - complex data vector containing the power spectrum.
     * \param[out] magFft - the real valued vector containing the 3 BIN summed values.
     * \param[in] fullSpectrum - keep the full spectrum.
     *
     * The real valued power spectrum stored in a complex vector is converted to
     * 3 BIN summed magnitude values. The negative half of the FFT is not needed
     * unless we want the full spectrum so if we only want the half spectrum
     * we use so only the first N/2 values.
     */
    static void To3BinSum(complex_vector const& powerSpectrum, real_vector& magFft,
                          bool fullSpectrum = false)
    {
        static constexpr auto       one     = static_cast<FloatType>(1);
        static constexpr auto       zero    = static_cast<FloatType>(0);
        static const complex_vector convVec = {{one, zero}, {one, zero}, {one, zero}};
        auto                        sqrt2   = dsp::SqrtTwo<FloatType>();
        auto           halfSize = fullSpectrum ? powerSpectrum.size() : powerSpectrum.size() / 2;
        auto           endItr   = std::next(powerSpectrum.begin(), static_cast<int32_t>(halfSize));
        complex_vector convRes(halfSize + convVec.size() - 1);
        magFft.resize(halfSize);

        Convolve(powerSpectrum.begin(), endItr, convVec.begin(), convVec.end(), convRes.begin());

        auto resBegin = std::next(convRes.begin());
        auto resEnd   = std::next(resBegin, halfSize);

        auto transformer = [sqrt2](auto const& z) { return std::sqrt(z.real()) * sqrt2; };

        std::transform(resBegin, resEnd, magFft.begin(), transformer);
    }

    /*!
     * \brief Convert power spectrum to 3 BIN summed magnitude spectrum.
     * \param[in] powerSpectrum - real valued data vector containing the power spectrum.
     * \param[out] magFft - the real valued vector containing the 3 BIN summed values.
     *
     * The real valued power spectrum stored in a real valued vector
     * is converted to 3 BIN summed magnitude values.
     */
    static void To3BinSum(real_vector const& powerSpectrum, real_vector& magFft)
    {
        static constexpr auto    one     = static_cast<FloatType>(1);
        static const real_vector convVec = {one, one, one};
        auto                     sqrt2   = dsp::SqrtTwo<FloatType>();
        auto                     size    = powerSpectrum.size();
        magFft.resize(size + convVec.size() - 1);

        Convolve(powerSpectrum.begin(),
                 powerSpectrum.end(),
                 convVec.begin(),
                 convVec.end(),
                 magFft.begin());

        std::rotate(magFft.begin(), std::next(magFft.begin()), magFft.end());
        magFft.resize(size);

        // Convert from RMS power to peak magnitude.
        auto transformer = [sqrt2](auto x) { return std::sqrt(x) * sqrt2; };

        std::transform(magFft.begin(), magFft.end(), magFft.begin(), transformer);
    }

private:
    /*!
     * \brief Perform Cooley-Tukey algorithm to compute the FFT spectrum in-place.
     * \param[in] data - complex vector of data to have FFT perfromed on it.
     *
     * After this function is the data vector contains the FFT spectrum.
     */
    static void CooleyTukeyFFT(complex_vector& data)
    {
        static constexpr auto one = static_cast<FloatType>(1);
        static const auto     pi  = Pi<FloatType>();

        // DFT.
        auto N      = static_cast<uint32_t>(data.size());
        auto fN     = static_cast<FloatType>(N);
        auto k      = N;
        auto thetaT = pi / fN;
        auto phiT   = complex_type(std::cos(thetaT), std::sin(thetaT));

        while (k > 1)
        {
            auto n = k;
            k >>= 1;
            phiT   = phiT * phiT;
            auto T = complex_type(one);

            for (uint32_t l = 0; l < k; ++l)
            {
                for (uint32_t a = l; a < N; a += n)
                {
                    uint32_t b = a + k;
                    auto     t = data[a] - data[b];
                    data[a] += data[b];
                    data[b] = t * T;
                }

                T *= phiT;
            }
        }

        // Decimate.
        auto m = static_cast<uint32_t>(std::log2(fN));

        for (uint32_t a = 0; a < N; ++a)
        {
            auto b = a;

            // Reverse bits.
            b = (((b & 0xaaaaaaaa) >> 1) | ((b & 0x55555555) << 1));
            b = (((b & 0xcccccccc) >> 2) | ((b & 0x33333333) << 2));
            b = (((b & 0xf0f0f0f0) >> 4) | ((b & 0x0f0f0f0f) << 4));
            b = (((b & 0xff00ff00) >> 8) | ((b & 0x00ff00ff) << 8));
            b = ((b >> 16) | (b << 16)) >> (32 - m);

            if (b > a)
            {
                std::swap(data[a], data[b]);
            }
        }
    }
};

/*! \brief Class to generate real-valued FFT spectrum using 3 BIN summing. */
template <typename FloatType> class ThreeBinSumFft final
{
    /*! \brief Typedef to window function. */
    using window_fn = WindowFunction<FloatType>;
    /*! \brief Typedef to complex floating point value */
    using complex_type = std::complex<FloatType>;
    /*! \brief Typedef to complex floating point vector */
    using complex_vector = std::vector<complex_type>;
    /*! \brief Typedef to real floating point vector */
    using real_vector = std::vector<FloatType>;
    /*! \brief Convenience typedef to ComplexFFT<FloatType>. */
    using complex_fft = ComplexFFT<FloatType>;

public:
    /*! \brief Default constructor. */
    ThreeBinSumFft() = default;
    /*! \brief Destructor. */
    ~ThreeBinSumFft() = default;
    /*! \brief Default copy constructor. */
    ThreeBinSumFft(ThreeBinSumFft const&) = default;
    /*! \brief Default move constructor. */
    ThreeBinSumFft(ThreeBinSumFft&&) = default;
    /*! \brief Default copy assignment operator. */
    ThreeBinSumFft& operator=(ThreeBinSumFft const&) = default;
    /*! \brief Default move assignment operator. */
    ThreeBinSumFft& operator=(ThreeBinSumFft&&) = default;

    /*!
     * \brief Initialisation constructor.
     * \param[in] generator - the window function generator functor.
     * \param[in] fftSize - the number of samples used by the FFT.
     */
    template <typename Generator>
    ThreeBinSumFft(Generator const& generator, size_t fftSize)
        : m_windowFunction(generator, fftSize + 1, true)
        , m_workspace(fftSize)
    {
        DSP_ASSERT_THROW(IsPowerOf2(fftSize), "FFT size not a power of 2");
    }

    /*!
     * \brief Initialiser method.
     * \param[in] generator - the window function generator functor.
     * \param[in] fftSize - the number of samples used by the FFT.
     */
    template <typename Generator> void Initialise(Generator const& generator, size_t fftSize)
    {
        *this = std::move(ThreeBinSumFft(generator, fftSize));
    }

    /*!
     * \brief Perform the FFT on the real-valued signal data samples.
     * \param[in] signalFirst - iterator to first sample of signal data.
     * \param[in] signalLast - iterator to one past the last sample of signal data.
     * \param[out] realSpectra - real-valued vector to hold result spectrum.
     * \param[in] fullSpectrum - keep the full spectrum.
     * \param[out] phases - (optional) bin by bin spectrum phases.
     */
    template <typename Iter>
    void operator()(Iter signalFirst, Iter signalLast, real_vector& realSpectra,
                    bool fullSpectrum = false, real_vector* phases = nullptr)
    {
        auto signalLength = std::distance(signalFirst, signalLast);
        DSP_ASSERT_THROW(signalLength == static_cast<decltype(signalLength)>(m_workspace.size()),
                         "signal length is incorrect");

        // Window the data.
        m_windowFunction(signalFirst, signalLast, realSpectra.begin());

        // Copy data into workspace.
        auto convertToComplex = [](auto x) { return complex_type{x, 0}; };
        std::transform(
            realSpectra.begin(), realSpectra.end(), m_workspace.begin(), convertToComplex);

        // Perform the forward FFT on the windowed data.
        complex_fft::Forward(m_workspace);

        // Typically you would normalise the FFT spectrum at this stage using:
        // complex_fft::Normalise(m_workspace);
        // To minimise the amount of looping over the spectrum we will wrap up
        // normalisation and window gain correction into one step.

        // Convert FFT spectrum to power spectrum, converting to real valued vector
        // as this speeds up amount of processing that needs doing below.
        complex_fft::ToPower(m_workspace, realSpectra, fullSpectrum);

        // Compute phases if required.
        if (nullptr != phases)
        {
            phases->resize(realSpectra.size());

            for (size_t i = 0; i < phases->size(); ++i)
            {
                (*phases)[i] = std::arg(m_workspace[i]);
            }
        }

        // Apply window gain correction and normalisation to power spectrum.
        auto normScalar = static_cast<FloatType>(m_workspace.size() * m_workspace.size());
        window_fn::ApplyGainCorrection(realSpectra.begin(),
                                       realSpectra.end(),
                                       realSpectra.begin(),
                                       m_windowFunction.CombinedGain() * normScalar);

        // Perform 3-BIN summing of corrected power spectrum.
        complex_fft::To3BinSum(realSpectra);
    }

    /*!
     * \brief Perform the FFT on the complex signal data samples.
     * \param[in] signal - vector of complex signal data.
     * \param[out] realSpectra - real-valued vector to hold result spectrum.
     * \param[in] fullSpectrum - keep the full spectrum.
     * \param[out] phases - (optional) bin by bin spectrum phases.
     */
    void operator()(complex_vector const& signal, real_vector& realSpectra,
                    bool fullSpectrum = false, real_vector* phases = nullptr)
    {
        DSP_ASSERT_THROW(signal.size() == m_workspace.size(), "signal length is incorrect");

        // Window the data.
        m_windowFunction(signal.begin(), signal.end(), m_workspace.begin());

        // Perform the forward FFT on the windowed data.
        complex_fft::Forward(m_workspace);

        // Typically you would normalise the FFT spectrum at this stage using:
        // complex_fft::Normalise(m_workspace);
        // To minimise the amount of looping over the spectrum we will wrap up
        // normalisation and window gain correction into one step.

        // Convert FFT spectrum to power spectrum, converting to real valued vector
        // as this speeds up amount of processing that needs doing below.
        complex_fft::ToPower(m_workspace, realSpectra, fullSpectrum);

        // Compute phases if required.
        if (nullptr != phases)
        {
            phases->resize(realSpectra.size());

            for (size_t i = 0; i < phases->size(); ++i)
            {
                (*phases)[i] = std::arg(m_workspace[i]);
            }
        }

        // Apply window gain correction and normalisation to power spectrum.
        auto normScalar = static_cast<FloatType>(m_workspace.size() * m_workspace.size());
        window_fn::ApplyGainCorrection(realSpectra.begin(),
                                       realSpectra.end(),
                                       realSpectra.begin(),
                                       m_windowFunction.CombinedGain() * normScalar);

        // Perform 3-BIN summing of corrected power spectrum.
        complex_fft::To3BinSum(realSpectra);
    }

private:
    /*! \brief Window function. */
    window_fn m_windowFunction{};
    /*! \brief Complex data vector workspace. */
    complex_vector m_workspace{};
};

/*! \brief Class to generate real-valued (magnitude) FFT spectrum from real valued data. */
template <typename FloatType> class MagnitudeFft final
{
    /*! \brief Typedef to window function. */
    using window_fn = WindowFunction<FloatType>;
    /*! \brief Typedef to complex floating point value */
    using complex_type = std::complex<FloatType>;
    /*! \brief Typedef to complex floating point vector */
    using complex_vector = std::vector<complex_type>;
    /*! \brief Typedef to real floating point vector */
    using real_vector = std::vector<FloatType>;
    /*! \brief Convenience typedef to ComplexFFT<FloatType>. */
    using complex_fft = ComplexFFT<FloatType>;

public:
    /*! \brief Default constructor. */
    MagnitudeFft() = default;
    /*! \brief Destructor. */
    ~MagnitudeFft() = default;
    /*! \brief Default copy constructor. */
    MagnitudeFft(MagnitudeFft const&) = default;
    /*! \brief Default move constructor. */
    MagnitudeFft(MagnitudeFft&&) = default;
    /*! \brief Default copy assignment operator. */
    MagnitudeFft& operator=(MagnitudeFft const&) = default;
    /*! \brief Default move assignment operator. */
    MagnitudeFft& operator=(MagnitudeFft&&) = default;

    /*!
     * \brief Initialisation constructor.
     * \param[in] generator - the window function generator functor.
     * \param[in] fftSize - the number of samples used by the FFT.
     */
    template <typename Generator>
    MagnitudeFft(Generator const& generator, size_t fftSize)
        : m_workspace(fftSize)
        , m_windowFunction(generator, fftSize + 1, true)
    {
        DSP_ASSERT_THROW(IsPowerOf2(fftSize), "FFT size not a power of 2");
    }

    /*!
     * \brief Initialiser method.
     * \param[in] generator - the window function generator functor.
     * \param[in] fftSize - the number of samples used by the FFT.
     */
    template <typename Generator> void Initialise(Generator const& generator, size_t fftSize)
    {
        *this = std::move(MagnitudeFft(generator, fftSize));
    }

    /*!
     * \brief Perform the FFT on the real-valued signal data samples.
     * \param[in] signalFirst - iterator to first sample of signal data.
     * \param[in] signalLast - iterator to one past the last sample of signal data.
     * \param[out] realSpectra - real-valued vector to hold result spectrum.
     * \param[in] fullSpectrum - keep the full spectrum.
     * \param[out] phases - (optional) bin by bin spectrum phases.
     */
    template <typename Iter>
    void operator()(Iter signalFirst, Iter signalLast, real_vector& realSpectra,
                    bool fullSpectrum = false, real_vector* phases = nullptr)
    {
        auto signalLength = std::distance(signalFirst, signalLast);
        DSP_ASSERT_THROW(signalLength == static_cast<decltype(signalLength)>(m_workspace.size()),
                         "signal length is incorrect");

        // Window the data.
        m_windowFunction(signalFirst, signalLast, realSpectra.begin());

        // Copy data into workspace.
        auto convertToComplex = [](auto x) { return complex_type{x, 0}; };
        std::transform(
            realSpectra.begin(), realSpectra.end(), m_workspace.begin(), convertToComplex);

        // Perform the forward FFT on the windowed data.
        complex_fft::Forward(m_workspace);

        // Typically you would normalise the FFT spectrum at this stage using:
        // complex_fft::Normalise(m_workspace);
        // To minimise the amount of looping over the spectrum we will wrap up
        // normalisation and window gain correction into one step.

        // Convert FFT spectrum to power spectrum, converting to real valued vector
        // as this speeds up amount of processing that needs doing below.
        complex_fft::ToMagnitude(m_workspace, realSpectra, fullSpectrum);

        // Compute phases if required.
        if (nullptr != phases)
        {
            phases->resize(realSpectra.size());

            for (size_t i = 0; i < phases->size(); ++i)
            {
                (*phases)[i] = std::arg(m_workspace[i]);
            }
        }

        // Apply window gain correction and normalisation to power spectrum.
        auto normScalar = static_cast<FloatType>(m_workspace.size());
        window_fn::ApplyGainCorrection(realSpectra.begin(),
                                       realSpectra.end(),
                                       realSpectra.begin(),
                                       m_windowFunction.CoherentGain() * normScalar);
    }

    /*!
     * \brief Perform the FFT on the complex signal data samples.
     * \param[in] signal - vector of complex signal data.
     * \param[out] realSpectra - real-valued vector to hold result spectrum.
     * \param[in] fullSpectrum - keep the full spectrum.
     * \param[out] phases - (optional) bin by bin spectrum phases.
     */
    void operator()(complex_vector const& signal, real_vector& realSpectra,
                    bool fullSpectrum = false, real_vector* phases = nullptr)
    {
        DSP_ASSERT_THROW(signal.size() == m_workspace.size(), "signal length is incorrect");

        // Window the data.
        m_windowFunction(signal.begin(), signal.end(), m_workspace.begin());

        // Perform the forward FFT on the windowed data.
        complex_fft::Forward(m_workspace);

        // Typically you would normalise the FFT spectrum at this stage using:
        // complex_fft::Normalise(m_workspace);
        // To minimise the amount of looping over the spectrum we will wrap up
        // normalisation and window gain correction into one step.

        // Convert FFT spectrum to power spectrum, converting to real valued vector
        // as this speeds up amount of processing that needs doing below.
        complex_fft::ToMagnitude(m_workspace, realSpectra, fullSpectrum);

        // Compute phases if required.
        if (nullptr != phases)
        {
            phases->resize(realSpectra.size());

            for (size_t i = 0; i < phases->size(); ++i)
            {
                (*phases)[i] = std::arg(m_workspace[i]);
            }
        }

        // Apply window gain correction and normalisation to power spectrum.
        auto normScalar = static_cast<FloatType>(m_workspace.size());
        window_fn::ApplyGainCorrection(realSpectra.begin(),
                                       realSpectra.end(),
                                       realSpectra.begin(),
                                       m_windowFunction.CoherentGain() * normScalar);
    }

private:
    /*! \brief Window function. */
    window_fn m_windowFunction{};
    /*! \brief Complex data vector workspace. */
    complex_vector m_workspace{};
};

/*! \brief Class to perform convolution of 2 ranges using FFTs. */
template <typename FloatType> class FftConvolve final
{
    /*! \brief Typedef to complex floating point value */
    using complex_type = std::complex<FloatType>;
    /*! \brief Typedef to complex floating point vector */
    using complex_vector = std::vector<complex_type>;
    /*! \brief Typedef to real floating point vector */
    using real_vector = std::vector<FloatType>;
    /*! \brief Convenience typedef to ComplexFFT<FloatType>. */
    using complex_fft = ComplexFFT<FloatType>;

public:
    /*! \brief Default constructor. */
    FftConvolve() = default;
    /*! \brief Destructor. */
    ~FftConvolve() = default;
    /*! \brief Default copy constructor. */
    FftConvolve(FftConvolve const&) = default;
    /*! \brief Default move constructor. */
    FftConvolve(FftConvolve&&) = default;
    /*! \brief Default copy assignment operator. */
    FftConvolve& operator=(FftConvolve const&) = default;
    /*! \brief Default move assignment operator. */
    FftConvolve& operator=(FftConvolve&&) = default;

    /*!
     * \brief Initialisation constructor.
     * \param[in] signalLength - data length for range of signal's samples.
     * \param[in] kernelLength - data length for filter kernel.
     */
    FftConvolve(size_t signalLength, size_t kernelLength)
    {
        DSP_ASSERT_THROW(signalLength > 0, "signalLength <= 0");
        DSP_ASSERT_THROW(kernelLength > 0, "kernelLength <= 0");
        m_discreteConvolutionLength = signalLength + kernelLength - 1;
        auto powerOf2 =
            static_cast<size_t>(std::floor(std::log2(m_discreteConvolutionLength)) + 0.5);
        size_t workspaceLength = static_cast<size_t>(1) << (powerOf2 + 1);
        m_workspace1.resize(workspaceLength, 0);
        m_workspace2.resize(workspaceLength, 0);
    }

    /*!
     * \brief Initialiser method.
     * \param[in] signalLength - data length for range of signal's samples.
     * \param[in] kernelLength - data length for filter kernel.
     */
    void Initialise(size_t signalLength, size_t kernelLength)
    {
        *this = std::move(FftConvolve(signalLength, kernelLength));
    }

    /*!
     * \brief Perform convolution on 2 data ranges.
     * \param[in] signalFirst - Start of the signal sample range, length N.
     * \param[in] signalLast - End of the signal sample range, length N.
     * \param[in] kernelFirst - Start of the filter kernel range, length M.
     * \param[in] kernelLast - End of the filter kernel range, length M.
     * \param[out] result - Start of the result range must be length N + M - 1.
     *
     * Performs convolution using FFTs.
     *
     * 1. Size workspaces to power of 2 greater than N + M - 1.
     * 2. Fill work space with data ranges adding appropriate 0 padding to the end of the
     * workspaces.
     * 3. Perform forward FFTs on the workspaces.
     * 4. Piecewise multiplication of 2 spectras complex values.
     * 5. Perform inverse FFT on multiplied result.
     * 6. Copy real values out of the complex result.
     */
    template <typename InIter1, typename InIter2, typename OutIter>
    void operator()(InIter1 signalFirst, InIter1 signalLast, InIter2 kernelFirst,
                    InIter2 kernelLast, OutIter result)
    {
        // TODO: This all works but rewrite to use overlap and add method where
        // if the input signal lengths are particularly big then we don't want the
        // FFT convolution object to create large workspaces (and large FFTs) so
        // using overlap and add breaks the problem down into smaller chunks of
        // source data requiring smaller sized FFTs and workspace and stitches
        // the results back together. So we may want an FFT size limit under
        // which we use a single FFT for the convolution but if data length requires
        // an FFT above the size limit we break down into a series of smaller FFTs
        // based on picking an FFT size bigger than the shortest data range but with
        // enough spare samples to sensibly convolve with the other data range.
        // We probably in this case want to say range 1 should always be the
        // source data and range 2 is the kernel to be convolved with, e.g. the filter
        // kernel.
        //
        // The current method below is at least 3 times faster than using discrete
        // convolution.
        auto signalLength = std::distance(signalFirst, signalLast);
        DSP_ASSERT_THROW(signalLength <= static_cast<decltype(signalLength)>(m_workspace1.size()),
                         "signal length is incorrect");
        auto kernelLength = std::distance(kernelFirst, kernelLast);
        DSP_ASSERT_THROW(kernelLength <= static_cast<decltype(kernelLength)>(m_workspace2.size()),
                         "kernel length is incorrect");

        std::fill(
            std::next(m_workspace1.begin(), signalLength), m_workspace1.end(), complex_type{0, 0});
        std::fill(
            std::next(m_workspace2.begin(), kernelLength), m_workspace2.end(), complex_type{0, 0});

        auto convertToComplex = [](auto x) { return complex_type{x, 0}; };

        std::transform(signalFirst, signalLast, m_workspace1.begin(), convertToComplex);
        std::transform(kernelFirst, kernelLast, m_workspace2.begin(), convertToComplex);

        complex_fft::Forward(m_workspace1);
        complex_fft::Forward(m_workspace2);

        auto crossMultiply = [](auto const& x, auto const& y) {
            auto re = x.real() * y.real() - x.imag() * y.imag();
            auto im = x.imag() * y.real() + x.real() * y.imag();
            return complex_type{re, im};
        };

        std::transform(m_workspace1.begin(),
                       m_workspace1.end(),
                       m_workspace2.begin(),
                       m_workspace1.begin(),
                       crossMultiply);

        complex_fft::Inverse(m_workspace1);

        auto outFirst = m_workspace1.begin();
        auto outLast  = std::next(outFirst, static_cast<int>(m_discreteConvolutionLength));

        auto convertToReal = [](auto const& z) { return z.real(); };
        std::transform(outFirst, outLast, result, convertToReal);
    }

private:
    /*! \brief Length we'd expected due to discrete convolution. */
    size_t m_discreteConvolutionLength{0};
    /*! \brief Complex data vector workspace 1. */
    complex_vector m_workspace1{};
    /*! \brief Complex data vector workspace 2. */
    complex_vector m_workspace2{};
};

/*! \brief Convenience typedef to ComplexFFT<float>. */
using complex_fft_f = ComplexFFT<float>;
/*! \brief Convenience typedef to ComplexFFT<double>. */
using complex_fft_d = ComplexFFT<double>;
/*! \brief Convenience typedef to ComplexFFT<long double>. */
using complex_fft_ld = ComplexFFT<long double>;
/*! \brief Convenience typedef to ThreeBinSumFft<float>. */
using three_bin_sum_fft_f = ThreeBinSumFft<float>;
/*! \brief Convenience typedef to ThreeBinSumFft<double>. */
using three_bin_sum_fft_d = ThreeBinSumFft<double>;
/*! \brief Convenience typedef to ThreeBinSumFft<long double>. */
using three_bin_sum_fft_ld = ThreeBinSumFft<long double>;
/*! \brief Convenience typedef to MagnitudeFft<float>. */
using magnitude_fft_f = MagnitudeFft<float>;
/*! \brief Convenience typedef to MagnitudeFft<double>. */
using magnitude_fft_d = MagnitudeFft<double>;
/*! \brief Convenience typedef to MagnitudeFft<long double>. */
using magnitude_fft_ld = MagnitudeFft<long double>;
/*! \brief Convenience typedef to MagnitudeFft<float>. */
using fft_convolve_f = FftConvolve<float>;
/*! \brief Convenience typedef to MagnitudeFft<double>. */
using fft_convolve_d = FftConvolve<double>;
/*! \brief Convenience typedef to MagnitudeFft<long double>. */
using fft_convolve_ld = FftConvolve<long double>;

} // namespace dsp

#endif // DSP_FFT_HPP
