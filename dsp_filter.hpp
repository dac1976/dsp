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
 * \file dsp_filter.hpp
 * \brief File containing FIR filter functions.
 *
 * These filter functins were inspired by the examples on http://www.iowahills.com/.
 */

#ifndef DSP_FILTER_HPP
#define DSP_FILTER_HPP

#include "dsp_fft.hpp"

/*! \brief dsp namespace */
namespace dsp
{

/*!
 * \brief FIR Low Pass Filter
 * \param[in] numTaps - Number of taps, i.e. number of filter coefficients.
 * \param[in] cutoffFreqHz - Cut-off frequency in Hz above which we'll filter out the signal.
 * \param[in] samplingFreqHz - Frequency the signal data has been sampled at in Hz.
 * \param[in] windowGenerator - A window generator functor to apply to the filter coefficients.
 * \return A vector filter coefficients.
 *
 * It is recommended to use a Kaiser window, e.g. KaiserGenerator. You can specify an odd or
 * even number of taps  but it is often better to use an odd number of taps to have a single
 * centre point around which the filter is symmetrical.
 */
template <typename FloatType, typename WinGenType>
std::vector<FloatType> FirLowPassFilter(size_t numTaps, FloatType cutoffFreqHz,
                                        FloatType samplingFreqHz, WinGenType const& windowGenerator)
{
    DSP_ASSERT_THROW(numTaps > 2, "numTaps too small");
    DSP_ASSERT_THROW(cutoffFreqHz > 0, "cutoffFreqHz <= 0");
    DSP_ASSERT_THROW(samplingFreqHz > 0, "samplingFreqHz  <= 0");
    auto availableBandwidthHz = samplingFreqHz / 2;
    DSP_ASSERT_THROW(cutoffFreqHz <= availableBandwidthHz, "cutoffFreqHz too high");
    static const auto      PI = Pi<FloatType>();
    std::vector<FloatType> filterCoeffs(numTaps);
    auto                   normalizedCutoffFreq   = cutoffFreqHz / availableBandwidthHz;
    auto                   numTapsMinusOneOverTwo = static_cast<FloatType>(numTaps - 1) / 2;
    size_t                 i                      = 0;

    for (auto& filterCoeff : filterCoeffs)
    {
        auto arg    = static_cast<FloatType>(i) - numTapsMinusOneOverTwo;
        filterCoeff = normalizedCutoffFreq * Sinc(normalizedCutoffFreq * arg * PI);
        ++i;
    }

    WindowFunction<FloatType> window(windowGenerator, numTaps, false);
    window(filterCoeffs.begin(), filterCoeffs.end(), filterCoeffs.begin());
    return filterCoeffs;
}

/*!
 * \brief FIR High Pass Filter
 * \param[in] numTaps - Number of taps, i.e. number of filter coefficients.
 * \param[in] cutoffFreqHz - Cut-off frequency in Hz below which we'll filter out the signal.
 * \param[in] samplingFreqHz - Frequency the signal data has been sampled at in Hz.
 * \param[in] windowGenerator - A window generator functor to apply to the filter coefficients.
 * \return A vector filter coefficients.
 *
 * It is recommended to use a Kaiser window, e.g. KaiserGenerator. The number of taps should be
 * odd for a high pass filter else the resultant filter will have undesirable zeroes and the
 * filtered signal will be attenuated.
 */
template <typename FloatType, typename WinGenType>
std::vector<FloatType> FirHighPassFilter(size_t numTaps, FloatType cutoffFreqHz,
                                         FloatType         samplingFreqHz,
                                         WinGenType const& windowGenerator)
{
    DSP_ASSERT_THROW(numTaps > 2, "numTaps too small");
    DSP_ASSERT_THROW(cutoffFreqHz > 0, "cutoffFreqHz <= 0");
    DSP_ASSERT_THROW(samplingFreqHz > 0, "samplingFreqHz  <= 0");
    auto availableBandwidthHz = samplingFreqHz / 2;
    DSP_ASSERT_THROW(cutoffFreqHz <= availableBandwidthHz, "cutoffFreqHz too high");
    bool isOdd = numTaps % 2 == 1;
    DSP_ASSERT_THROW(isOdd, "numTaps should be an odd number for high pass filter");
    static const auto      PI = Pi<FloatType>();
    std::vector<FloatType> filterCoeffs(numTaps);
    auto                   normalizedCutoffFreq   = cutoffFreqHz / availableBandwidthHz;
    auto                   numTapsMinusOneOverTwo = static_cast<FloatType>(numTaps - 1) / 2;
    size_t                 i                      = 0;

    for (auto& filterCoeff : filterCoeffs)
    {
        auto arg = static_cast<FloatType>(i) - numTapsMinusOneOverTwo;
        filterCoeff =
            Sinc(arg * PI) - (normalizedCutoffFreq * Sinc(normalizedCutoffFreq * arg * PI));
        ++i;
    }

    WindowFunction<FloatType> window(windowGenerator, numTaps, false);
    window(filterCoeffs.begin(), filterCoeffs.end(), filterCoeffs.begin());
    return filterCoeffs;
}

/*!
 * \brief FIR Band Pass Filter
 * \param[in] numTaps - Number of taps, i.e. number of filter coefficients.
 * \param[in] centreFreqHz - Centre frequency in Hz around which we'll retain the signal content.
 * \param[in] bandwidthHz - Bandwidth in Hz used with the centre frequency to define the pass band.
 * \param[in] samplingFreqHz - Frequency the signal data has been sampled at in Hz.
 * \param[in] windowGenerator - A window generator functor to apply to the filter coefficients.
 * \return A vector filter coefficients.
 *
 * It is recommended to use a Kaiser window, e.g. KaiserGenerator. You can specify an odd or
 * even number of taps  but it is often better to use an odd number of taps to have a single
 * centre point around which the filter is symmetrical.
 */
template <typename FloatType, typename WinGenType>
std::vector<FloatType> FirBandPassFilter(size_t numTaps, FloatType centreFreqHz,
                                         FloatType bandwidthHz, FloatType samplingFreqHz,
                                         WinGenType const& windowGenerator)
{
    DSP_ASSERT_THROW(numTaps > 2, "numTaps too small");
    DSP_ASSERT_THROW(centreFreqHz > 0, "cutoffFreqHz <= 0");
    DSP_ASSERT_THROW(samplingFreqHz > 0, "samplingFreqHz  <= 0");
    auto availableBandwidthHz = samplingFreqHz / 2;
    DSP_ASSERT_THROW(centreFreqHz <= availableBandwidthHz, "cutoffFreqHz too high");
    DSP_ASSERT_THROW(bandwidthHz > 0, "bandwidthHz  <= 0");
    DSP_ASSERT_THROW(bandwidthHz <= availableBandwidthHz, "bandwidthHz too high");
    static const auto      PI = Pi<FloatType>();
    std::vector<FloatType> filterCoeffs(numTaps);
    auto                   normalizedCentreFreq = centreFreqHz / availableBandwidthHz;
    auto                   normalizedBandwidth  = bandwidthHz / availableBandwidthHz;
    auto   normalizedLowCutoffFreq              = normalizedCentreFreq - (normalizedBandwidth / 2);
    auto   normalizeHighCutoffFreq              = normalizedCentreFreq + (normalizedBandwidth / 2);
    auto   numTapsMinusOneOverTwo               = static_cast<FloatType>(numTaps - 1) / 2;
    size_t i                                    = 0;

    for (auto& filterCoeff : filterCoeffs)
    {
        auto arg    = static_cast<FloatType>(i) - numTapsMinusOneOverTwo;
        filterCoeff = std::abs(arg) < FloatType(1.e-3) ? FloatType(0.)
                                                       : (cos(normalizedLowCutoffFreq * arg * PI) -
                                                          cos(normalizeHighCutoffFreq * arg * PI)) /
                                                             PI / arg;
        ++i;
    }

    WindowFunction<FloatType> window(windowGenerator, numTaps, false);
    window(filterCoeffs.begin(), filterCoeffs.end(), filterCoeffs.begin());
    return filterCoeffs;
}

/*!
 * \brief FIR Notch Filter
 * \param[in] numTaps - Number of taps, i.e. number of filter coefficients.
 * \param[in] centreFreqHz - Centre frequency in Hz around which we'll retain the signal content.
 * \param[in] bandwidthHz - Bandwidth in Hz used with the centre frequency to define the rejection
 * band.
 * \param[in] samplingFreqHz - Frequency the signal data has been sampled at in Hz.
 * \param[in] windowGenerator - A window generator functor to apply to the filter coefficients.
 * \return A vector filter coefficients.
 *
 * It is recommended to use a Kaiser window. The number of taps should be odd for a high pass
 * filter else the resultant filter will have undesirable zeroes and the filtered signal will
 * be attenuated.
 */
template <typename FloatType, typename WinGenType>
std::vector<FloatType> FirNotchFilter(size_t numTaps, FloatType centreFreqHz, FloatType bandwidthHz,
                                      FloatType samplingFreqHz, WinGenType const& windowGenerator)
{
    DSP_ASSERT_THROW(numTaps > 2, "numTaps too small");
    DSP_ASSERT_THROW(centreFreqHz > 0, "cutoffFreqHz <= 0");
    DSP_ASSERT_THROW(samplingFreqHz > 0, "samplingFreqHz  <= 0");
    auto availableBandwidthHz = samplingFreqHz / 2;
    DSP_ASSERT_THROW(centreFreqHz <= availableBandwidthHz, "cutoffFreqHz too high");
    DSP_ASSERT_THROW(bandwidthHz > 0, "bandwidthHz  <= 0");
    DSP_ASSERT_THROW(bandwidthHz <= availableBandwidthHz, "bandwidthHz too high");
    static const auto      PI = Pi<FloatType>();
    std::vector<FloatType> filterCoeffs(numTaps);
    auto                   normalizedCentreFreq = centreFreqHz / availableBandwidthHz;
    auto                   normalizedBandwidth  = bandwidthHz / availableBandwidthHz;
    auto   normalizedLowCutoffFreq              = normalizedCentreFreq - (normalizedBandwidth / 2);
    auto   normalizeHighCutoffFreq              = normalizedCentreFreq + (normalizedBandwidth / 2);
    auto   numTapsMinusOneOverTwo               = static_cast<FloatType>(numTaps - 1) / 2;
    size_t i                                    = 0;

    for (auto& filterCoeff : filterCoeffs)
    {
        auto arg    = static_cast<FloatType>(i) - numTapsMinusOneOverTwo;
        filterCoeff = Sinc(arg * PI) -
                      (normalizeHighCutoffFreq * Sinc(normalizeHighCutoffFreq * arg * PI)) -
                      (normalizedLowCutoffFreq * Sinc(normalizedLowCutoffFreq * arg * PI));
        ++i;
    }

    WindowFunction<FloatType> window(windowGenerator, numTaps, false);
    window(filterCoeffs.begin(), filterCoeffs.end(), filterCoeffs.begin());
    return filterCoeffs;
}

/*! \brief Filter holder class to manage applying filter coefficients to a signal. */
template <typename FloatType> class FilterHolder final
{
    /*! \brief Typedef to FftConvolve functor. */
    using fft_convolve_t = FftConvolve<FloatType>;

public:
    /*! \brief Default constructor. */
    FilterHolder() = default;
    /*! \brief Destructor. */
    ~FilterHolder() = default;
    /*! \brief Default copy constructor. */
    FilterHolder(FilterHolder const&) = default;
    /*! \brief Default move constructor. */
    FilterHolder(FilterHolder&&) = default;
    /*! \brief Default copy assignment operator. */
    FilterHolder& operator=(FilterHolder const&) = default;
    /*! \brief Default move assignment operator. */
    FilterHolder& operator=(FilterHolder&&) = default;

    /*!
     * \brief Initialisation constructor.
     * \param[in] signalLength - th length of the signal samples to be filtered.
     * \param[in] filterCoeffs - the filter coefficients to hold and apply.
     * \param[in] useFastConvolution - choose whether to use fast FFT based convolution or not.
     */
    FilterHolder(size_t signalLength, std::vector<FloatType> const& filterCoeffs,
                 bool useFastConvolution)
        : m_signalLength(signalLength)
        , m_filterCoeffs(filterCoeffs)
        , m_useFastConvolution(useFastConvolution)
        , m_filteredSignal(signalLength + m_filterCoeffs.size() - 1)
    {
        DSP_ASSERT_THROW(m_signalLength > 2, "signalLength is too small");
        DSP_ASSERT_THROW(!m_filterCoeffs.empty(), "filterCoeffs is empty");

        if (m_useFastConvolution)
        {
            m_fftConvolve = fft_convolve_t(signalLength, filterCoeffs.size());
        }
    }

    /*!
     * \brief Initialiser method.
     * \param[in] signalLength - th length of the signal samples to be filtered.
     * \param[in] filterCoeffs - the filter coefficients to hold and apply.
     * \param[in] useFastConvolution - choose whether to use fast FFT based convolution or not.
     */
    void Initialise(size_t signalLength, std::vector<FloatType> const& filterCoeffs,
                    bool useFastConvolution)
    {
        *this = std::move(FilterHolder(signalLength, filterCoeffs, useFastConvolution));
    }

    /*!
     * \brief Function operator used to apply filter to signal.
     * \param[in] signalFirst - Iterator to first signal sample.
     * \param[in] signalLast - Iterator to one past last actual signal sample.
     * \param[in] removeDelay - When true the delay introduced while applying the filter will ne
     * removed.
     * \param[out] resultFirst - Iterator to start of correctly sized result container.
     *
     * If the user opts to not remove the delay then the returned vector will contain N + M - 1
     * samples, where N is the number of signal samples and M is the number of filter coefficients.
     * If the user opts to remove the delay then the returned vector will contain N samples.
     */
    template <typename IterIn, typename OutIter>
    void operator()(IterIn signalFirst, IterIn signalLast, OutIter resultFirst, bool removeDelay)
    {
        auto signalLen = std::distance(signalFirst, signalLast);
        DSP_ASSERT_THROW(signalLen == static_cast<decltype(signalLen)>(m_signalLength),
                         "signal sample range incorrect");

        if (m_useFastConvolution)
        {
            m_fftConvolve(signalFirst,
                          signalLast,
                          m_filterCoeffs.begin(),
                          m_filterCoeffs.end(),
                          m_filteredSignal.begin());
        }
        else
        {

            Convolve(signalFirst,
                     signalLast,
                     m_filterCoeffs.begin(),
                     m_filterCoeffs.end(),
                     m_filteredSignal.begin());
        }

        if (removeDelay)
        {
            auto offset = m_filteredSignal.size() - m_signalLength;
            offset >>= 1;
            auto first = std::next(m_filteredSignal.begin(), static_cast<int>(offset));
            auto last  = std::next(first, signalLen);

            std::copy(first, last, resultFirst);
        }
    }

private:
    /*! \brief Copy of the filter coefficients. */
    size_t m_signalLength{0};
    /*! \brief Copy of the filter coefficients. */
    std::vector<FloatType> m_filterCoeffs{};
    /*! \brief Use fast FFT convolution flag. */
    bool m_useFastConvolution{true};
    /*! \brief FFT convolution functor. */
    fft_convolve_t m_fftConvolve{};
    /*! \brief Workspace vector holding the filtered signal. */
    std::vector<FloatType> m_filteredSignal{};
};

/*! \brief Convenience typedef to FilterHolder<float>. */
using filter_hldr_f = FilterHolder<float>;
/*! \brief Convenience typedef to FilterHolder<double>. */
using filter_hldr_d = FilterHolder<double>;
/*! \brief Convenience typedef to FilterHolder<long double>. */
using filter_holder_ld = FilterHolder<long double>;

} // namespace dsp

#endif // DSP_FILTER_HPP
