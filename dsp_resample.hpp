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
 * \file dsp_resample.hpp
 * \brief File containing generic signal resampling functions.
 */

#ifndef DSP_RESAMPLE_HPP
#define DSP_RESAMPLE_HPP

#include <utility>
#include <limits>
#include "dsp_filter.hpp"
#include "dsp_signals.hpp"

/*! \brief dsp namespace */
namespace dsp
{

/*!
 * \brief Round a floating type's value down to the integer value below, converting to an integer
 *        type.
 * \param[in] value - Floating point value to round.
 * \return Integer value.
 */
template <typename FloatT, typename IntegerT> IntegerT FloatToInt(FloatT value)
{
    static_assert(std::is_floating_point<FloatT>::value, "invalid floating point type");
    static_assert(std::is_integral<IntegerT>::value, "invalid integral point type");

    return (value < FloatT(0)) ? IntegerT(std::ceil(value) - 0.5)
                               : IntegerT(std::floor(value) + 0.5);
}

/*!
 * \brief Resample the input range into the output range using linear interpolation.
 * \param[in] first - First item of source range.
 * \param[in] last - Last item of source range, same role as end iterator on an STL container.
 * \param[out] targetFirst - First item of target range.
 * \param[out] targetLast - Last item of target range, same role as end iterator on an STL
 * container.
 *
 * This function works out the resampling factor from the ratio of the size of the
 * source range to the target range.
 */
template <typename InIter, typename OutIter>
void ResampleRange(InIter first, InIter last, OutIter targetFirst, OutIter targetLast)
{
    // Do we have a valid source range?
    auto tmpSourceSize = std::distance(first, last);

    DSP_ASSERT_THROW(tmpSourceSize > 0, "std::distance(first, last) <= 0");

    auto tmpTargetSize = std::distance(targetFirst, targetLast);

    DSP_ASSERT_THROW(tmpTargetSize > 0, "std::distance(targetFirst, targetLast) <= 0");

    // Do we need to do any downsampling? If not copy
    // source range into target vector and return early.
    const auto sourceSize = static_cast<size_t>(tmpSourceSize);
    const auto targetSize = static_cast<size_t>(tmpTargetSize);

    if (sourceSize == targetSize)
    {
        std::copy(first, last, targetFirst);
        return;
    }

    // Work out the exact real-valued sample stride.
    // This means we can find the exact sample position
    // in the original data where our resampled value
    // should be taken from for the downsampled data.
    // The exact sample position is likely to lie between
    // 2 given samples of the source data.
    const auto sampleStride =
        static_cast<double>(sourceSize - 1) / static_cast<double>(targetSize - 1);

    double exactSamplePos = 0.0;
    auto   finalItem      = std::next(first, sourceSize - 1);

    // Tracking iterator and sample pos.
    auto   itrBefore    = first;
    size_t sampleBefore = 0;
    size_t pos          = 0;
    using out_type_t    = typename std::iterator_traits<OutIter>::value_type;

    for (auto outIter = targetFirst; outIter != targetLast; std::advance(outIter, 1), ++pos)
    {
        out_type_t interpolatedSample;

        // Keep endpoints as they are.
        if (0 == pos)
        {
            interpolatedSample = *first;
        }
        else if (pos == targetSize - 1)
        {
            interpolatedSample = *finalItem;
        }
        // Downsample into smaller buffer by using
        // linear interpolation to maintain more
        // accurate amplitude values.
        else
        {
            // Remember previous sample before last exact pos.
            const size_t prevSample = sampleBefore;
            // Find the sample in the source data just
            // before our exact sample position.
            sampleBefore = FloatToInt<double, size_t>(exactSamplePos);
            // Work our the distance of our exact sample position
            // from the sample before as a ratio.
            const double ratio = exactSamplePos - double(sampleBefore);
            std::advance(itrBefore, sampleBefore - prevSample);
            auto itrAfter = std::next(itrBefore);
            // Work out an interpolated correction factor to
            // add to the sample value of the sample before
            // our exact sample position.
            const double correction = (*itrAfter - *itrBefore) * ratio;
            interpolatedSample      = *itrBefore + out_type_t(correction);
        }

        *outIter = interpolatedSample;
        exactSamplePos += sampleStride;
    }
}

/*!
 * \brief Compute closest resample up and down factors given a real valued resample factor.
 * \param[in] requiredResampleFactor - The resample factor we ideally want to use, should be != 1.
 * \param[in] maxNumerator - The maximum numerator to allow in the final result.
 * \param[in] maxDenominator - The maximum denominator to allow in the final result.
 * \return A pair containing the upsample factor and the downsample factor.
 *
 * This function uses the concept of mediants to interval bisect between an upper and lower
 * bound either side of the requiredResampleFactor until the termination conditions are
 * met and a suitable pair of resample factors is found.
 */
inline std::pair<size_t, size_t> ComputeResampleFactors(double requiredResampleFactor,
                                                        size_t maxNumerator   = 128,
                                                        size_t maxDenominator = 128)
{
    DSP_ASSERT_THROW(requiredResampleFactor > 0, "requiredResampleFactor <= 0");
    std::pair<size_t, size_t> factors(0, 0);
    auto                      n_a   = static_cast<size_t>(std::floor(requiredResampleFactor) + 0.5);
    size_t                    d_a   = 1;
    auto                      n_b   = static_cast<size_t>(std::ceil(requiredResampleFactor) + 0.5);
    size_t                    d_b   = 1;
    double                    error = std::numeric_limits<double>::max();

    while (true)
    {
        auto n_m = n_a + n_b;
        auto d_m = d_a + d_b;

        // Make sure we reduce the numerator and denominator as much as possible.
        // Although in practice we may never need this here as n_m and d_m are
        // likely already in their lowest form.
        auto g = Gcd(n_m, d_m);

        if (g > 1)
        {
            n_m /= g;
            d_m /= g;
        }

        if ((n_m > maxNumerator) || (d_m > maxDenominator))
        {
            break;
        }

        double m       = static_cast<double>(n_m) / static_cast<double>(d_m);
        auto   absDiff = std::abs(m - requiredResampleFactor);

        if (absDiff < error)
        {
            error   = absDiff;
            factors = {n_m, d_m};
        }

        if (m <= requiredResampleFactor)
        {
            n_a = n_m;
            d_a = d_m;
        }
        else
        {
            n_b = n_m;
            d_b = d_m;
        }
    }

    return factors;
}

/*! \brief Resampling class uses FIR filter and Kaiser window.
 *
 * This class is designed to be created upfront and then reused repeatedly to
 * perform ongoing resampling of a given signal and so it creates its internal
 * objects, such as filters, window function and workspace buffers upfront.
 *
 * Workspace buffers are allocated by this class the first time it resamples signal
 * samples and will only resize the workspace if required on subsequent calls to
 * resample more signal samples.
 *
 * Therefore creation and first call of this class to perform resampling is slower
 * than subsequent calls due to setting everything up.
 *
 * The general algorithm is as follows:
 * 1. Fill upsample buffer with signal samples spaced with zero padding samples.
 * 2. Correct for signal attenuation caused by zero padding upsampling.
 * 3. Low pass filter the upsample buffer.
 * 4. Downsample by skipping not required samples.
 * 5. Return the resampled signal.
 */
template <typename FloatType> class Resample final
{
    /*! \brief Typedef to FilterHolder. */
    using filter_hldr_t = FilterHolder<FloatType>;

public:
    /*! \brief Default constructor. */
    Resample() = default;
    /*! \brief Destructor. */
    ~Resample() = default;
    /*! \brief Default copy constructor. */
    Resample(Resample const&) = default;
    /*! \brief Default move constructor. */
    Resample(Resample&&) = default;
    /*! \brief Default copy assignment operator. */
    Resample& operator=(Resample const&) = default;
    /*! \brief Default move assignment operator. */
    Resample& operator=(Resample&&) = default;

    /*!
     * \brief Initialisation constructor.
     * \param[in] signalLength - The number of signal samples to be resampled each time.
     * \param[in] upsampleFactor - Factor by which we need to upsample by >= 1.
     * \param[in] downsampleFactor - Factor by which we need to downsample by >= 1.
     * \param[in] samplingFreqHz - Sampling frequency of signal in Hz which we wish to resample.
     * \param[in] maxCutoffFreqHz - The low pass filter max cutoff freq, typically new sample rate /
     * 2 when start freq is 0Hz.
     * \param[in] numFilterTaps - The number of filter coefficients we
     * require.
     * \param[in] kaiserWindowBeta - Beta parameter controlling side lobes of Kaiser
     * window.
     * \param[in] useFastConvolution - choose whether to use fast FFT based convolution or
     * not.
     */
    Resample(size_t signalLength, size_t upsampleFactor, size_t downsampleFactor,
             FloatType samplingFreqHz, FloatType maxCutoffFreqHz, size_t numFilterTaps,
             double kaiserWindowBeta, bool useFastConvolution)
        : m_signalLength(signalLength)
        , m_upsampleFactor(upsampleFactor)
        , m_downsampleFactor(downsampleFactor)
    {
        DSP_ASSERT_THROW(m_signalLength > 0, "too few signal samples");
        DSP_ASSERT_THROW(m_upsampleFactor > 0, "invalid upsample factor");
        DSP_ASSERT_THROW(m_downsampleFactor > 0, "invalid downsample factor");

        ComputeResampledSize();
        ResizeWorkspace();

        auto upsampleLength  = m_upsampleFactor * m_signalLength;
        auto upsampledFreqHz = samplingFreqHz * static_cast<FloatType>(m_upsampleFactor);
        auto resampledFreqHz = samplingFreqHz * static_cast<FloatType>(upsampleFactor) /
                               static_cast<FloatType>(downsampleFactor);
        auto cutoffFreqHz = std::min(samplingFreqHz, resampledFreqHz) / 2;

        if (upsampleFactor > downsampleFactor)
        {
            maxCutoffFreqHz = std::min(cutoffFreqHz, maxCutoffFreqHz);
        }
        else
        {
            maxCutoffFreqHz = std::max(cutoffFreqHz, maxCutoffFreqHz);
        }

        m_filterHolder = filter_hldr_t(
            upsampleLength,
            FirLowPassFilter(
                numFilterTaps, maxCutoffFreqHz, upsampledFreqHz, KaiserGenerator(kaiserWindowBeta)),
            useFastConvolution);
    }

    /*!
     * \brief Initialisation method.
     * \param[in] signalLength - The number of signal samples to be resampled each time.
     * \param[in] upsampleFactor - Factor by which we need to upsample by >= 1.
     * \param[in] downsampleFactor - Factor by which we need to downsample by >= 1.
     * \param[in] samplingFreqHz - Sampling frequency of signal in Hz which we wish to resample.
     * \param[in] maxCutoffFreqHz - The low pass filter cutoff freq, typically new sample rate / 2
     * when start freq is 0Hz.
     * \param[in] numFilterTaps - The number of filter coefficients we
     * require. \param[in] kaiserWindowBeta - Beta parameter controlling side lobes of Kaiser
     * window. \param[in] useFastConvolution - choose whether to use fast FFT based convolution or
     * not.
     */
    void Initialise(size_t signalLength, size_t upsampleFactor, size_t downsampleFactor,
                    FloatType samplingFreqHz, FloatType maxCutoffFreqHz, size_t numFilterTaps,
                    double kaiserWindowBeta, bool useFastConvolution)
    {
        *this = std::move(Resample(signalLength,
                                   upsampleFactor,
                                   downsampleFactor,
                                   samplingFreqHz,
                                   maxCutoffFreqHz,
                                   numFilterTaps,
                                   kaiserWindowBeta,
                                   useFastConvolution));
    }

    /*!
     * \brief Get the original data size.
     * \return The size of the signal data.
     */
    size_t DataSize() const
    {
        return m_signalLength;
    }

    /*!
     * \brief Get the resampled data size.
     * \return The size required for the resampled container.
     */
    size_t ResampledSize() const
    {
        return m_resampledLength;
    }

    /*!
     * \brief Function operator to perform resampling.
     * \param[in] signalFirst - Iterator to first signal sample.
     * \param[in] signalLast - Iterator to one past the last signal sample.
     * \param[out] resultFirst - Iterator to start of correctly sized result container.
     *
     * The signal samples should be from a uniformly sampled source.
     *
     * This class will perform upsampling and downsampling as specifed and as required.
     *
     * The number of resampled samples returned in a vector will be:
     * floor(N * U / D), where N is number of signal samples, U is the upsample factor
     * and D is the downsample factor.
     *
     * The sample rate of the resampled data is equivalent to:
     * S * U / D, where S is the signals original sample rate in Hz.
     *
     * To compute S to allocate a correctly size result vector call
     * ResampledSize method.
     */
    template <typename IterIn, typename OutIter>
    void operator()(IterIn signalFirst, IterIn signalLast, OutIter resultFirst)
    {
        auto signalLength = std::distance(signalFirst, signalLast);
        DSP_ASSERT_THROW(signalLength == static_cast<decltype(signalLength)>(m_signalLength),
                         "sample length is incorrect");

        if (m_upsampleFactor > 1)
        {
            // Make sure upsample buffer starts zeroed.
            std::fill(
                m_workspaceBuffer.begin(), m_workspaceBuffer.end(), static_cast<FloatType>(0));

            // Fill upsample buffer with signal samples with zero padding.
            auto upIter = m_workspaceBuffer.begin();

            for (auto itr = signalFirst; (itr != signalLast) && (upIter < m_workspaceBuffer.end());
                 std::advance(itr, 1), std::advance(upIter, static_cast<int>(m_upsampleFactor)))
            {
                // Note: we correct for attenuation by upsample factor, caused by inserting zeroes.
                *upIter = *itr * static_cast<FloatType>(m_upsampleFactor);
            }

            // Low pass filter the workspace buffer.
            m_filterHolder(m_workspaceBuffer.begin(),
                           m_workspaceBuffer.end(),
                           m_workspaceBuffer.begin(),
                           true);

            // If required perform downsampling.
            if (m_downsampleFactor > 1)
            {
                // Fill the result buffer by dropping samples as required.
                auto upIter     = m_workspaceBuffer.begin();
                auto resultIter = resultFirst;

                for (auto i = 0; (i < m_resampledLength) && (upIter < m_workspaceBuffer.end()); ++i,
                          std::advance(upIter, static_cast<int>(m_downsampleFactor)),
                          std::advance(resultIter, 1))
                {
                    *resultIter = *upIter;
                }
            }
            // Else return upsample buffer as is copied ito the result buffer.
            else
            {
                std::copy(m_workspaceBuffer.begin(), m_workspaceBuffer.end(), resultFirst);
            }
        }
        else
        {
            // Low pass filter the signal data.
            m_filterHolder(signalFirst, signalLast, m_workspaceBuffer.begin(), true);

            // Fill the result buffer by dropping samples as required.
            auto   resultIter    = resultFirst;
            size_t resampleCount = 0;

            for (auto itr = m_workspaceBuffer.begin();
                 (itr < m_workspaceBuffer.end()) && (resampleCount < m_resampledLength);
                 std::advance(itr, static_cast<int>(m_downsampleFactor)),
                      std::advance(resultIter, 1),
                      ++resampleCount)
            {
                *resultIter = *itr;
            }
        }
    }

private:
    /*! \brief Compute the resampled data size. */
    void ComputeResampledSize()
    {
        auto upsampleLength = m_upsampleFactor * m_signalLength;
        m_resampledLength =
            static_cast<size_t>(std::floor((static_cast<double>(upsampleLength) /
                                            static_cast<double>(m_downsampleFactor))) +
                                0.5);
    }
    /*! \brief Resize the workspace buffers. */
    void ResizeWorkspace()
    {
        if (m_upsampleFactor > 1)
        {
            auto upsampleLength = m_upsampleFactor * m_signalLength;
            m_workspaceBuffer.resize(upsampleLength);
        }
        else
        {
            m_workspaceBuffer.resize(m_signalLength);
        }
    }

private:
    /*! \brief The number of signal samples to resample. */
    size_t m_signalLength{0};
    /*! \brief The upsample factor. */
    size_t m_upsampleFactor{1};
    /*! \brief The downsample factor. */
    size_t m_downsampleFactor{1};
    /*! \brief The number of samples in the resampled result. */
    size_t m_resampledLength{0};
    /*! \brief Filter holder for our kaiser window low pass FIR filter. */
    filter_hldr_t m_filterHolder{};
    /*! \brief Workspace buffer. */
    std::vector<FloatType> m_workspaceBuffer{};
};

/*! \brief Convenience typedef to Resample<float>. */
using resample_f = Resample<float>;
/*! \brief Convenience typedef to Resample<double>. */
using resample_d = Resample<double>;
/*! \brief Convenience typedef to Resample<long double>. */
using resample_ld = Resample<long double>;

} // namespace dsp

#endif // DSP_RESAMPLE_HPP
