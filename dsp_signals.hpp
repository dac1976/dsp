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
 * \file dsp_signals.hpp
 * \brief File containing generic signal generators.
 */

#ifndef DSP_SIGNALS_HPP
#define DSP_SIGNALS_HPP

#include "dsp_math.hpp"

/*! \brief dsp namespace */
namespace dsp
{

/*! \brief Struct holding params to define a sinusoidal tone. */
struct ToneParams
{
    /*! \brief Peak amplitude of sine wave. */
    double amplitude{0};
    /*! \brief Signal frequency in Hz. */
    double frequency{0};
    /*! \brief Phase offset in radians. */
    double phase{0};
    /*! \brief Amplitude (DC) offset. */
    double offset{0};
};

/*!
 * \brief Sinusoidal single tone generator.
 * \param[in] params - Struct holding params to define a sinusoidal signal.
 * \param[in] sampleRate - Sample rate in Hz of generated samples.
 * \param[in] count - Number of samples to generate.
 * \return Series of samples representing a sinusoisal signal.
 *
 * The sampleRate should be >= 2*frequency so that Nyquist's
 * theorem holds.
 */
template <typename FloatType>
std::vector<FloatType> Tone(ToneParams const& params, double sampleRate, size_t count)
{
    static_assert(std::is_floating_point<FloatType>::value, "invalid floating point type");
    std::vector<FloatType> samples(count, static_cast<FloatType>(0));
    auto                   timeInterval = static_cast<FloatType>(1. / sampleRate);
    auto                   time         = static_cast<FloatType>(0);

    for (auto& sample : samples)
    {
        sample = Sine(static_cast<FloatType>(params.amplitude),
                      time,
                      static_cast<FloatType>(params.frequency),
                      static_cast<FloatType>(params.phase),
                      static_cast<FloatType>(params.offset));
        time += timeInterval;
    }

    return samples;
}

/*!
 * \brief Sinusoidal multi-tone generator.
 * \param[in] allParams - Struct holding params to define a sinusoidal signal.
 * \param[in] sampleRate - Sample rate in Hz of generated samples.
 * \param[in] count - Number of samples to generate.
 * \return Series of samples representing a sinusoisal signal.
 *
 * The sampleRate should be >= 2*max(frequency) so that Nyquist's
 * theorem holds.
 */
template <typename FloatType>
std::vector<FloatType> MultiTone(std::vector<ToneParams> const& allParams, double sampleRate,
                                 size_t count)
{
    static_assert(std::is_floating_point<FloatType>::value, "invalid floating point type");
    std::vector<FloatType> samples(count, static_cast<FloatType>(0));

    auto timeInterval = static_cast<FloatType>(1. / sampleRate);
    auto time         = static_cast<FloatType>(0);

    for (auto& sample : samples)
    {
        for (auto const& params : allParams)
        {
            sample += Sine(static_cast<FloatType>(params.amplitude),
                           time,
                           static_cast<FloatType>(params.frequency),
                           static_cast<FloatType>(params.phase),
                           static_cast<FloatType>(params.offset));
        }

        time += timeInterval;
    }

    return samples;
}

} // namespace dsp

#endif // DSP_SIGNALS_HPP
