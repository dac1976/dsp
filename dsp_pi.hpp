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
 * \file dsp_pi.hpp
 * \brief File containing generic accurate definitions for Pi for any floating point type.
 */
#ifndef DSP_PI_HPP
#define DSP_PI_HPP

#include <type_traits>
#include <cmath>

/*! \brief dsp namespace */
namespace dsp
{

/*!
 * \brief Function to compute Pi using acos(-1).
 * \return Result represented as a floating point value of type T.
 *
 * This function only has computational cost when first called. All
 * subsequent calls simply return the static const result from
 * the first time the function was called.
 *
 * If called with anything other than a flaoting point type for T
 * then a static_assert will fire reporting an invalid type.
 */
template <typename FloatType> FloatType Pi()
{
    static_assert(std::is_floating_point<FloatType>::value, "invalid floating point type");
    static const auto result = std::acos(static_cast<FloatType>(-1));
    return result;
}

/*!
 * \brief Function to compute Pi/2.
 * \return Result represented as a floating point value of type T.
 *
 * This function only has computational cost when first called. All
 * subsequent calls simply return the static const result from
 * the first time the function was called.
 *
 * If called with anything other than a flaoting point type for T
 * then a static_assert will fire reporting an invalid type.
 */
template <typename FloatType> FloatType HalfPi()
{
    static const auto result = Pi<FloatType>() / static_cast<FloatType>(2);
    return result;
}

/*!
 * \brief Function to compute Pi/4.
 * \return Result represented as a floating point value of type T.
 *
 * This function only has computational cost when first called. All
 * subsequent calls simply return the static const result from
 * the first time the function was called.
 *
 * If called with anything other than a flaoting point type for T
 * then a static_assert will fire reporting an invalid type.
 */
template <typename FloatType> FloatType QuarterPi()
{
    static const auto result = HalfPi<FloatType>() / static_cast<FloatType>(2);
    return result;
}

/*!
 * \brief Function to compute 2Pi.
 * \return Result represented as a floating point value of type T.
 *
 * This function only has computational cost when first called. All
 * subsequent calls simply return the static const result from
 * the first time the function was called.
 *
 * If called with anything other than a flaoting point type for T
 * then a static_assert will fire reporting an invalid type.
 */
template <typename FloatType> FloatType TwoPi()
{
    static const auto result = static_cast<FloatType>(2) * Pi<FloatType>();
    return result;
}

/*!
 * \brief Function to compute 1/Pi.
 * \return Result represented as a floating point value of type T.
 *
 * This function only has computational cost when first called. All
 * subsequent calls simply return the static const result from
 * the first time the function was called.
 *
 * If called with anything other than a flaoting point type for T
 * then a static_assert will fire reporting an invalid type.
 */
template <typename FloatType> FloatType OneOverPi()
{
    static const auto result = static_cast<FloatType>(1) / Pi<FloatType>();
    return result;
}

/*!
 * \brief Function to compute 2/Pi.
 * \return Result represented as a floating point value of type T.
 *
 * This function only has computational cost when first called. All
 * subsequent calls simply return the static const result from
 * the first time the function was called.
 *
 * If called with anything other than a flaoting point type for T
 * then a static_assert will fire reporting an invalid type.
 */
template <typename FloatType> FloatType TwoOverPi()
{
    static const auto result = static_cast<FloatType>(2) * OneOverPi<FloatType>();
    return result;
}

/*!
 * \brief Function to compute 3*Pi/2.
 * \return Result represented as a floating point value of type T.
 *
 * This function only has computational cost when first called. All
 * subsequent calls simply return the static const result from
 * the first time the function was called.
 *
 * If called with anything other than a flaoting point type for T
 * then a static_assert will fire reporting an invalid type.
 */
template <typename FloatType> FloatType ThreeOverTwoPi()
{
    static const auto result = static_cast<FloatType>(3. / 2.) * Pi<FloatType>();
    return result;
}

/*!
 * \brief Function to compute 2/sqrt(Pi).
 * \return Result represented as a floating point value of type T.
 *
 * This function only has computational cost when first called. All
 * subsequent calls simply return the static const result from
 * the first time the function was called.
 */
template <typename FloatType> FloatType TwoOverSqrtPi()
{
    static const auto result = static_cast<FloatType>(2) / std::sqrt(Pi<FloatType>());
    return result;
}

} // namespace dsp

#endif // DSP_PI_HPP
