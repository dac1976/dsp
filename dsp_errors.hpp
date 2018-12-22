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
 * \file dsp_errors.hpp
 * \brief File containing error related macros, functions definitions.
 */

#ifndef DSP_ERRORS_HPP
#define DSP_ERRORS_HPP

#include <cassert>
#include <stdexcept>

/*! \brief Macro to use to assert a value.*/
#define DSP_ASSERT_THROW(x, s)                                                                     \
    do                                                                                             \
    {                                                                                              \
        assert(x);                                                                                 \
        if (!(x))                                                                                  \
            throw std::runtime_error(s);                                                           \
    } while (false)

#endif // DSP_ERRORS_HPP
