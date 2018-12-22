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

#ifndef TEST_TIMER_HPP
#define TEST_TIMER_HPP

#include <chrono>

/*! \brief namespace test */
namespace test
{

/*! \brief Timer class to help with test benchmarking.*/
class Timer final
{
    /*! \brief Typedef to high resolution clock.*/
    using high_res_clock_t = std::chrono::high_resolution_clock;
    /*! \brief Typedef to high resolution clock time point.*/
    using time_point_t     = std::chrono::high_resolution_clock::time_point;
    /*! \brief Typedef to timer duration.*/
    using duration_t       = std::chrono::duration<double>;

public:
    /*!
     * \brief Compute elapsed time of timer.
     * \return Elapsed time of the timer as a double.
     */
    double Elapsed() const
    {
        return std::chrono::duration_cast<duration_t>(high_res_clock_t::now() - m_begin).count();
    }

    /*! \brief Reset the timers start time point.*/
    void Reset()
    {
        m_begin = high_res_clock_t::now();
    }

private:
	/*! \brief Initialise begin time point to now as soon as Timer is created.*/
    time_point_t m_begin{high_res_clock_t::now()};
};

} // namespace test

#endif // TEST_TIMER_HPP
