//  Math helper functions tests
//  Copyright (C) 2022 Alessandro Lo Cuoco (alessandro.locuoco@gmail.com)

//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.

//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.

//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <iostream> // cout, endl
#include <cassert>
#include <cmath> // sin, cos, sqrt, fabs, exp, erfc, erf

/*

Compilation:
g++ helper.cpp -o helper -std=c++20 -Wall -Wextra -pedantic -Ofast -pthread -fmax-errors=1

*/

#include "../../math/helper.hpp"

constexpr double error_threshold = 1e-15;

void test_deg2rad()
// unit test showing math::deg2rad(360.) == 2 pi
{
	using std::fabs;
	assert(fabs(math::deg2rad(360.) - math::two_pi<double>) < error_threshold);
}

void test_rad2deg()
// unit test showing math::rad2deg(pi / 4) == 45
{
	using std::fabs;
	assert(fabs(math::rad2deg(math::pi_4<double>) - 45) < error_threshold);
}

void test_fastexp(double x, double max_relerr)
// test relative error of `math::fastexp` against `std::exp` for argument `x`.
// `max_relerr` is the maximum relative error tolerated for this test.
{
	using std::exp;
	double ref = exp(x);
	double relerr = fabs((math::fastexp(x) - ref)/ref);
	assert(relerr <= max_relerr);
}

void test_fasterfc(double x)
// test absolute error of `math::fasterfc` against `std::erfc` for argument `x`.
{
	using std::erfc;
	assert(fabs(math::fasterfc(x) - erfc(x)) < 1e-6);
}

void test_fasterf(double x)
// test absolute error of `math::fasterf` against `std::erf` for argument `x`.
{
	using std::erf;
	assert(fabs(math::fasterf(x) - erf(x)) < 1e-6);
}

void test_mod1()
// test that math::mod(-0.5, 4.) == 3.5
{
	assert(math::mod(-0.5, 4.) == 3.5);
}

void test_mod2()
// test that math::mod(-3, 4) == 1
{
	assert(math::mod(-3, 4) == 1);
}

void test_sinc_denormal()
// unit test showing that `math::sinc` of a denormal is 1.
// Note that sinc(x) ~ 1 - x^2 / 6  for x ~ 0. Since floating point arithmetic
// has finite precision, for any `x` smaller than a threshold, 1 is the correctly
// rounded result.
// This test should pass even when compiled with -Ofast in GCC
{
	assert(math::sinc(1.e-310) == 1);
}

void test_clamp1()
// test that math::clamp(-1, 0, 1) == 0
{
	assert(math::clamp(-1, 0, 1) == 0);
}

void test_clamp2()
// test that math::clamp(2, -3, 1) == 1
{
	assert(math::clamp(2., -3., 1.) == 1);
}

int main()
{
	test_deg2rad();
	test_rad2deg();
	test_fastexp(1e-6, 1e-7);
	test_fastexp(-.01, 1e-6);
	test_fastexp(.1, 1e-6);
	test_fastexp(1, 1e-8);
	test_fastexp(2, 1e-7);
	test_fastexp(-10, 1e-7);
	// fastexp becomes less and less accurate far from 0
	test_fastexp(-200, 1e-5);
	test_fastexp(100, 1e-5);
	test_fasterfc(0);
	test_fasterfc(-1);
	test_fasterfc(10);
	test_fasterf(0);
	test_fasterf(-1);
	test_fasterf(10);
	test_mod1();
	test_mod2();
	test_sinc_denormal();
	test_clamp1();
	test_clamp2();

	std::cout << "All tests passed successfully!" << std::endl;

	return 0;
}
























