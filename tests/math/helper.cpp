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
#include <cmath> // sin, cos, sqrt, abs, exp, erfc, erf
#include <random> // mt19937_64
#include <complex>

/*

Compilation:
g++ helper.cpp -o helper -std=c++20 -Wall -Wextra -pedantic -Ofast -fmax-errors=1

*/

#include "../../math/helper.hpp"

constexpr double error_threshold = 1e-15;

std::mt19937_64 mersenne_twister(1234); // seed set to 1234
// mersenne_twister will behave in the same way for all compilers/runs

void test_pow_of_2()
// test that `is_pow_of_2` check powers of 2 correctly for some arguments.
// Since the assertion is checked at compile time and the function is already
// instantiated (it is not templated), the function does not need to be called.
{
	static_assert(math::is_pow_of_2(0) == false);
	static_assert(math::is_pow_of_2(1) == true);
	static_assert(math::is_pow_of_2(2) == true);
	static_assert(math::is_pow_of_2(3) == false);
	static_assert(math::is_pow_of_2(16) == true);
	static_assert(math::is_pow_of_2(257) == false);
	// ~0000 = 1111,  1111 >> 1 = 0111,  0111 + 1 = 1000 -> a power of 2
	static_assert(math::is_pow_of_2((~0ull >> 1) + 1) == true);
}

void test_powm1(int n)
// test math::powm1(n) = (-1)^n
{
	int res = (n % 2 == 0) ? 1 : -1;
	assert(math::powm1(n) == res);
}

void test_rms()
{
	std::valarray<double> x;
	std::valarray<std::complex<double>> z;
}

void test_deg2rad()
// unit test showing math::deg2rad(360.) == 2 pi
{
	assert(std::abs(math::deg2rad(360.) - math::two_pi<double>) < error_threshold);
}

void test_rad2deg()
// unit test showing math::rad2deg(pi / 4) == 45
{
	assert(std::abs(math::rad2deg(math::pi_4<double>) - 45) < error_threshold);
}

void test_fastexp(double x, double max_relerr)
// test relative error of `math::fastexp` against `std::exp` for argument `x`.
// `max_relerr` is the maximum relative error tolerated for this test.
{
	double ref = std::exp(x);
	double relerr = std::abs((math::fastexp(x) - ref)/ref);
	assert(relerr <= max_relerr);
}

void test_fasterfc(double x)
// test absolute error of `math::fasterfc` against `std::erfc` for argument `x`.
{
	assert(std::abs(math::fasterfc(x) - std::erfc(x)) < 1e-6);
}

void test_fasterf(double x)
// test absolute error of `math::fasterf` against `std::erf` for argument `x`.
{
	assert(std::abs(math::fasterf(x) - std::erf(x)) < 1e-6);
}

void test_mod1()
// test that math::mod(-0.5, 4.) == 3.5
{
	assert(math::mod(-0.5, 4.) == 3.5);
}

void test_mod2()
// test that math::mod(-3, 4) == 1.
// Since the assertion is checked at compile time and the function is already
// instantiated (it is not templated), the function does not need to be called.
{
	static_assert(math::mod(-3, 4) == 1);
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
// test that math::clamp(-1, 0, 1) == 0.
// Since the assertion is checked at compile time and the function is already
// instantiated (it is not templated), the function does not need to be called.
{
	static_assert(math::clamp(-1, 0, 1) == 0);
}

void test_clamp2()
// test that math::clamp(2, -3, 1) == 1.
// Since the assertion is checked at compile time and the function is already
// instantiated (it is not templated), the function does not need to be called.
{
	static_assert(math::clamp(2., -3., 1.) == 1);
}

int main()
{
	test_powm1(2);
	test_powm1(5);
	test_powm1(9);
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
	test_sinc_denormal();

	std::cout << "All tests passed successfully!" << std::endl;

	return 0;
}
























