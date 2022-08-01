//  Math helper functions
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

#ifndef MATH_HELPER_H
#define MATH_HELPER_H

#include <concepts> // floating_point, assignable_from
#include <cmath> // sqrt, fmod, sin
#include <type_traits> // is_integral_v
#include <numbers> // numbers::pi_v
#include <valarray>
#include <complex>

namespace math
{
	template <typename T>
	inline constexpr T half_pi = 1.5707963267948966192313216916398L; // pi/2

	template <typename T>
	inline constexpr T pi_4 = 0.78539816339744830961566084581988L; // pi/4

	template <typename T>
	inline constexpr T two_pi = 6.283185307179586476925286766559L; // 2 pi

	template <typename T>
	inline constexpr T sqrtpi = 1.7724538509055160272981674833411L; // square root of pi

	template <typename T>
	inline constexpr T cbrtpi = 1.4645918875615232630201425272638L; // cube root of pi

	template <typename T>
	inline constexpr T inv_sqrt2 = 0.70710678118654752440084436210485L; // 1/sqrt(2)

	template <typename T>
	inline constexpr T two1_6 = 1.1224620483093729814335330496792L; // 2^(1/6)

	template <typename T>
	inline constexpr bool is_complex = false;
	
	template <typename T>
	inline constexpr bool is_complex<std::complex<T>> = true;

	template <std::integral T>
	constexpr bool is_even(T num) noexcept
	// return true if and only if `num` is even
	{
		return (num & 1) == 0;
		// if the least significant bit is 0 the number is even.
	}
	template <std::integral T>
	constexpr bool is_odd(T num) noexcept
	// return true if and only if `num` is odd
	{
		return !is_even(num);
	}

	constexpr bool is_pow_of_2(unsigned long long num) noexcept
	// return true if and only if `num` is a power of 2
	{
		return (num & (num-1)) == 0 && num != 0;
		// the preceding number of any number in binary has at least one common bit
		// (the msb) unless it is 0 or a power of 2.
	}

	constexpr int powm1(int expo) noexcept
	// calculate (-1)^expo (power of minus 1)
	{
		return 1-2*(expo&1);
	}

	template <typename Vec>
	requires (std::assignable_from<std::valarray<typename Vec::value_type>&, Vec> && std::floating_point<typename Vec::value_type>)
	auto rms(const Vec& v)
	// calculate root mean square of a valarray of real numbers `v`
	{
		using std::sqrt;
		return sqrt((v*v).sum()/v.size());
	}

	template <typename Vec>
	requires (std::assignable_from<std::valarray<typename Vec::value_type>&, Vec> && is_complex<typename Vec::value_type>)
	auto rms(const Vec& v)
	// calculate root mean square of a valarray of complex numbers `v`, using the formula:
	// rms(v) = sqrt( sum_i |v_i|^2 / n)
	// where n is the number of complex numbers.
	{
		using std::sqrt;
		auto sqr = [](const auto& z) { return std::complex(norm(z), 0.); };
		return sqrt(v.apply(sqr).sum().real()/v.size());
	}

	template <std::floating_point T>
	constexpr T deg2rad(T deg) noexcept
	// convert degrees to radians
	{
		return deg*std::numbers::pi_v<T>/180;
	}

	template <std::floating_point T>
	constexpr T rad2deg(T rad) noexcept
	// convert radians to degrees
	{
		return rad*180/std::numbers::pi_v<T>;
	}

	template <std::floating_point T, std::size_t Exponent = sizeof(T)*4>
	constexpr T fastexp(T x) noexcept
	// based on the fundamental limit:
	// (1 + x/n)^n -> e^x for n -> inf
	// and exponentiation by squaring,
	// with n = 2^Exponent
	// if `Exponent` is too big, underflows/truncation errors are likely to occur
	// if it is too small, the result will be inaccurate for big values of |x|.
	// These limitations make this algorithm not suited to obtain a correct
	// result up to machine precision, but it is much faster than std::exp
	{
		x = 1 + x / (1ull << Exponent);
		for (std::size_t i = 0; i < Exponent; ++i)
			x *= x;
		return x;
	}

	template <std::floating_point T>
	constexpr T fasterfc(T x) noexcept
	// From Abramowitz and Stegun (1964)
	// calculate the complementary error function
	// max relative error should be around 10^-7
	{
		if (x < 0)
			return 2 - fasterfc(-x);
		T res = 1;
		res += T(0.0705230784L)*x;
		T xp = x*x;
		res += T(0.0422820123L)*xp;
		xp *= x;
		res += T(0.0092705272L)*xp;
		xp *= x;
		res += T(0.0001520143L)*xp;
		xp *= x;
		res += T(0.0002765672L)*xp;
		xp *= x;
		res += T(0.0000430638L)*xp;
		res *= res;
		res *= res;
		res *= res;
		res *= res;
		return 1/res;
	}

	template <std::floating_point T>
	constexpr T fasterf(T x) noexcept
	// calculate the error function
	// max relative error should be around 10^-7
	{
		return 1 - fasterfc(x);
	}

	template <typename T>
	constexpr T mod(T x, T y) noexcept
	// remaps x in the [0, y) range, i.e.:
	//	res = x mod y
	// Note that x % y = n such that x = trunc(x/y)*y + n.
	// On the other hand, x mod y = n such that x = floor(x/y)*y + n.
	// The value returned by this function is thus always positive
	// (fmod behaves as % but with floating-point types)
	{
		using std::fmod;
		if constexpr (std::is_integral_v<T>)
			return (x %= y) < 0 ? x+y : x;
		else
			return (x = fmod(x, y)) < 0 ? x+y : x;
	}

	template <typename T>
	T sinc(T x) noexcept
	// calculate cardinal sine function:
	// sinc(0) = 1
	// sinc(x) = sin(x)/x  for x != 0
	{
		using std::sin;
		return x*x == 0 ? 1 : sin(x)/x;
		// by multiplying x by itself, I force denormals to flush to zero
		// In this case, 1 is the correct result
		// (if x is a denormal (or zero) and -ffast-math is enabled, then this
		// avoids division by 0)
	}

	template <typename T>
	constexpr T clamp(T x, T a, T b) noexcept
	// clamp `x` between `a` and `b`
	{
		return x>b?b:x<a?a:x;
	}

}

#endif // MATH_HELPER_H




















