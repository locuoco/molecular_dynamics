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

#include <limits> // numeric_limits
#include <cmath> // sqrt, fmod, sin
#include <type_traits> // is_integral_v

namespace math
{
	template <typename T>
	constexpr T pi() noexcept // in C++20 we have std::numbers::pi_v<T>
	{
		return (T)3.1415926535897932384626433832795L;
	}

	template <typename T>
	constexpr T half_pi() noexcept
	{
		return (T)1.5707963267948966192313216916398L;
	}

	template <typename T>
	constexpr T pi_4() noexcept // pi/4
	{
		return (T)0.78539816339744830961566084581988L;
	}

	template <typename T>
	constexpr T two_pi() noexcept
	{
		return (T)6.283185307179586476925286766559L;
	}
	template <typename T>
	constexpr T sqrtpi() noexcept
	{
		return (T)1.7724538509055160272981674833411L;
	}

	template <typename T>
	constexpr T sqrt2() noexcept // in C++20 we have std::numbers::sqrt2_v<T>
	{
		return (T)1.4142135623730950488016887242097L;
	}

	template <typename T>
	constexpr T sqrt2_() noexcept // 1/sqrt(2)
	{
		return (T)0.70710678118654752440084436210485L;
	}

	template <typename T>
	constexpr T eps() noexcept
	{
		return std::numeric_limits<T>::epsilon();
	}

	template <typename T>
	constexpr T one_minus_eps() noexcept
	{
		return 1 - std::numeric_limits<T>::epsilon();
	}

	template <typename T>
	constexpr T one_plus_eps() noexcept
	{
		return 1 + std::numeric_limits<T>::epsilon();
	}

	template <typename T>
	void quadratic_solver(T &x1, T &x2, T a, T b, T c) noexcept
	{
		using std::sqrt;
		T sign = b > 0 ? 1 : b < 0 ? -1 : 0;
		T u = -b - sign * sqrt(b * b - 4 * a * c);
		x1 = u / (2 * a);
		x2 = 2 * c / u;
	}

	template <typename T>
	constexpr T deg2rad(T deg) noexcept
	{
		return deg*pi<T>()/180;
	}

	template <typename T>
	constexpr T rad2deg(T rad) noexcept
	{
		return rad/pi<T>()*180;
	}

	template <typename T, std::size_t exponent = sizeof(T)*4>
	T fastexp(T x) noexcept
	// based on the fundamental limit:
	// (1 + x/n)^n -> e^x for n -> inf
	// and exponentiation by squaring,
	// with n = 2^exponent
	// if exponent is too big, underflows/truncation errors are likely to occur
	// if it is too small, the result will be inaccurate for big values of |x|.
	// These limitations make this algorithm not suited to obtain a correct
	// result up to machine precision, but it is really fast
	{
		x = 1 + x / (1ull << exponent);
		for (std::size_t i = 0; i < exponent; ++i)
			x *= x;
		return x;
	}

	template <typename T>
	T fasterfc(T x) noexcept
	// From Abramowitz and Stegun (1964)
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

	template <typename T>
	T fasterf(T x) noexcept
	{
		return 1 - fasterfc(x);
	}

	template <typename T>
	T mod(T x, T y) noexcept
	// Note that x % y = n such that x = trunc(x/y)*y + n
	// on the other hand, x mod y = n such that x = floor(x/y)*y + n
	// the value returned by this function is thus always positive
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
	{
		using std::sin;
		return x*x == 0 ? 1 : sin(x)/x;
		// by multiplying x by itself, I force denormals to flush to zero
		// In this case, 1 is the correct result
		// (if x is a denormal (or zero) and -ffast-math is enabled, then this
		// avoids division by 0)
	}

	template <typename T>
	T clamp(T x, T a, T b) noexcept
	// clamp x between a and b
	{
		return x>b?b:x<a?a:x; // c:
	}

}

#endif // MATH_HELPER_H




















