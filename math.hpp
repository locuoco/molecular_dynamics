#ifndef MATH_H
#define MATH_H

#include <limits>
#include <cmath>

#include <GL/glew.h>

#include <glm/ext.hpp>

namespace math
{
	template <typename T>
	constexpr T pi()
	{
		return (T)3.1415926535897932384626433832795L;
	}

	template <typename T>
	constexpr T half_pi()
	{
		return (T)1.5707963267948966192313216916398L;
	}

	template <typename T>
	constexpr T pi_4()
	{
		return (T)0.78539816339744830961566084581988L;
	}

	template <typename T>
	constexpr T two_pi()
	{
		return (T)6.283185307179586476925286766559L;
	}

	template <typename T>
	constexpr T sqrt2()
	{
		return (T)1.4142135623730950488016887242097L;
	}

	template <typename T>
	constexpr T sqrt2_()
	{
		return (T)0.70710678118654752440084436210485L;
	}

	template <typename T>
	constexpr T eps()
	{
		return std::numeric_limits<T>::epsilon();
	}

	template <typename T>
	constexpr T one_minus_eps()
	{
		return 1 - std::numeric_limits<T>::epsilon();
	}

	template <typename T>
	constexpr T one_plus_eps()
	{
		return 1 + std::numeric_limits<T>::epsilon();
	}

	template <typename T>
	void quadratic_solver(T &x1, T &x2, T a, T b, T c)
	{
		using std::sqrt;
		T u = -b - glm::sign(b) * sqrt(b * b - 4 * a * c);
		x1 = u / (2 * a);
		x2 = 2 * c / u;
	}

	template <typename T>
	constexpr T deg2rad(T deg)
	{
		return deg*pi<T>()/180;
	}

	template <typename T>
	constexpr T rad2deg(T rad)
	{
		return rad/pi<T>()*180;
	}
}

#endif




















