//  Data structures for use in other parts of the program
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

#ifndef PHYSICS_POINT_H
#define PHYSICS_POINT_H

#include <concepts> // floating_point, integral
#include <valarray>
#include <vector>
#include <array>
#include <utility> // index_sequence
#include <cmath> // sqrt

namespace physics
{
	namespace detail_
	{
		template <typename R, typename T, std::size_t ... Ns>
		constexpr R make_filled_impl(T scal, std::index_sequence<Ns...>) noexcept
		{
			return {(static_cast<void>(Ns), scal)...};
		}
	}
	template <typename R, std::size_t N, typename T>
	constexpr R make_filled(T scal) noexcept
	{
		return detail_::make_filled_impl<R, T>(scal, std::make_index_sequence<N>{});
	}

	template <std::floating_point T, std::size_t M, std::size_t N = M>
	requires (M >= 1 && N >= 1)
	struct mat;

	template <typename T, std::size_t Dim>
	using point = mat<T, Dim, 1>;

	template <std::floating_point T, std::size_t Dim>
	constexpr point<T, Dim> make_point(T scal) noexcept
	{
		return make_filled<point<T, Dim>, Dim>(scal);
	}

	template <std::floating_point T, std::size_t M, std::size_t N>
	requires (M >= 1 && N >= 1)
	struct mat : std::array<T, M*N>
	{
		private:
			using base = std::array<T, M*N>;

		public:

		constexpr mat() noexcept = default;

		template <typename ... Ts>
		requires (sizeof...(Ts) <= M*N && sizeof...(Ts) >= 2 && ((std::floating_point<Ts> || std::integral<Ts>) && ...))
		constexpr mat(Ts ... ts) noexcept : base{static_cast<T>(ts)...} {}

		constexpr mat(T t) noexcept : base(make_filled<mat<T, M, N>, M*N>(t)) {}

		constexpr T& operator()(std::size_t i, std::size_t j) noexcept
		{
			return base::operator[](i*N + j);
		}

		constexpr const T& operator()(std::size_t i, std::size_t j) const noexcept
		{
			return base::operator[](i*N + j);
		}

		constexpr const mat& operator+() const noexcept
		{
			return *this;
		}
		constexpr mat operator-() const noexcept
		{
			mat res;
			for (std::size_t i = 0; i < M*N; ++i)
				res[i] = -base::operator[](i);
			return res;
		}
		constexpr mat& operator+=(const mat& other) noexcept
		{
			for (std::size_t i = 0; i < M*N; ++i)
				base::operator[](i) += other[i];
			return *this;
		}
		constexpr mat& operator-=(const mat& other) noexcept
		{
			for (std::size_t i = 0; i < M*N; ++i)
				base::operator[](i) -= other[i];
			return *this;
		}
		constexpr mat& operator*=(const mat& other) noexcept
		{
			for (std::size_t i = 0; i < M*N; ++i)
				base::operator[](i) *= other[i];
			return *this;
		}
		constexpr mat& operator/=(const mat& other) noexcept
		{
			for (std::size_t i = 0; i < M*N; ++i)
				base::operator[](i) /= other[i];
			return *this;
		}
		constexpr mat& operator*=(T scal) noexcept
		{
			for (std::size_t i = 0; i < M*N; ++i)
				base::operator[](i) *= scal;
			return *this;
		}
		constexpr mat& operator/=(T scal) noexcept
		{
			for (std::size_t i = 0; i < M*N; ++i)
				base::operator[](i) /= scal;
			return *this;
		}
	};

	template <std::floating_point T, std::size_t M, std::size_t N>
	constexpr mat<T, M, N> operator+(mat<T, M, N> lhs, const mat<T, M, N>& rhs) noexcept
	{
		return lhs += rhs;
	}
	template <std::floating_point T, std::size_t M, std::size_t N>
	constexpr mat<T, M, N> operator-(mat<T, M, N> lhs, const mat<T, M, N>& rhs) noexcept
	{
		return lhs -= rhs;
	}
	template <std::floating_point T, std::size_t M, std::size_t N>
	constexpr mat<T, M, N> operator*(mat<T, M, N> lhs, const mat<T, M, N>& rhs) noexcept
	{
		return lhs *= rhs;
	}
	template <std::floating_point T, std::size_t M, std::size_t N>
	constexpr mat<T, M, N> operator/(mat<T, M, N> lhs, const mat<T, M, N>& rhs) noexcept
	{
		return lhs /= rhs;
	}
	template <std::floating_point T, std::size_t M, std::size_t N>
	constexpr mat<T, M, N> operator*(mat<T, M, N> lhs, typename mat<T, M, N>::value_type scal) noexcept
	{
		return lhs *= scal;
	}
	template <std::floating_point T, std::size_t M, std::size_t N>
	constexpr mat<T, M, N> operator*(typename mat<T, M, N>::value_type scal, mat<T, M, N> rhs) noexcept
	{
		return rhs * scal;
	}
	template <std::floating_point T, std::size_t M, std::size_t N>
	constexpr mat<T, M, N> operator/(mat<T, M, N> lhs, typename mat<T, M, N>::value_type scal) noexcept
	{
		return lhs /= scal;
	}
	template <std::floating_point T, std::size_t M, std::size_t N>
	constexpr mat<T, M, N> operator/(typename mat<T, M, N>::value_type scal, mat<T, M, N> rhs) noexcept
	{
		return rhs / scal;
	}

	template <std::floating_point T, std::size_t M, std::size_t N, std::size_t K>
	constexpr mat<T, M, K> operator%(const mat<T, M, N>& lhs, const mat<T, N, K>& rhs) noexcept
	{
		using std::size_t;
		mat<T, M, K> res(0);
		for (size_t i = 0; i < M; ++i)
			for (size_t j = 0; j < N; ++j)
				for (size_t k = 0; k < K; ++k)
					res(i, k) += lhs(i, j) * rhs(j, k);
		return res;
	}

	template <std::floating_point T, std::size_t M, std::size_t N>
	constexpr mat<T, M, N> remainder(mat<T, M, N> A, T modulo) noexcept
	{
		using std::round;
		for (std::size_t i = 0; i < M*N; ++i)
			A[i] = A[i] - round(A[i]/modulo)*modulo;
		return A;
	}

	template <std::floating_point T, std::size_t M, std::size_t N>
	constexpr mat<T, M, N> outer(const point<T, M>& lhs, const point<T, N>& rhs) noexcept
	{
		using std::size_t;
		mat<T, M, N> res;
		for (size_t i = 0; i < M; ++i)
			for (size_t j = 0; j < N; ++j)
				res(i, j) = lhs[i] * rhs[j];
		return res;
	}

	template <std::floating_point T, std::size_t Dim>
	constexpr T trace(mat<T, Dim> A) noexcept
	{
		T res(0);
		for (std::size_t i = 0; i < Dim; ++i)
			res += A(i, i);
		return res;
	}

	template <std::floating_point T, std::size_t Dim>
	constexpr T dot(const point<T, Dim>& lhs, const point<T, Dim>& rhs) noexcept
	{
		T res(0);
		for (std::size_t i = 0; i < Dim; ++i)
			res += lhs[i] * rhs[i];
		return res;
	}

	template <std::floating_point T, std::size_t Dim>
	constexpr T norm(const point<T, Dim>& pnt) noexcept
	{
		using std::sqrt;
		return sqrt(dot(pnt, pnt));
	}
	template <std::floating_point T, std::size_t Dim>
	constexpr point<T, Dim> normalize(const point<T, Dim>& pnt) noexcept
	{
		T n = norm(pnt);
		if (n == 0)
			return 0;
		return pnt / n;
	}

	template <std::floating_point T>
	using point3 = point<T, 3>;

	template <std::floating_point T>
	constexpr point3<T> cross(const point3<T>& lhs, const point3<T>& rhs) noexcept
	{
		return {lhs[1] * rhs[2] - lhs[2] * rhs[1],
				lhs[2] * rhs[0] - lhs[0] * rhs[2],
				lhs[0] * rhs[1] - lhs[1] * rhs[0]};
	}

	template <std::floating_point T, std::size_t Dim>
	using state = std::valarray<point<T, Dim>>; // a "state" is an alias for a valarray of points

	template <std::floating_point T, std::size_t M, std::size_t N>
	constexpr mat<T, M, N> sum_outer(const state<T, M>& lhs, const state<T, N>& rhs) noexcept
	{
		mat<T, M, N> res(0);
		for (std::size_t i = 0; i < lhs.size(); ++i)
			res += outer(lhs[i], rhs[i]);
		return res;
	}

	template <std::floating_point T, std::size_t Dim>
	constexpr T dot(const state<T, Dim>& lhs, const state<T, Dim>& rhs) noexcept
	{
		T res(0);
		for (std::size_t i = 0; i < lhs.size(); ++i)
			res += dot(lhs[i], rhs[i]);
		return res;
	}

	template <std::floating_point T, std::size_t Dim>
	constexpr T norm(const state<T, Dim>& st) noexcept
	{
		using std::sqrt;
		return sqrt(dot(st, st));
	}
	template <std::floating_point T, std::size_t Dim>
	constexpr state<T, Dim> normalize(const state<T, Dim>& st) noexcept
	{
		T n = norm(st);
		if (n == 0)
			return 0;
		return st / n;
	}

	template <std::floating_point T, std::size_t Dim>
	constexpr state<T, Dim> remainder(state<T, Dim> st, T modulo) noexcept
	{
		for (std::size_t i = 0; i < Dim; ++i)
			st[i] = remainder(st[i], modulo);
		return st;
	}

	template <std::size_t N>
	class fixed_list : std::array<unsigned int, N>
	{
		using base = std::array<unsigned int, N>;

		public:

			unsigned int n;

			constexpr fixed_list() noexcept : n(0) {}

			template <typename ... Ts>
			requires (sizeof...(Ts) <= N)
			constexpr fixed_list(Ts ... ts) noexcept
				: base{static_cast<unsigned int>(ts)...}, n(sizeof...(Ts)) {}

			constexpr void push_back(unsigned int i) noexcept
			{
				if (n < N)
					base::operator[](n++) = i;
			}

			using base::operator[];
	};

#define PHYSICS_GEN_MATN_ALIAS(n) \
	template <std::floating_point T> \
	using mat##n = mat<T, n>

	PHYSICS_GEN_MATN_ALIAS(2);
	PHYSICS_GEN_MATN_ALIAS(3);
	PHYSICS_GEN_MATN_ALIAS(4);

	template <std::floating_point T>
	mat3<T> rotation_x(T angle) noexcept
	{
		using std::sin;
		using std::cos;
		T s = sin(angle), c = cos(angle);
		return
		{
			1,	0,	0,
			0,	c,	-s,
			0,	s,	c,
		};
	}
	template <std::floating_point T>
	mat3<T> rotation_y(T angle) noexcept
	{
		using std::sin;
		using std::cos;
		T s = sin(angle), c = cos(angle);
		return
		{
			c,	0,	s,
			0,	1,	0,
			-s,	0,	c,
		};
	}
	template <std::floating_point T>
	mat3<T> rotation_z(T angle) noexcept
	{
		using std::sin;
		using std::cos;
		T s = sin(angle), c = cos(angle);
		return
		{
			c,	-s,	0,
			s,	c,	0,
			0,	0,	1,
		};
	}
	template <std::floating_point T>
	mat3<T> rotation_yaw_pitch_roll(T yaw, T pitch, T roll) noexcept
	{
		return rotation_z(yaw) % rotation_y(pitch) % rotation_x(roll);
	}

} // namespace physics

#endif // PHYSICS_POINT_H
































