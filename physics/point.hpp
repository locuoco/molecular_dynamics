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

#include <concepts> // floating_point, integral, convertible_to
#include <valarray>
#include <vector>
#include <array>
#include <utility> // make_index_sequence, index_sequence, swap
#include <cmath> // sqrt, remainder, floor, round, trunc

#include "../math/helper.hpp" // mod

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

	template <typename T, std::size_t ... Ns>
	requires ((Ns * ...) >= 1)
	struct tens; // tensor

	template <typename T, std::size_t M, std::size_t N = M>
	using mat = tens<T, M, N>;

	template <typename T, std::size_t Dim>
	using point = tens<T, Dim, 1>;

	template <typename T, std::size_t Dim>
	constexpr point<T, Dim> make_filled_point(T scal) noexcept
	{
		return make_filled<point<T, Dim>, Dim>(scal);
	}

	template <typename T, std::size_t Dim>
	constexpr mat<T, Dim> make_identity_matrix()
	{
		using std::size_t;
		mat<T, Dim> res(0);
		for (size_t i = 0; i < Dim; ++i)
			res(i, i) = 1;
		return res;
	}

	// static_multi_index
	// multi-index to single index conversion
	namespace detail_
	{
		template <std::size_t N = 1, std::size_t ... Ns, typename ... Ts>
		constexpr std::size_t static_multi_index_impl(std::size_t i = 0, std::size_t j = 0, Ts ... is) noexcept
		{
			return static_multi_index_impl<Ns...>(i*N + j, is...);
		}
		template <>
		constexpr std::size_t static_multi_index_impl<>(std::size_t i, std::size_t) noexcept
		{
			return i;
		}
	}
	template <std::size_t N = 1, std::size_t ... Ns, typename ... Ts>
	constexpr std::size_t static_multi_index(Ts ... is) noexcept
	{
		return detail_::static_multi_index_impl<Ns...>(is...); // discard first template parameter
	}

	template <typename T, std::size_t ... Ns>
	requires ((Ns * ...) >= 1)
	struct tens : std::array<T, (Ns * ...)>
	{
		private:
			using base = std::array<T, (Ns * ...)>;

		public:

		constexpr tens() noexcept = default;

		template <typename ... Ts>
		requires (sizeof...(Ts) <= (Ns * ...) && sizeof...(Ts) >= 2 && ((std::floating_point<Ts> || std::integral<Ts>) && ...))
		constexpr tens(Ts ... ts) noexcept : base{static_cast<T>(ts)...} {}

		constexpr tens(T t) noexcept : base(make_filled<tens<T, Ns...>, (Ns * ...)>(t)) {}

		template <std::convertible_to<T> S>
		constexpr tens(const tens<S, Ns...>& other)
		{
			for (std::size_t i = 0; i < (Ns * ...); ++i)
				base::operator[](i) = other[i];
		}

		template <std::convertible_to<std::size_t> ... Ts>
		requires (sizeof...(Ts) == sizeof...(Ns) && sizeof...(Ns) >= 2)
		constexpr T& operator()(Ts ... is)
		{
			return base::operator[](static_multi_index<Ns...>(is...));
		}

		template <std::convertible_to<std::size_t> ... Ts>
		requires (sizeof...(Ts) == sizeof...(Ns) && sizeof...(Ns) >= 2)
		constexpr const T& operator()(Ts ... is) const
		{
			return base::operator[](static_multi_index<Ns...>(is...));
		}

		constexpr const tens& operator+() const noexcept
		{
			return *this;
		}
		constexpr tens operator-() const
		{
			tens res;
			for (std::size_t i = 0; i < (Ns * ...); ++i)
				res[i] = -base::operator[](i);
			return res;
		}
		constexpr tens& operator+=(const tens& other)
		{
			for (std::size_t i = 0; i < (Ns * ...); ++i)
				base::operator[](i) += other[i];
			return *this;
		}
		constexpr tens& operator-=(const tens& other)
		{
			for (std::size_t i = 0; i < (Ns * ...); ++i)
				base::operator[](i) -= other[i];
			return *this;
		}
		constexpr tens& operator*=(const tens& other)
		{
			for (std::size_t i = 0; i < (Ns * ...); ++i)
				base::operator[](i) *= other[i];
			return *this;
		}
		constexpr tens& operator/=(const tens& other)
		{
			for (std::size_t i = 0; i < (Ns * ...); ++i)
				base::operator[](i) /= other[i];
			return *this;
		}
		constexpr tens& operator+=(T scal)
		{
			for (std::size_t i = 0; i < (Ns * ...); ++i)
				base::operator[](i) += scal;
			return *this;
		}
		constexpr tens& operator-=(T scal)
		{
			for (std::size_t i = 0; i < (Ns * ...); ++i)
				base::operator[](i) -= scal;
			return *this;
		}
		constexpr tens& operator*=(T scal)
		{
			for (std::size_t i = 0; i < (Ns * ...); ++i)
				base::operator[](i) *= scal;
			return *this;
		}
		constexpr tens& operator/=(T scal)
		{
			for (std::size_t i = 0; i < (Ns * ...); ++i)
				base::operator[](i) /= scal;
			return *this;
		}
	};

	template <typename T, std::size_t ... Ns>
	constexpr tens<T, Ns...> operator+(tens<T, Ns...> lhs, const tens<T, Ns...>& rhs)
	{
		return lhs += rhs;
	}
	template <typename T, std::size_t ... Ns>
	constexpr tens<T, Ns...> operator-(tens<T, Ns...> lhs, const tens<T, Ns...>& rhs)
	{
		return lhs -= rhs;
	}
	template <typename T, std::size_t ... Ns>
	constexpr tens<T, Ns...> operator*(tens<T, Ns...> lhs, const tens<T, Ns...>& rhs)
	{
		return lhs *= rhs;
	}
	template <typename T, std::size_t ... Ns>
	constexpr tens<T, Ns...> operator/(tens<T, Ns...> lhs, const tens<T, Ns...>& rhs)
	{
		return lhs /= rhs;
	}
	template <typename T, std::size_t ... Ns>
	constexpr tens<T, Ns...> operator+(tens<T, Ns...> lhs, typename tens<T, Ns...>::value_type scal)
	{
		return lhs += scal;
	}
	template <typename T, std::size_t ... Ns>
	constexpr tens<T, Ns...> operator+(typename tens<T, Ns...>::value_type scal, tens<T, Ns...> rhs)
	{
		return rhs + scal;
	}
	template <typename T, std::size_t ... Ns>
	constexpr tens<T, Ns...> operator-(tens<T, Ns...> lhs, typename tens<T, Ns...>::value_type scal)
	{
		return lhs -= scal;
	}
	template <typename T, std::size_t ... Ns>
	constexpr tens<T, Ns...> operator-(typename tens<T, Ns...>::value_type scal, tens<T, Ns...> rhs)
	{
		return rhs - scal;
	}
	template <typename T, std::size_t ... Ns>
	constexpr tens<T, Ns...> operator*(tens<T, Ns...> lhs, typename tens<T, Ns...>::value_type scal)
	{
		return lhs *= scal;
	}
	template <typename T, std::size_t ... Ns>
	constexpr tens<T, Ns...> operator*(typename tens<T, Ns...>::value_type scal, tens<T, Ns...> rhs)
	{
		return rhs * scal;
	}
	template <typename T, std::size_t ... Ns>
	constexpr tens<T, Ns...> operator/(tens<T, Ns...> lhs, typename tens<T, Ns...>::value_type scal)
	{
		return lhs /= scal;
	}
	template <typename T, std::size_t ... Ns>
	constexpr tens<T, Ns...> operator/(typename tens<T, Ns...>::value_type scal, tens<T, Ns...> rhs)
	{
		return rhs / scal;
	}

	template <std::floating_point T, std::size_t M, std::size_t N, std::size_t K>
	constexpr mat<T, M, K> operator%(const mat<T, M, N>& lhs, const mat<T, N, K>& rhs)
	{
		using std::size_t;
		mat<T, M, K> res(0);
		for (size_t i = 0; i < M; ++i)
			for (size_t j = 0; j < N; ++j)
				for (size_t k = 0; k < K; ++k)
					res(i, k) += lhs(i, j) * rhs(j, k);
		return res;
	}

	template <std::floating_point T, std::size_t ... Ns>
	constexpr tens<T, Ns...> remainder(tens<T, Ns...> t, typename tens<T, Ns...>::value_type modulo)
	{
		using std::round;
		for (std::size_t i = 0; i < (Ns * ...); ++i)
			t[i] = t[i] - round(t[i]/modulo)*modulo;
		return t;
	}
	template <std::floating_point T, std::size_t ... Ns>
	constexpr tens<T, Ns...> mod(tens<T, Ns...> t, typename tens<T, Ns...>::value_type modulo)
	{
		using math::mod;
		for (std::size_t i = 0; i < (Ns * ...); ++i)
			t[i] = mod(t[i], modulo);
		return t;
	}
	template <std::floating_point T, std::size_t ... Ns>
	constexpr tens<T, Ns...> floor(tens<T, Ns...> t)
	{
		using std::floor;
		for (std::size_t i = 0; i < (Ns * ...); ++i)
			t[i] = floor(t[i]);
		return t;
	}
	template <std::floating_point T, std::size_t ... Ns>
	constexpr tens<T, Ns...> round(tens<T, Ns...> t)
	{
		using std::round;
		for (std::size_t i = 0; i < (Ns * ...); ++i)
			t[i] = round(t[i]);
		return t;
	}
	template <std::floating_point T, std::size_t ... Ns>
	constexpr tens<T, Ns...> trunc(tens<T, Ns...> t)
	{
		using std::trunc;
		for (std::size_t i = 0; i < (Ns * ...); ++i)
			t[i] = trunc(t[i]);
		return t;
	}

	template <typename T, std::size_t M, std::size_t N>
	constexpr mat<T, M, N> outer(const point<T, M>& lhs, const point<T, N>& rhs)
	{
		using std::size_t;
		mat<T, M, N> res;
		for (size_t i = 0; i < M; ++i)
			for (size_t j = 0; j < N; ++j)
				res(i, j) = lhs[i] * rhs[j];
		return res;
	}

	template <typename T, std::size_t Dim>
	constexpr T trace(mat<T, Dim> A)
	{
		T res(0);
		for (std::size_t i = 0; i < Dim; ++i)
			res += A(i, i);
		return res;
	}

	template <typename T, std::size_t Dim>
	constexpr T dot(const point<T, Dim>& lhs, const point<T, Dim>& rhs)
	{
		T res(0);
		for (std::size_t i = 0; i < Dim; ++i)
			res += lhs[i] * rhs[i];
		return res;
	}

	template <std::floating_point T, std::size_t Dim>
	constexpr T norm(const point<T, Dim>& pnt)
	{
		using std::sqrt;
		return sqrt(dot(pnt, pnt));
	}
	template <std::floating_point T, std::size_t Dim>
	constexpr point<T, Dim> normalize(const point<T, Dim>& pnt)
	{
		T n = norm(pnt);
		if (n == 0)
			return 0;
		return pnt / n;
	}

#define PHYSICS_GEN_POINTN_ALIAS(N) \
	template <typename T> \
	using point##N = point<T, N>

	PHYSICS_GEN_POINTN_ALIAS(2);
	PHYSICS_GEN_POINTN_ALIAS(3);
	PHYSICS_GEN_POINTN_ALIAS(4);

#define PHYSICS_GEN_POINTNT_ALIAS(TYPE, SUFFIX) \
	using point2##SUFFIX = point2<TYPE>; \
	using point3##SUFFIX = point3<TYPE>; \
	using point4##SUFFIX = point4<TYPE>

	PHYSICS_GEN_POINTNT_ALIAS(int, i);
	PHYSICS_GEN_POINTNT_ALIAS(float, f);
	PHYSICS_GEN_POINTNT_ALIAS(double, d);
	PHYSICS_GEN_POINTNT_ALIAS(long double, ld);

#define PHYSICS_GEN_MATN_ALIAS(n) \
	template <typename T> \
	using mat##n = mat<T, n>

	PHYSICS_GEN_MATN_ALIAS(2);
	PHYSICS_GEN_MATN_ALIAS(3);
	PHYSICS_GEN_MATN_ALIAS(4);

#define PHYSICS_GEN_MATNT_ALIAS(TYPE, SUFFIX) \
	using mat2##SUFFIX = mat2<TYPE>; \
	using mat3##SUFFIX = mat3<TYPE>; \
	using mat4##SUFFIX = mat4<TYPE>

	PHYSICS_GEN_MATNT_ALIAS(int, i);
	PHYSICS_GEN_MATNT_ALIAS(float, f);
	PHYSICS_GEN_MATNT_ALIAS(double, d);
	PHYSICS_GEN_MATNT_ALIAS(long double, ld);

	template <typename T>
	constexpr point3<T> cross(const point3<T>& lhs, const point3<T>& rhs)
	{
		return {lhs[1] * rhs[2] - lhs[2] * rhs[1],
				lhs[2] * rhs[0] - lhs[0] * rhs[2],
				lhs[0] * rhs[1] - lhs[1] * rhs[0]};
	}

	template <typename T, std::size_t Dim>
	constexpr void transpose(mat<T, Dim>& A)
	{
		using std::size_t;
		using std::swap;
		for (size_t i = 0; i < Dim; ++i)
			for (size_t j = 0; j < i; ++j)
				swap(A(i, j), A(j, i));
	}

	template <typename T, std::size_t Dim>
	using state = std::valarray<point<T, Dim>>; // a "state" is an alias for a valarray of points

	template <typename T, std::size_t M, std::size_t N>
	constexpr mat<T, M, N> sum_outer(const state<T, M>& lhs, const state<T, N>& rhs)
	{
		mat<T, M, N> res(0);
		for (std::size_t i = 0; i < lhs.size(); ++i)
			res += outer(lhs[i], rhs[i]);
		return res;
	}

	template <typename T, std::size_t Dim>
	constexpr T dot(const state<T, Dim>& lhs, const state<T, Dim>& rhs)
	{
		T res(0);
		for (std::size_t i = 0; i < lhs.size(); ++i)
			res += dot(lhs[i], rhs[i]);
		return res;
	}

	template <std::floating_point T, std::size_t Dim>
	T norm(const state<T, Dim>& st)
	{
		using std::sqrt;
		return sqrt(dot(st, st));
	}
	template <std::floating_point T, std::size_t Dim>
	state<T, Dim> normalize(const state<T, Dim>& st)
	{
		T n = norm(st);
		if (n == 0)
			return 0;
		return st / n;
	}

	template <std::floating_point T, std::size_t Dim>
	constexpr state<T, Dim> remainder(state<T, Dim> st, T modulo)
	{
		for (std::size_t i = 0; i < Dim; ++i)
			st[i] = remainder(st[i], modulo);
		return st;
	}
	template <std::floating_point T, std::size_t Dim>
	constexpr state<T, Dim> mod(state<T, Dim> st, T modulo)
	{
		for (std::size_t i = 0; i < Dim; ++i)
			st[i] = mod(st[i], modulo);
		return st;
	}
	template <std::floating_point T, std::size_t Dim>
	constexpr state<T, Dim> floor(state<T, Dim> st)
	{
		for (std::size_t i = 0; i < Dim; ++i)
			st[i] = floor(st[i]);
		return st;
	}
	template <std::floating_point T, std::size_t Dim>
	constexpr state<T, Dim> round(state<T, Dim> st)
	{
		for (std::size_t i = 0; i < Dim; ++i)
			st[i] = round(st[i]);
		return st;
	}
	template <std::floating_point T, std::size_t Dim>
	constexpr state<T, Dim> trunc(state<T, Dim> st)
	{
		for (std::size_t i = 0; i < Dim; ++i)
			st[i] = trunc(st[i]);
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

			constexpr void push_back(unsigned int i)
			{
				if (n < N)
					base::operator[](n++) = i;
			}

			using base::operator[];
	};

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
	mat3<T> rotation_yaw_pitch_roll(T yaw, T pitch, T roll)
	{
		return rotation_z(yaw) % rotation_y(pitch) % rotation_x(roll);
	}

	template <std::floating_point T>
	constexpr mat4<T> orthographic_projection(T left, T right, T bottom, T top) noexcept
	{
		return
		{
			2/(right-left), 0, 0, (right+left)/(left-right),
			0, 2/(top-bottom), 0, (top+bottom)/(bottom-top),
			0, 0, -1, 0,
			0, 0, 0, 1,
		};
	}

	template <std::floating_point T>
	constexpr mat4<T> orthographic_projection_col(T left, T right, T bottom, T top)
	{
		mat4<T> res(orthographic_projection(left, right, bottom, top));
		transpose(res);
		return res;
	}

	template <std::floating_point T>
	mat4<T> look_at(const point3<T>& eye, const point3<T>& center, const point3<T>& up)
	{
		point3<T> f(normalize(center - eye));
		point3<T> s(normalize(cross(f, up)));
		point3<T> u(cross(s, f));

		return
		{
			s[0], s[1], s[2], -dot(s, eye),
			u[0], u[1], u[2], -dot(u, eye),
			-f[0], -f[1], -f[2], dot(f, eye),
			0, 0, 0, 1,
		};
	}

	template <std::floating_point T>
	mat4<T> look_at_col(const point3<T>& eye, const point3<T>& center, const point3<T>& up)
	{
		mat4<T> res(look_at(eye, center, up));
		transpose(res);
		return res;
	}

} // namespace physics

#endif // PHYSICS_POINT_H
































