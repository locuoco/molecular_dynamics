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

#ifndef PHYSICS_TENSOR_H
#define PHYSICS_TENSOR_H

#include <concepts> // floating_point, integral, convertible_to
#include <valarray>
#include <array>
#include <utility> // make_index_sequence, index_sequence, swap
#include <cmath> // sqrt, remainder, floor, round, trunc
#include <stdexcept> // runtime_error

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
	// create a sequence with `N` elements filled with `scal` values.
	// `R` is the return type, which must accept a brace-initializer list.
	{
		return detail_::make_filled_impl<R, T>(scal, std::make_index_sequence<N>{});
	}

	template <typename T, std::size_t ... Ns>
	requires ((Ns * ...) >= 1)
	struct tensor;

	template <typename T, std::size_t M, std::size_t N = M>
	using mat = tensor<T, M, N>;

	template <typename T, std::size_t Dim>
	using vec = tensor<T, Dim, 1>;

	template <typename T, std::size_t Dim>
	using row_vec = tensor<T, 1, Dim>;

	// tensor, mat, vec traits

	template <typename T>
	inline constexpr bool is_tensor = false;

	template <typename T, std::size_t ... Ns>
	inline constexpr bool is_tensor<tensor<T, Ns...>> = true;

	template <typename T>
	inline constexpr bool is_mat = false;

	template <typename T, std::size_t M, std::size_t N>
	inline constexpr bool is_mat<mat<T, M, N>> = true;

	template <typename T>
	inline constexpr bool is_vec = false;

	template <typename T, std::size_t N>
	inline constexpr bool is_vec<vec<T, N>> = true;

	template <typename T>
	inline constexpr bool is_row_vec = false;

	template <typename T, std::size_t N>
	inline constexpr bool is_row_vec<row_vec<T, N>> = true;

	template <typename T, std::size_t Dim>
	constexpr vec<T, Dim> make_filled_vec(T scal) noexcept
	// create a vector with `Dim` elements filled with `scal` values.
	{
		return make_filled<vec<T, Dim>, Dim>(scal);
	}

	template <typename T, std::size_t Dim>
	constexpr mat<T, Dim> make_identity()
	// create a Dim x Dim identity matrix
	{
		using std::size_t;
		mat<T, Dim> res(0);
		for (size_t i = 0; i < Dim; ++i)
			res(i, i) = 1;
		return res;
	}

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
	// static_multi_index:
	// remap the tensor indices to a plain array (single) index
	// respecting row-major ordering.
	{
		return detail_::static_multi_index_impl<Ns...>(is...); // discard first template parameter
	}

	template <typename T, std::size_t ... Ns>
	requires ((Ns * ...) >= 1)
	struct tensor : std::array<T, (Ns * ...)>
	// a tensor inherits from `std::array` all his member functions (except constructors).
	// `T` is the underlying (scalar) type and `Ns...` is the tensor shape.
	{
		private:
			using base = std::array<T, (Ns * ...)>;

		public:

		constexpr tensor() noexcept = default;
		// trivial (non-user-defined) default constructor
		// since std::array is a standard-layout, tensor is also a standard-layout (since C++17)
		// and so tensor is a POD

		template <typename ... Ts>
		requires (sizeof...(Ts) <= (Ns * ...) && sizeof...(Ts) >= 2)
		constexpr tensor(Ts ... ts) noexcept : base{static_cast<T>(ts)...} {}
		// constructor:
		// the elements of the array are initialized with ts...
		// this overload is choosen only if the number of variadic arguments is 2 or greater.

		constexpr tensor(T t) noexcept : base(make_filled<tensor<T, Ns...>, (Ns * ...)>(t)) {}
		// constructor:
		// the array is filled with `t`'s value.

		template <std::convertible_to<T> S>
		constexpr tensor(const tensor<S, Ns...>& other)
		// template copy constructor
		// it copies another tensor with different scalar type (`S`) but same shape (`Ns...`)
		{
			for (std::size_t i = 0; i < (Ns * ...); ++i)
				(*this)[i] = other[i];
		}

		template <std::convertible_to<std::size_t> ... Ts>
		requires (sizeof...(Ts) == sizeof...(Ns) && sizeof...(Ns) >= 2)
		constexpr T& operator()(Ts ... is)
		// access operator
		// access tensor element with index `is...` respecting row-major ordering.
		// return a reference to the element
		{
			return (*this)[static_multi_index<Ns...>(is...)];
		}

		template <std::convertible_to<std::size_t> ... Ts>
		requires (sizeof...(Ts) == sizeof...(Ns) && sizeof...(Ns) >= 2)
		constexpr const T& operator()(Ts ... is) const
		// access operator
		// access tensor element with index `is...` respecting row-major ordering.
		// return a reference to const of the element.
		{
			return (*this)[static_multi_index<Ns...>(is...)];
		}

		constexpr bool all(T t) const
		// check if all elements are equal to `t`
		{
			for (auto& e : *this)
				if (e != t)
					return false;
			return true;
		}

		constexpr bool any(T t) const
		// check if any element is equal to `t`
		{
			for (auto& e : *this)
				if (e == t)
					return true;
			return false;
		}

		constexpr T sum() const
		// sum all elements of the tensor.
		{
			T accum = 0;
			for (auto& e : *this)
				accum += e;
			return accum;
		}

		constexpr const tensor& operator+() const noexcept
		// unary plus operator
		{
			return *this;
		}
		constexpr tensor operator-() const
		// negate operator
		{
			tensor res;
			for (std::size_t i = 0; i < (Ns * ...); ++i)
				res[i] = -(*this)[i];
			return res;
		}
		constexpr tensor& operator+=(const tensor& other)
		// add-assignment operator with an`other` tensor
		{
			for (std::size_t i = 0; i < (Ns * ...); ++i)
				(*this)[i] += other[i];
			return *this;
		}
		constexpr tensor& operator-=(const tensor& other)
		// subtract-assignment operator with an`other` tensor
		{
			for (std::size_t i = 0; i < (Ns * ...); ++i)
				(*this)[i] -= other[i];
			return *this;
		}
		constexpr tensor& operator*=(const tensor& other)
		// (element-wise) multiply-assignment operator with an`other` tensor
		{
			for (std::size_t i = 0; i < (Ns * ...); ++i)
				(*this)[i] *= other[i];
			return *this;
		}
		constexpr tensor& operator/=(const tensor& other)
		// (element-wise) divide-assignment operator with an`other` tensor
		{
			for (std::size_t i = 0; i < (Ns * ...); ++i)
				(*this)[i] /= other[i];
			return *this;
		}
		constexpr tensor& operator+=(T scal)
		// add-assignment operator with a `scal`ar
		{
			for (auto& e : *this)
				e += scal;
			return *this;
		}
		constexpr tensor& operator-=(T scal)
		// subtract-assignment operator with a `scal`ar
		{
			for (auto& e : *this)
				e -= scal;
			return *this;
		}
		constexpr tensor& operator*=(T scal)
		// multiply-assignment operator with a `scal`ar
		{
			for (auto& e : *this)
				e *= scal;
			return *this;
		}
		constexpr tensor& operator/=(T scal)
		// divide-assignment operator with a `scal`ar
		{
			for (auto& e : *this)
				e /= scal;
			return *this;
		}
	};

	template <typename T, std::size_t ... Ns>
	constexpr tensor<T, Ns...> operator+(tensor<T, Ns...> lhs, const tensor<T, Ns...>& rhs)
	// add operator between two `tensor`s
	{
		return lhs += rhs;
	}
	template <typename T, std::size_t ... Ns>
	constexpr tensor<T, Ns...> operator-(tensor<T, Ns...> lhs, const tensor<T, Ns...>& rhs)
	// subtract operator between two `tensor`s
	{
		return lhs -= rhs;
	}
	template <typename T, std::size_t ... Ns>
	constexpr tensor<T, Ns...> operator*(tensor<T, Ns...> lhs, const tensor<T, Ns...>& rhs)
	// (element-wise) multiply operator between two `tensor`s
	{
		return lhs *= rhs;
	}
	template <typename T, std::size_t ... Ns>
	constexpr tensor<T, Ns...> operator/(tensor<T, Ns...> lhs, const tensor<T, Ns...>& rhs)
	// (element-wise) divide operator between two `tensor`s
	{
		return lhs /= rhs;
	}
	template <typename T, std::size_t ... Ns>
	constexpr tensor<T, Ns...> operator+(tensor<T, Ns...> lhs, typename tensor<T, Ns...>::value_type scal)
	// add operator between a tensor and a scalar
	{
		return lhs += scal;
	}
	template <typename T, std::size_t ... Ns>
	constexpr tensor<T, Ns...> operator+(typename tensor<T, Ns...>::value_type scal, tensor<T, Ns...> rhs)
	// add operator between a scalar and a tensor
	{
		return rhs + scal;
	}
	template <typename T, std::size_t ... Ns>
	constexpr tensor<T, Ns...> operator-(tensor<T, Ns...> lhs, typename tensor<T, Ns...>::value_type scal)
	// subtract operator between a tensor and a scalar
	{
		return lhs -= scal;
	}
	template <typename T, std::size_t ... Ns>
	constexpr tensor<T, Ns...> operator-(typename tensor<T, Ns...>::value_type scal, tensor<T, Ns...> rhs)
	// subtract operator between a scalar and a tensor
	{
		return rhs - scal;
	}
	template <typename T, std::size_t ... Ns>
	constexpr tensor<T, Ns...> operator*(tensor<T, Ns...> lhs, typename tensor<T, Ns...>::value_type scal)
	// multiply operator between a tensor and a scalar
	{
		return lhs *= scal;
	}
	template <typename T, std::size_t ... Ns>
	constexpr tensor<T, Ns...> operator*(typename tensor<T, Ns...>::value_type scal, tensor<T, Ns...> rhs)
	// multiply operator between a scalar and a tensor
	{
		return rhs * scal;
	}
	template <typename T, std::size_t ... Ns>
	constexpr tensor<T, Ns...> operator/(tensor<T, Ns...> lhs, typename tensor<T, Ns...>::value_type scal)
	// divide operator between a tensor and a scalar
	{
		return lhs /= scal;
	}
	template <typename T, std::size_t ... Ns>
	constexpr tensor<T, Ns...> operator/(typename tensor<T, Ns...>::value_type scal, tensor<T, Ns...> rhs)
	// (element-wise) divide operator between a scalar and a tensor
	{
		return rhs / scal;
	}

	template <typename T, std::size_t M, std::size_t N, std::size_t K>
	constexpr mat<T, M, K> matmul(const mat<T, M, N>& m1, const mat<T, N, K>& m2)
	// matrix multiplication between `m1` and `m2`.
	// It can be used also for matrix-vector (K=1), vector-matrix multiplication (M=1)
	// and dot product between a `row_vec` and a `vec` (M=K=1).
	{
		using std::size_t;
		mat<T, M, K> res(0);
		for (size_t i = 0; i < M; ++i)
			for (size_t j = 0; j < N; ++j)
				for (size_t k = 0; k < K; ++k)
					res(i, k) += m1(i, j) * m2(j, k);
		return res;
	}

	template <std::floating_point T, std::size_t M, std::size_t N, std::size_t K>
	constexpr mat<T, M, K> operator%(const mat<T, M, N>& lhs, const mat<T, N, K>& rhs)
	// matrix multiplication operator, defined only for floating point scalar types.
	// Same as matmul(lhs, rhs).
	{
		return matmul(lhs, rhs);
	}

	template <std::floating_point T, std::size_t ... Ns>
	constexpr tensor<T, Ns...> remainder(tensor<T, Ns...> t, typename tensor<T, Ns...>::value_type modulo)
	// remainder function applied to tensor `t` element-wise
	{
		using std::round;
		for (auto& e : t)
			e = e - round(e/modulo)*modulo;
		return t;
	}
	template <typename T, std::size_t ... Ns>
	constexpr tensor<T, Ns...> mod(tensor<T, Ns...> t, typename tensor<T, Ns...>::value_type modulo)
	// mod function applied to tensor `t` element-wise
	{
		using math::mod;
		for (auto& e : t)
			e = mod(e, modulo);
		return t;
	}
	template <std::floating_point T, std::size_t ... Ns>
	constexpr tensor<T, Ns...> floor(tensor<T, Ns...> t)
	// floor function applied to tensor `t` element-wise
	{
		using std::floor;
		for (auto& e : t)
			e = floor(e);
		return t;
	}
	template <std::floating_point T, std::size_t ... Ns>
	constexpr tensor<T, Ns...> round(tensor<T, Ns...> t)
	// round function applied to tensor `t` element-wise
	{
		using std::round;
		for (auto& e : t)
			e = round(e);
		return t;
	}
	template <std::floating_point T, std::size_t ... Ns>
	constexpr tensor<T, Ns...> trunc(tensor<T, Ns...> t)
	// trunc function applied to tensor `t` element-wise
	{
		using std::trunc;
		for (auto& e : t)
			e = trunc(e);
		return t;
	}

	template <typename T, std::size_t M, std::size_t N>
	constexpr mat<T, M, N> outer(const vec<T, M>& u, const vec<T, N>& v)
	// outer product between two vecs `u` and `v`
	// return a `mat` with dimensions corresponding to the two vectors.
	{
		using std::size_t;
		mat<T, M, N> res;
		for (size_t i = 0; i < M; ++i)
			for (size_t j = 0; j < N; ++j)
				res(i, j) = u[i] * v[j];
		return res;
	}

	template <typename T, std::size_t Dim>
	constexpr T trace(mat<T, Dim> A)
	// trace of a matrix `A`
	{
		T res(0);
		for (std::size_t i = 0; i < Dim; ++i)
			res += A(i, i);
		return res;
	}

	template <typename T, std::size_t Dim>
	constexpr T dot(const vec<T, Dim>& u, const vec<T, Dim>& v)
	// scalar product between two vectors `u` and `v`
	// return a scalar of type `T`
	{
		T res(0);
		for (std::size_t i = 0; i < Dim; ++i)
			res += u[i] * v[i];
		return res;
	}

	template <std::floating_point T, std::size_t Dim>
	T norm(const vec<T, Dim>& v)
	// norm of a vector
	{
		using std::sqrt;
		return sqrt(dot(v, v));
	}
	template <std::floating_point T, std::size_t Dim>
	vec<T, Dim> normalize(const vec<T, Dim>& v)
	// return the normalized unit vector of `v`
	{
		T n = norm(v);
		if (n == 0)
			return 0;
		return v / n;
	}

	// generate vec aliases vec2<T>, vec3<T>, ... using a macro

#define PHYSICS_GEN_VECN_ALIAS(N) \
	template <typename T> \
	using vec##N = vec<T, N>

	PHYSICS_GEN_VECN_ALIAS(2);
	PHYSICS_GEN_VECN_ALIAS(3);
	PHYSICS_GEN_VECN_ALIAS(4);

	// generate vec aliases vec2f, vec3d, ... using a macro

#define PHYSICS_GEN_VECNT_ALIAS(TYPE, SUFFIX) \
	using vec2##SUFFIX = vec2<TYPE>; \
	using vec3##SUFFIX = vec3<TYPE>; \
	using vec4##SUFFIX = vec4<TYPE>

	PHYSICS_GEN_VECNT_ALIAS(int, i);
	PHYSICS_GEN_VECNT_ALIAS(float, f);
	PHYSICS_GEN_VECNT_ALIAS(double, d);
	PHYSICS_GEN_VECNT_ALIAS(long double, ld);

	// generate mat aliases mat2<T>, mat3<T>, ... using a macro

#define PHYSICS_GEN_MATN_ALIAS(n) \
	template <typename T> \
	using mat##n = mat<T, n>

	PHYSICS_GEN_MATN_ALIAS(2);
	PHYSICS_GEN_MATN_ALIAS(3);
	PHYSICS_GEN_MATN_ALIAS(4);

	// generate mat aliases mat2f, mat3d, ... using a macro

#define PHYSICS_GEN_MATNT_ALIAS(TYPE, SUFFIX) \
	using mat2##SUFFIX = mat2<TYPE>; \
	using mat3##SUFFIX = mat3<TYPE>; \
	using mat4##SUFFIX = mat4<TYPE>

	PHYSICS_GEN_MATNT_ALIAS(int, i);
	PHYSICS_GEN_MATNT_ALIAS(float, f);
	PHYSICS_GEN_MATNT_ALIAS(double, d);
	PHYSICS_GEN_MATNT_ALIAS(long double, ld);

	template <typename T>
	constexpr vec3<T> cross(const vec3<T>& u, const vec3<T>& v)
	// cross product between two vectors `u` and `v`
	{
		return
		{
			u[1] * v[2] - u[2] * v[1],
			u[2] * v[0] - u[0] * v[2],
			u[0] * v[1] - u[1] * v[0],
		};
	}

	template <typename T, std::size_t Dim>
	constexpr void inplace_transpose(mat<T, Dim>& m)
	// transpose matrix `m` in-place.
	// Only for square Dim x Dim matrices.
	{
		using std::size_t;
		using std::swap;
		for (size_t i = 0; i < Dim; ++i)
			for (size_t j = 0; j < i; ++j)
				swap(m(i, j), m(j, i));
	}

	template <typename T, std::size_t Dim>
	constexpr mat<T, Dim> transpose(mat<T, Dim> m)
	// transpose matrix `m` out-of-place.
	// Return the transposition without modifying the original matrix.
	// Only for square Dim x Dim matrices.
	{
		inplace_transpose(m);
		return m;
	}

	template <typename T, std::size_t M, std::size_t N>
	constexpr bool is_transpose(const mat<T, M, N>& m1, const mat<T, N, M>& m2)
	// check if `m1` is the transpose of `m2` (or, equivalently, viceversa)
	// assuming exact arithmetic (may not be good for certain tests).
	{
		using std::size_t;
		for (size_t i = 0; i < M; ++i)
			for (size_t j = 0; j < N; ++j)
				if (m1(i, j) != m2(j, i))
					return false;
		return true;
	}

	template <typename T, std::size_t Dim>
	using state = std::valarray<vec<T, Dim>>; // a "state" is an alias for a valarray of vecs
	// a valarray defines +,-,*,/ operators and some math functions

	template <typename T, std::size_t M, std::size_t N>
	constexpr mat<T, M, N> sum_outer(const state<T, M>& s1, const state<T, N>& s2)
	// calculate outer product between states s1 and s2 elements and sum them.
	// If s1.size() != s2.size(), the behaviour is undefined.
	// Return a M x N matrix with the sum of the outer products.
	{
		mat<T, M, N> res(0);
		for (std::size_t i = 0; i < s1.size(); ++i)
			res += outer(s1[i], s2[i]);
		return res;
	}

	template <typename T, std::size_t Dim>
	constexpr T dot(const state<T, Dim>& s1, const state<T, Dim>& s2)
	// calculate dot product between states s1 and s2.
	// If s1.size() != s2.size(), the behaviour is undefined.
	// Return a scalar of type `T`.
	{
		T res(0);
		for (std::size_t i = 0; i < s1.size(); ++i)
			res += dot(s1[i], s2[i]);
		return res;
	}

	template <std::floating_point T, std::size_t Dim>
	T norm(const state<T, Dim>& st)
	// calculate norm of state `st`
	{
		using std::sqrt;
		return sqrt(dot(st, st));
	}
	template <std::floating_point T, std::size_t Dim>
	state<T, Dim> normalize(const state<T, Dim>& st)
	// return the normalized unit vector of `st`
	{
		T n = norm(st);
		if (n == 0)
			return 0;
		return st / n;
	}

	// element-wise functions

	template <std::floating_point T, std::size_t Dim>
	constexpr state<T, Dim> remainder(state<T, Dim> st, T modulo)
	{
		for (auto& v : st)
			v = remainder(v, modulo);
		return st;
	}
	template <std::floating_point T, std::size_t Dim>
	constexpr state<T, Dim> mod(state<T, Dim> st, T modulo)
	{
		for (auto& v : st)
			v = mod(v, modulo);
		return st;
	}
	template <std::floating_point T, std::size_t Dim>
	constexpr state<T, Dim> floor(state<T, Dim> st)
	{
		for (auto& v : st)
			v = floor(v);
		return st;
	}
	template <std::floating_point T, std::size_t Dim>
	constexpr state<T, Dim> round(state<T, Dim> st)
	{
		for (auto& v : st)
			v = round(v);
		return st;
	}
	template <std::floating_point T, std::size_t Dim>
	constexpr state<T, Dim> trunc(state<T, Dim> st)
	{
		for (auto& v : st)
			v = trunc(v);
		return st;
	}

	template <std::size_t N>
	class fixed_list : std::array<unsigned, N>
	// a wrapper of a `std::array<unsigned, N>` with a push_back method and a counter `n`.
	// used in `molecule.hpp`
	{
		using base = std::array<unsigned, N>;

		public:

			unsigned n;

			constexpr fixed_list() noexcept : n(0) {}
			// default constructor:
			// set the counter to 0

			template <typename ... Ts>
			requires (sizeof...(Ts) <= N)
			constexpr fixed_list(Ts ... ts) noexcept
			// constructor:
			// put `ts...` elements inside the list and set the counter accordingly
				: base{static_cast<unsigned>(ts)...}, n(sizeof...(Ts))
			{}

			constexpr void push_back(unsigned i)
			// insert element `i` and increment the counter.
			// If n >= N, throw a `std::runtime_error` exception.
			{
				if (n >= N)
					throw std::runtime_error("Cannot put other elements inside fixed_list!");
				(*this)[n++] = i;
			}

			constexpr std::size_t size() const noexcept
			{
				return n;
			}

			// past-the-end iterators
			constexpr auto end() noexcept
			{
				return base::begin()+n;
			}
			constexpr auto end() const noexcept
			{
				return base::begin()+n;
			}
			constexpr auto cend() const noexcept
			{
				return base::cbegin()+n;
			}
			// reversed begin iterators
			constexpr auto rbegin() noexcept
			{
				return base::rend()-n;
			}
			constexpr auto rbegin() const noexcept
			{
				return base::rend()-n;
			}
			constexpr auto crbegin() const noexcept
			{
				return base::crend()-n;
			}

			using base::operator[];
			using base::begin;
			using base::cbegin;
			using base::rend;
			using base::crend;
	};

	// rotation matrices are used in graphics

	template <std::floating_point T>
	mat3<T> rotation_x(T angle) noexcept
	// rotation matrix along x-axis
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
	// rotation matrix along y-axis
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
	// rotation matrix along z-axis
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
	// rotation matrix using yaw, pitch, roll angles
	{
		return rotation_z(yaw) % rotation_y(pitch) % rotation_x(roll);
	}

	template <std::floating_point T>
	constexpr mat4<T> orthographic_projection(T left, T right, T bottom, T top) noexcept
	// orthographic projection matrix (row-major order)
	// `left`, `right`, `bottom` and `top` are the positions of the clipping planes
	// (near and far planes, along z-axis, are set to -inf and +inf)
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
	// orthographic projection matrix for column-major order
	// `left`, `right`, `bottom` and `top` are the positions of the clipping planes
	{
		mat4<T> res(orthographic_projection(left, right, bottom, top));
		inplace_transpose(res);
		return res;
	}

	template <std::floating_point T>
	mat4<T> look_at(const vec3<T>& eye, const vec3<T>& center, const vec3<T>& up)
	// look-at matrix (see gluLookAt Khronos reference page)
	// `eye` is the position of the camera, `center` the position it points to and
	// `up` is the up-vector
	{
		vec3<T> f(normalize(center - eye));
		vec3<T> s(normalize(cross(f, up)));
		vec3<T> u(cross(s, f));

		return
		{
			s[0], s[1], s[2], -dot(s, eye),
			u[0], u[1], u[2], -dot(u, eye),
			-f[0], -f[1], -f[2], dot(f, eye),
			0, 0, 0, 1,
		};
	}

	template <std::floating_point T>
	mat4<T> look_at_col(const vec3<T>& eye, const vec3<T>& center, const vec3<T>& up)
	// look-at matrix (see gluLookAt Khronos reference page) for column-major order
	// `eye` is the position of the camera, `center` the position it points to and
	// `up` is the up-vector
	{
		mat4<T> res(look_at(eye, center, up));
		inplace_transpose(res);
		return res;
	}

	template <typename Vec>
	requires (std::assignable_from<std::valarray<typename Vec::value_type>&, Vec> && is_tensor<typename Vec::value_type>)
	auto rms(const Vec& v)
	// calculate root mean square of a valarray of tensors `v`, using the formula:
	// rms(v) = sqrt( sum_i ||v_i||^2 / n)
	// where n is the number of tensors (i.e. v.size()).
	{
		using std::sqrt;
		auto sqr = [](physics::vec3d x) { return x*x; };
		return sqrt(v.apply(sqr).sum().sum() / v.size());
	}

} // namespace physics

#endif // PHYSICS_TENSOR_H
































