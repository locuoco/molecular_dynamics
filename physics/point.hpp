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

	template <std::floating_point T, std::size_t Dim>
	requires (Dim >= 1)
	struct point;

	template <std::floating_point T, std::size_t Dim>
	constexpr point<T, Dim> make_point(T scal) noexcept
	{
		return make_filled<point<T, Dim>, Dim>(scal);
	}

	template <std::floating_point T, std::size_t Dim>
	requires (Dim >= 1)
	struct point : std::array<T, Dim> // a "point" is an array of floating points with arithmetic operators
	{
		private:
			using base = std::array<T, Dim>;

		public:

		constexpr point() noexcept = default;

		template <typename ... Ts>
		requires (sizeof...(Ts) <= Dim && sizeof...(Ts) >= 2 && ((std::floating_point<Ts> || std::integral<Ts>) && ...))
		constexpr point(Ts ... ts) noexcept : base{static_cast<T>(ts)...} {}

		constexpr point(T t) noexcept : base(make_point<T, Dim>(t)) {}

		constexpr const point& operator+() const noexcept
		{
			return *this;
		}
		constexpr point operator-() const noexcept
		{
			point res;
			for (std::size_t i = 0; i < Dim; ++i)
				res[i] = -base::operator[](i);
			return res;
		}
		constexpr point& operator+=(const point& other) noexcept
		{
			for (std::size_t i = 0; i < Dim; ++i)
				base::operator[](i) += other[i];
			return *this;
		}
		constexpr point& operator-=(const point& other) noexcept
		{
			for (std::size_t i = 0; i < Dim; ++i)
				base::operator[](i) -= other[i];
			return *this;
		}
		constexpr point& operator*=(const point& other) noexcept
		{
			for (std::size_t i = 0; i < Dim; ++i)
				base::operator[](i) *= other[i];
			return *this;
		}
		constexpr point& operator/=(const point& other) noexcept
		{
			for (std::size_t i = 0; i < Dim; ++i)
				base::operator[](i) /= other[i];
			return *this;
		}
		constexpr point& operator*=(T scal) noexcept
		{
			for (std::size_t i = 0; i < Dim; ++i)
				base::operator[](i) *= scal;
			return *this;
		}
		constexpr point& operator/=(T scal) noexcept
		{
			for (std::size_t i = 0; i < Dim; ++i)
				base::operator[](i) /= scal;
			return *this;
		}
	};
	template <std::floating_point T, std::size_t Dim>
	constexpr point<T, Dim> operator+(point<T, Dim> lhs, const point<T, Dim>& rhs) noexcept
	{
		return lhs += rhs;
	}
	template <std::floating_point T, std::size_t Dim>
	constexpr point<T, Dim> operator-(point<T, Dim> lhs, const point<T, Dim>& rhs) noexcept
	{
		return lhs -= rhs;
	}
	template <std::floating_point T, std::size_t Dim>
	constexpr point<T, Dim> operator*(point<T, Dim> lhs, const point<T, Dim>& rhs) noexcept
	{
		return lhs *= rhs;
	}
	template <std::floating_point T, std::size_t Dim>
	constexpr point<T, Dim> operator/(point<T, Dim> lhs, const point<T, Dim>& rhs) noexcept
	{
		return lhs /= rhs;
	}
	template <std::floating_point T, std::size_t Dim>
	constexpr point<T, Dim> operator*(point<T, Dim> lhs, typename point<T, Dim>::value_type scal) noexcept
	{
		return lhs *= scal;
	}
	template <std::floating_point T, std::size_t Dim>
	constexpr point<T, Dim> operator*(typename point<T, Dim>::value_type scal, point<T, Dim> rhs) noexcept
	{
		return rhs * scal;
	}
	template <std::floating_point T, std::size_t Dim>
	constexpr point<T, Dim> operator/(point<T, Dim> lhs, typename point<T, Dim>::value_type scal) noexcept
	{
		return lhs /= scal;
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

} // namespace physics

#endif // PHYSICS_POINT_H
































