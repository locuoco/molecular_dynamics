//  Runge-Kutta methods
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

#ifndef PHYSICS_RUNGE_KUTTA_H
#define PHYSICS_RUNGE_KUTTA_H

#include <array>

#include "integrator_base.hpp"

namespace physics
{
	template <having_coordinates System>
	struct runge_kutta_base : virtual integrator_base<System>
	// all Runge-Kutta methods inherit from this class
	{
		virtual ~runge_kutta_base() = default;
	};

	template <typename Integ>
	concept runge_kutta_method = std::is_base_of_v<runge_kutta_base<typename Integ::system_type>, Integ>;

	template <having_coordinates System, std::size_t Stages>
	struct explicit_runge_kutta_base : virtual runge_kutta_base<System>
	// all explicit Runge-Kutta methods inherit from this class
	// `Stages` is the number of stages of the method (number of force evaluations).
	{
		using runge_kutta_base<System>::first_step;

		static constexpr std::size_t num_stages = Stages;

		System& s;

		template <typename ... Ts>
		requires (sizeof...(Ts) == Stages*Stages)
		explicit_runge_kutta_base(System& s, Ts ... pars) : s(s), pars{scalar_type_of<System>(pars)...} {}
		// constructor:
		// 	`s` is the system to integrate.
		//	`pars...` is a variadic argument which contains the parameters for the explicit Runge-Kutta method.

		virtual ~explicit_runge_kutta_base() = default;

		void step(scalar_type_of<System> dt) override
		{
			using std::size_t;

			prev_x = s.x;
			prev_p = s.p;
			prev_t = s.t;

			kx[0] = s.vel();
			kp[0] = s.force(first_step);
			for (size_t i = 1; i < Stages; ++i)
			{
				s.x = 0;
				s.p = 0;
				// sum smaller contributes (in dt) first to minimize rounding errors
				for (size_t j = 0; j < i; ++j)
				{
					scalar_type_of<System> mult = pars[(i-1)*Stages + j+1] * dt;
					s.x += kx[j] * mult;
					s.p += kp[j] * mult;
				}
				s.x += prev_x;
				s.p += prev_p;
				s.t = prev_t + pars[(i-1)*Stages] * dt;

				kx[i] = s.vel();
				kp[i] = s.force();
			}

			s.x = 0;
			s.p = 0;
			for (size_t i = 0; i < Stages; ++i)
			{
				scalar_type_of<System> mult = pars[(Stages-1)*Stages + i] * dt;
				s.x += kx[i] * mult;
				s.p += kp[i] * mult;
			}
			s.x += prev_x;
			s.p += prev_p;
			s.t = prev_t + dt;
			s.force();
			first_step = false;
		}

		// parameters of the Runge-Kutta method
		const std::array<scalar_type_of<System>, Stages*Stages> pars;

		private:

			state_type_of<System> kx[Stages], kp[Stages], prev_x, prev_p;
			scalar_type_of<System> prev_t;
	};

	template <typename Integ>
	concept explicit_runge_kutta_method = std::is_base_of_v<explicit_runge_kutta_base<typename Integ::system_type, Integ::num_stages>, Integ>;

	template <having_coordinates System>
	struct euler : explicit_runge_kutta_base<System, 1>
	// Euler method (1st order, 1 stage)
	{
		euler(System& s)
			: explicit_runge_kutta_base<System, 1>{s, 1}
		{}
	};

	template <having_coordinates System>
	struct rk2 : explicit_runge_kutta_base<System, 2>
	// Parametrized Runge-Kutta 2 method (2nd order, 2 stages)
	// All second-order explicit Runge-Kutta methods can be written
	// as setting the parameter of this method properly
	{
		rk2(System& s, long double a)
		// constructor: `a` is the parameter of the parametrized Runge-Kutta 2 method
			: explicit_runge_kutta_base<System, 2>
			{
				s,
				a, a,
				1 - 1/(2*a), 1/(2*a)
			}
		{}
	};

	template <having_coordinates System>
	struct midpoint : rk2<System>
	// Midpoint method (2nd order, 2 stages)
	{
		midpoint(System& s) : rk2<System>(s, .5L)
		{}
	};

	template <having_coordinates System>
	struct heun2 : rk2<System>
	// Heun method (2nd order, 2 stages)
	{
		heun2(System& s) : rk2<System>(s, 1)
		{}
	};

	template <having_coordinates System>
	struct ralston2 : rk2<System>
	// Ralston method (2nd order, 2 stages)
	{
		ralston2(System& s) : rk2<System>(s, 2/3.L)
		{}
	};

	template <having_coordinates System>
	struct rk4 : explicit_runge_kutta_base<System, 4>
	// Classical Runge-Kutta 4 method (4th order, 4 stages)
	{
		rk4(System& s)
			: explicit_runge_kutta_base<System, 4>
			{
				s,
				.5L, .5L, 0, 0,
				.5L, 0, .5L, 0,
				1, 0, 0, 1,
				1/6.L, 1/3.L, 1/3.L, 1/6.L
			}
		{}
	};

	template <having_coordinates System>
	struct rk4_3_8 : explicit_runge_kutta_base<System, 4>
	// 3/8-rule Runge-Kutta 4 method (4th order, 4 stages)
	{
		rk4_3_8(System& s)
			: explicit_runge_kutta_base<System, 4>
			{
				s,
				1/3.L, 1/3.L, 0, 0,
				2/3.L, -1/3.L, 1, 0,
				1, 1, -1, 1,
				1/8.L, 3/8.L, 3/8.L, 1/8.L
			}
		{}
	};

	template <having_coordinates System>
	struct ralston4 : explicit_runge_kutta_base<System, 4>
	// Ralston method (4th order, 4 stages)
	{
		ralston4(System& s)
			: explicit_runge_kutta_base<System, 4>
			{
				s,
				.4L, .4L, 0, 0,
				.45573725L, .29697761L, .15875964L, 0,
				1, .2181004L, -3.05096516L, 3.83286476L,
				.17476028L, -.55148066L, 1.2055356L, .17118478L
			}
		{}
	};

	template <having_coordinates System>
	struct butcher6 : explicit_runge_kutta_base<System, 7>
	// Butcher method (6th order, 7 stages)
	{
		butcher6(System& s)
			: explicit_runge_kutta_base<System, 7>
			{
				s,
				1/3.L, 1/3.L, 0, 0, 0, 0, 0,
				2/3.L, 0, 2/3.L, 0, 0, 0, 0,
				1/3.L, 1/12.L, 1/3.L, -1/12.L, 0, 0, 0,
				.5L, -1/16.L, 9/8.L, -3/16.L, -3/8.L, 0, 0,
				.5L, 0, 9/8.L, -3/8.L, -3/4.L, .5L, 0,
				1, 9/44.L, -9/11.L, 63/44.L, 18/11.L, -16/11.L, 0,
				11/120.L, 0, 27/40.L, 27/40.L, -4/15.L, -4/15.L, 11/120.L
			}
		{}
	};

	template <having_coordinates System>
	struct verner8 : explicit_runge_kutta_base<System, 11>
	// Verner method (8th order, 11 stages)
	// S. Hippolyte, A. K. Richard, "A New Eighth Order Runge-Kutta Family Method", Journal of Mathematics Research, 2019
	{
		verner8(System& s)
			: explicit_runge_kutta_base<System, 11>
			{
				s,
				.5L, .5L, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				.5L, .25L, .25L, 0, 0, 0, 0, 0, 0, 0, 0,
				(7+s21)/14, 1/7.L, (-7-3*s21)/98, (21+5*s21)/49, 0, 0, 0, 0, 0, 0, 0,
				(7+s21)/14, (11+s21)/84, 0, (18+4*s21)/63, (21-s21)/252, 0, 0, 0, 0, 0, 0,
				.5L, (5+s21)/48, 0, (9+s21)/36, (-231+14*s21)/360, (63-7*s21)/80, 0, 0, 0, 0, 0,
				(7-s21)/14, (10-s21)/42, 0, (-432+92*s21)/315, (633-145*s21)/90, (-504+115*s21)/70, (63-13*s21)/35, 0, 0, 0, 0,
				(7-s21)/14, 1/14.L, 0, 0, 0, (14-3*s21)/126, (13-3*s21)/63, 1/9.L, 0, 0, 0,
				.5L, 1/32.L, 0, 0, 0, (91-21*s21)/576, 11/72.L, (-385-75*s21)/1152, (63+13*s21)/128, 0, 0,
				(7+s21)/14, 1/14.L, 0, 0, 0, 1/9.L, (-733-147*s21)/2205, (515+111*s21)/504, (-51-11*s21)/56, (132+28*s21)/245, 0,
				1, 0, 0, 0, 0, (-42+7*s21)/18, (-18+28*s21)/45, (-273-53*s21)/72, (301+53*s21)/72, (28-28*s21)/45, (49-7*s21)/18,
				1/20.L, 0, 0, 0, 0, 0, 0, 49/180.L, 16/45.L, 49/180.L, 1/20.L
			}
		{}

		private:

			static constexpr long double s21 = 4.582575694955840006588047193728L; // sqrt(21)
	};

	// RUNGE-KUTTA-NYSTROM SCHEMES

	// H = p^2/2m + V(q)
	// only quadratic kinetic energy!

	template <having_coordinates System, std::size_t Stages>
	struct runge_kutta_nystrom_base : virtual runge_kutta_base<System>, virtual symplectic_integrator_base<System>
	// all Runge-Kutta-Nystrom methods inherit from this class
	// `Stages` is the number of stages of the method (number of force evaluations).
	{
		using runge_kutta_base<System>::first_step;

		static constexpr std::size_t num_stages = Stages;

		System& s;

		template <typename ... Ts>
		requires (sizeof...(Ts) == (Stages+1)*Stages)
		runge_kutta_nystrom_base(System& s, Ts ... pars) : s(s), pars{scalar_type_of<System>(pars)...} {}
		// constructor:
		// 	`s` is the system to integrate.
		//	`pars...` is a variadic argument which contains the parameters for the Runge-Kutta-Nystrom method.

		virtual ~runge_kutta_nystrom_base() = default;

		void step(scalar_type_of<System> dt) override
		{
			using std::size_t;

			prev_x = s.x;
			prev_v = s.vel();
			prev_p = s.p;

			k[0] = s.force(first_step);
			for (size_t i = 1; i < Stages; ++i)
			{
				s.x = 0;
				// sum smaller contributes (in dt^2 and in dt) first to minimize rounding errors
				for (size_t j = 0; j < i; ++j)
					s.x += k[j] * (pars[(i-1)*(Stages+1) + j+1] * dt * dt);
				s.x += prev_v * (pars[(i-1)*(Stages+1)] * dt);
				s.x += prev_x;

				k[i] = s.force();
			}

			s.x = 0;
			s.p = 0;
			for (size_t i = 0; i < Stages; ++i)
			{
				s.x += k[i] * (pars[(Stages-1)*Stages + i] * dt * dt);
				s.p += k[i] * (pars[Stages*Stages + i] * dt);
			}
			s.x += prev_v * dt;
			s.x += prev_x;
			s.p += prev_p;
			s.t += dt;
			s.force();
			first_step = false;
		}

		// parameters of the Runge-Kutta-Nystrom method
		const std::array<scalar_type_of<System>, (Stages+1)*Stages> pars;

		private:

			state_type_of<System> k[Stages], prev_x, prev_v, prev_p;
	};

	// TODO (not needed now)

} // namespace physics

#endif // PHYSICS_RUNGE_KUTTA_H























