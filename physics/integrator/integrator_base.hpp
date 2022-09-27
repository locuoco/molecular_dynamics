//  Integrator base
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

#ifndef PHYSICS_INTEGRATOR_INTEGRATOR_BASE_H
#define PHYSICS_INTEGRATOR_INTEGRATOR_BASE_H

#include <cmath> // sqrt

#include "../physical_system.hpp" // physical_system, thermodynamical_system

namespace physics
{
	template <physical_system System>
	struct integrator_base
	// `System` is the template argument of the integrator class, corresponding to
	// the class of the system to integrate.
	{
		using system_type = System;
		using state_type = state_type_of<System>;
		using scalar_type = scalar_type_of<System>;

		bool first_step = true;
		// `first_step` is a boolean flag that shall be set to true for the first integration step.

		virtual ~integrator_base() = default;

		virtual void step(System& s, scalar_type dt) = 0;
		// `s` is the system to integrate.
		// `dt` is the integration step.
		// All integrators deriving from `integrator_base` must override this method.

		void simulate(System& s, scalar_type dt, std::size_t n_steps)
		// Integrate for `n_steps` steps. Do not reset at the end.
		// `s` is the system to integrate.
		// `dt` is the integration step.
		// Same as:
		//	for (std::size_t i = 0; i < n_steps; ++i) step(s, dt);
		{
			for (std::size_t i = 0; i < n_steps; ++i)
				step(s, dt);
		}

		void reset()
		// reset integrator (for a new simulation)
		{
			first_step = true;
		}
	};

	template <typename System>
	struct symplectic_integrator_base : virtual integrator_base<System>
	// Inherit from this class if the following assumption about the problem
	// to be solved is thought to be valid:
	// 	H = T(p) + V(x), with f = -grad V(x) and v = grad T(p)
	// Furthermore, H must be explicitly independent on time
	{
		virtual ~symplectic_integrator_base() = default;
	};

	template <typename System>
	struct stochastic_integrator_base : virtual integrator_base<System>
	// Inherit from this class if the integrator is stochastic (which could
	// arise in the resolution of stochastic processes)
	{
		virtual ~stochastic_integrator_base() = default;
	};

	// concepts associated to these template classes (can be used as constraints in template
	// type arguments)
	template <typename Integ>
	concept integrator = std::is_base_of_v<integrator_base<typename Integ::system_type>, Integ>;

	template <typename Integ>
	concept symplectic_integrator = std::is_base_of_v<symplectic_integrator_base<typename Integ::system_type>, Integ>;

	template <typename Integ>
	concept stochastic_integrator = std::is_base_of_v<stochastic_integrator_base<typename Integ::system_type>, Integ>;

	template <integrator Integ>
	requires thermodynamical_system<typename Integ::system_type>
	struct thermodynamical_statistics
	// thermodynamical statistics calculator
	// used to calculate pressure, temperature, volume and density mean values
	// and standard deviations starting from instantaneous values
	{
		using scalar_type = scalar_type_of<typename Integ::system_type>;

		Integ integ;
		scalar_type pressure_sum = 0, pressure_sumsq = 0;
		scalar_type temperature_sum = 0, temperature_sumsq = 0;
		scalar_type volume_sum = 0, volume_sumsq = 0;
		scalar_type density_sum = 0, density_sumsq = 0;
		std::size_t n_steps = 0;

		thermodynamical_statistics(Integ integ = Integ()) : integ(integ) {}
		// constructor:
		//	`integ` is the integrator to use for the simulation

		void step(typename Integ::system_type& s, scalar_type dt)
		// `s` is the system to integrate and from which we want to calculate the statistics.
		// `dt` is the integration step.
		{
			integ.step(s, dt);

			scalar_type quantity;

			quantity = s.pressure();
			pressure_sum += quantity;
			pressure_sumsq += quantity*quantity;

			quantity = s.temperature();
			temperature_sum += quantity;
			temperature_sumsq += quantity*quantity;

			quantity = s.volume();
			volume_sum += quantity;
			volume_sumsq += quantity*quantity;

			quantity = s.density();
			density_sum += quantity;
			density_sumsq += quantity*quantity;

			++n_steps;
		}

		void simulate(typename Integ::system_type& s, scalar_type dt, std::size_t n_steps)
		// Integrate for `n_steps` steps. Do not reset at the end.
		// `s` is the system to integrate and from which we want to calculate the statistics.
		// `dt` is the integration step.
		// Same as:
		//	for (std::size_t i = 0; i < n_steps; ++i) step(s, dt);
		{
			for (std::size_t i = 0; i < n_steps; ++i)
				step(s, dt);
		}

		scalar_type pressure_mean()
		// return pressure mean
		{
			return pressure_sum / n_steps;
		}
		scalar_type pressure_sd()
		// return pressure standard deviation
		{
			using std::sqrt;
			scalar_type mean = pressure_mean();
			return sqrt((pressure_sumsq/n_steps - mean*mean)/(n_steps-1));
		}

		scalar_type temperature_mean()
		// return temperature mean
		{
			return temperature_sum / n_steps;
		}
		scalar_type temperature_sd()
		// return temperature standard deviation
		{
			using std::sqrt;
			scalar_type mean = temperature_mean();
			return sqrt((temperature_sumsq/n_steps - mean*mean)/(n_steps-1));
		}

		scalar_type volume_mean()
		// return volume mean
		{
			return volume_sum / n_steps;
		}
		scalar_type volume_sd()
		// return volume standard deviation
		{
			using std::sqrt;
			scalar_type mean = volume_mean();
			return sqrt((volume_sumsq/n_steps - mean*mean)/(n_steps-1));
		}

		scalar_type density_mean()
		// return density mean
		{
			return density_sum / n_steps;
		}
		scalar_type density_sd()
		// return density standard deviation
		{
			using std::sqrt;
			scalar_type mean = density_mean();
			return sqrt((density_sumsq/n_steps - mean*mean)/(n_steps-1));
		}

		void reset()
		// reset thermodynamical statistics and integrator (for a new simulation)
		{
			pressure_sum = pressure_sumsq = 0;
			temperature_sum = temperature_sumsq = 0;
			volume_sum = volume_sumsq = 0;
			density_sum = density_sumsq = 0;
			n_steps = 0;

			integ.reset();
		}
	};

} // namespace physics

#endif // PHYSICS_INTEGRATOR_INTEGRATOR_BASE_H






























