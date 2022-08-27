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

#include "../physical_system.hpp"

namespace physics
{
	template <typename System>
	struct integrator_base
	{
		bool first_step = true;

		virtual ~integrator_base() = default;

		virtual void step(System& s, scalar_type_of<System> dt) = 0;

		void simulate(System& s, scalar_type_of<System> dt, std::size_t n_steps)
		// Integrate for `n_steps` steps.
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

	// All integrators derive from `integrator_base` and will have a method
	// called `step` with the following signature:

	//		void step(System& s, scalar_type_of<System> dt);

	// where `System` is the template argument of the integrator class, corresponding to
	// the class of the system to integrate.
	// `dt` is the integration step and `first_step` is a boolean flag that shall
	// be set to true for the first integration step.

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
	template <typename Integ, typename System>
	concept integrator = std::is_base_of_v<integrator_base<System>, Integ>;

	template <typename Integ, typename System>
	concept symplectic_integrator = std::is_base_of_v<symplectic_integrator_base<System>, Integ>;

	template <typename Integ, typename System>
	concept stochastic_integrator = std::is_base_of_v<stochastic_integrator_base<System>, Integ>;

} // namespace physics

#endif // PHYSICS_INTEGRATOR_INTEGRATOR_BASE_H






























