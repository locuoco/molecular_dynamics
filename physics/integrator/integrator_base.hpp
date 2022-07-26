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
	template <typename T, typename State>
	struct integrator_base {};

	// All integrators derive from `integrator_base` and will have a template method
	// called `step` with the following signature:

	//		template <constraint<T, State> System>
	//		void step(System& s, T dt, bool first_step = false);

	// where `T` and `State` are template arguments of the integrator class,
	// while `constraint` is a template concept that puts constraints to the
	// generic `System` class, depending on the integrator (see also physical_system.hpp).
	// `dt` is the integration step and `first_step` is a boolean flag that shall
	// be set to true for the first integration step.

	template <typename T, typename State>
	struct symplectic_integrator_base : integrator_base<T, State> {};
	// Inherit from this class if the following assumption about the problem
	// to be solved is thought to be valid:
	// 	H = T(p) + V(x), with f = -grad V(x) and v = grad T(p)
	// Furthermore, H must be explicitly independent on time

	template <typename T, typename State>
	struct stochastic_integrator_base : integrator_base<T, State> {};
	// Inherit from this class if the integrator is stochastic (which could
	// arise in the resolution of stochastic processes)

	// concepts associated to these template classes (can be used as constraints in template
	// type arguments)
	template <typename Integ, typename T, typename State>
	concept integrator = std::is_base_of_v<integrator_base<T, State>, Integ>;

	template <typename Integ, typename T, typename State>
	concept symplectic_integrator = std::is_base_of_v<symplectic_integrator_base<T, State>, Integ>;

	template <typename Integ, typename T, typename State>
	concept stochastic_integrator = std::is_base_of_v<stochastic_integrator_base<T, State>, Integ>;

} // namespace physics

#endif // PHYSICS_INTEGRATOR_INTEGRATOR_BASE_H






























