//  Physical system
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

#ifndef PHYSICS_PHYSICAL_SYSTEM_H
#define PHYSICS_PHYSICAL_SYSTEM_H

#include <concepts> // convertible_to
#include <type_traits> // is_base_of

namespace physics
{
	template <typename T, typename State>
	struct physical_system_base
	// any class which wants to use any of the integrators implemented by this library
	// must derive from `physical_system_base` (an abstract class).
	// `force` and `vel` methods must be available with the required signature.
	// The `eval` argument will be true if a new evaluation of the forces or the velocities is probably required,
	// and false if it is not (because it would give the same result as the last evaluation). If the actual
	// values are stored internally, the implementation may take advantage of this.
	{
		using scalar_type = T;
		using state_type = State;

		virtual const State& force(bool eval = true) = 0;
		virtual const State& vel(bool eval = true) = 0;
		virtual ~physical_system_base() = default;
	};

	template <typename System>
	using scalar_type_of = typename System::scalar_type;

	template <typename System>
	using state_type_of = typename System::state_type;

	template <typename System>
	concept physical_system = std::is_base_of_v<physical_system_base<scalar_type_of<System>, state_type_of<System>>, System>;
	// A `physical_system` is any class that is derived from a `physical_system_base` template class.

	template <typename System>
	concept having_coordinates = physical_system<System>
		&& requires(System& s, scalar_type_of<System> dt)
	// A system `having_coordinates` requires to be a physical system
	// and has `x` (positions), `p` (momenta) and `t` (time coordinate) defined.
	// The following operations are required to be defined.
		{
			s.x += s.vel() * dt;
			s.p += s.force() * dt;
			s.t += dt;
		};

	template <typename System>
	concept thermodynamical_system = physical_system<System>
		&& requires(System& s, scalar_type_of<System> t)
	// A `thermodynamical_system` requires to be a physical system and has the following
	// methods defined, whose return values have to be convertible to the scalar type.
		{
			t += s.pressure();
			t += s.temperature();
			t += s.volume();
			t += s.density();
		};

	// abstract template class that can be inherited and used to satisfy the previous constraints.
	template <typename T, typename State>
	struct coordinates_system_base : physical_system_base<T, State>
	{
		State x, p;
		T t;
	};

} // namespace physics

#endif // PHYSICS_PHYSICAL_SYSTEM_H






























