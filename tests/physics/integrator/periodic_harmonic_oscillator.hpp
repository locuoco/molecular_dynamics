//  Periodic harmonic oscillator class
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

#ifndef TESTS_PHYSICS_INTEGRATOR_PERIODIC_HARMONIC_OSCILLATOR
#define TESTS_PHYSICS_INTEGRATOR_PERIODIC_HARMONIC_OSCILLATOR

#include <cmath> // remainder

#include "../../../physics/physical_system.hpp"
#include "../../../physics/tensor.hpp" // scalard

struct periodic_harmonic_oscillator : physics::physical_system_base<double, physics::scalard>
// 1-d harmonic oscillator with periodic boundaries
// having finite volume, the pressure will be greater than 0
{
	physics::scalard x = 0, p = 0, v = 0, f = 0; // oscillator coordinates and derivatives
	physics::scalard elastic_k = 1, m = 1; // elastic constant and mass
	double potential = 0, virial = 0; // potential and virial
	double t = 0; // time
	double side = 4; // length of the system (centered at the origin)
	double temperature_ref = 1, pressure_ref = 1; // reference temperature and pressure

	static constexpr std::size_t n = 1, dof = 1; // number of particles and degrees of freedom

	const physics::scalard& vel(bool = true) override
	// return velocity (calculated from momentum).
	// Boolean argument is ignored.
	{
		v = p / m;
		return v;
	}

	const physics::scalard& force(bool eval = true) override
	// return the force
	// if `eval` is false, reuse old value
	// otherwise, it also updates potential and virial
	{
		if (eval)
		{
			// remap x in the range [-side/2, side/2]
			x = std::remainder(x, side);
			f = -elastic_k * x;
			potential = elastic_k * x*x / 2;
			virial = f*x;
		}
		return f;
	}

	double kinetic_energy() const noexcept
	{
		return p*p / (2*m);
	}

	double total_energy() const noexcept
	// return the total energy
	// potential energy is the one calculated in the `force` method
	{
		return kinetic_energy() + potential;
	}

	double temperature() const noexcept
	// Instantaneous temperature of the oscillator
	// Boltzmann constant is assumed to be 1.
	{
		return 2*kinetic_energy();
	}

	double volume() const noexcept
	// Return the volume of the oscillator.
	// It is a length since we are in 1-d.
	{
		return side;
	}

	double pressure() const noexcept
	// Instantaneous pressure calculated with instantaneous temperature and virial
	{
		return (2*kinetic_energy() + virial)/volume();
	}

	double kT_ref() const noexcept
	// Return the reference temperature multiplied by the Boltzmann constant.
	// Boltzmann constant is assumed to be 1.
	{
		return temperature_ref;
	}

	double kT() const noexcept
	// Return the instantaneous temperature multiplied by the Boltzmann constant.
	// Boltzmann constant is assumed to be 1.
	{
		return temperature();
	}
};

#endif // TESTS_PHYSICS_INTEGRATOR_PERIODIC_HARMONIC_OSCILLATOR






















