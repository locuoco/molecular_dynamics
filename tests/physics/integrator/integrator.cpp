//  Symplectic and stochastic integrators tests
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

#include <iostream> // cout, endl
#include <cassert>
#include <cmath> // abs, sqrt, exp, erf, cos
#include <numbers> // numbers::pi

/*

Compilation:
g++ integrator.cpp -o integrator -std=c++20 -Wall -Wextra -pedantic -Ofast -fmax-errors=1

*/

#include "periodic_harmonic_oscillator.hpp"
#include "../../../physics/integrator/integrator.hpp"
#include "../../../math/helper.hpp" // twopi

void test_leapfrog_energy()
// Test conservation of energy of leapfrog method on a harmonic oscillator.
{
	periodic_harmonic_oscillator sys;
	physics::leapfrog<periodic_harmonic_oscillator> integ;

	// initial condition
	sys.x = 1;

	sys.force(); // needed to update potential energy
	double energy0 = sys.total_energy();

	integ.simulate(sys, 1e-3, 10'000'000);

	double energy1 = sys.total_energy();

	assert(std::abs(energy1 - energy0) < 1e-7);
}

void test_symplectic_euler_same()
// Test that symplectic Euler scheme is the same as leapfrog after correcting
// the first and the last step.
{
	periodic_harmonic_oscillator sys1, sys2;
	physics::symplectic_euler<periodic_harmonic_oscillator> integ1;
	physics::leapfrog<periodic_harmonic_oscillator> integ2;

	// initial conditions
	sys1.x = sys2.x = 1;

	unsigned n_steps = 10'000'000;
	double dt = 1e-3;

	// symplectic Euler + corrections
	sys1.p -= sys1.force() * dt/2;
	integ1.simulate(sys1, dt, n_steps);
	sys1.p += sys1.force() * dt/2;

	// leapfrog
	integ2.simulate(sys2, dt, n_steps);

	assert(std::abs(sys1.x - sys2.x) < 1e-12);
	assert(std::abs(sys1.p - sys2.p) < 1e-13);
}

void test_pefrl_order()
// Test that PEFRL method is fourth order.
// The test passes if, by halving the timestep, the accuracy of the final position
// of the numerically integrated harmonic oscillator is 16 times better.
{
	periodic_harmonic_oscillator sys;
	physics::pefrl<periodic_harmonic_oscillator> integ;

	double dt = 1e-2;
	unsigned n_steps = 100'000;
	double x0 = 1, p0 = 0;
	double x_ref = std::cos(std::sqrt(sys.elastic_k/sys.m)*dt*n_steps);
	// initial condition
	sys.x = x0;
	sys.p = p0;

	integ.simulate(sys, dt, n_steps);

	double error1 = std::abs(sys.x - x_ref);

	// halving timestep, doubling number of steps
	dt /= 2;
	n_steps *= 2;
	// initial condition
	sys.x = x0;
	sys.p = p0;

	integ.reset();
	integ.simulate(sys, dt, n_steps);

	double error2 = std::abs(sys.x - x_ref);

	assert(std::abs((error1 - 16*error2) / error1) < 1e-2);
}

void test_vefrl_order()
// Test that VEFRL method is fourth order.
// The test passes if, by halving the timestep, the accuracy of the final position
// of the numerically integrated harmonic oscillator is 16 times better.
{
	periodic_harmonic_oscillator sys;
	physics::vefrl<periodic_harmonic_oscillator> integ;

	double dt = 1e-2;
	unsigned n_steps = 100'000;
	double x0 = 1, p0 = 0;
	double x_ref = std::cos(std::sqrt(sys.elastic_k/sys.m)*dt*n_steps);
	// initial condition
	sys.x = x0;
	sys.p = p0;

	integ.simulate(sys, dt, n_steps);

	double error1 = std::abs(sys.x - x_ref);

	// halving timestep, doubling number of steps
	dt /= 2;
	n_steps *= 2;
	// initial condition
	sys.x = x0;
	sys.p = p0;

	integ.reset();
	integ.simulate(sys, dt, n_steps);

	double error2 = std::abs(sys.x - x_ref);

	assert(std::abs((error1 - 16*error2) / error1) < 1e-2);
}

void test_isokinetic_temperature()
// Test conservation of kinetic energy of the isokinetic integrator on a harmonic oscillator.
// In particular, since we have only 1 degree of freedom, the body should move at constant
// speed (effectively ignoring the potential). That the momentum is conserved is also checked.
{
	periodic_harmonic_oscillator sys;
	physics::isokinetic_leapfrog<periodic_harmonic_oscillator> integ;

	double p0 = 1;
	// initial condition
	sys.x = 1;
	sys.p = p0;

	double temp0 = sys.kinetic_energy();

	integ.simulate(sys, 1e-3, 10'000'000);

	double temp1 = sys.kinetic_energy();
	double p1 = sys.p;

	assert(std::abs(temp1 - temp0) < 1e-6);
	assert(std::abs(p1 - p0) < 1e-6);
}

void test_nose_hoover()
// simulate a canonical (VT) ensemble for the periodic harmonic oscillator using an integrator
// for Nosé-Hoover equations.
// The test passes if the average pressure of the harmonic oscillator corresponds to the
// one calculated analytically for a finite-volume harmonic oscillator (using the
// canonical ensemble), and if the average temperature corresponds to the reference
// temperature.
{
	periodic_harmonic_oscillator sys;
	physics::nose_hoover<periodic_harmonic_oscillator> integ(10, 7e-2);

	double a = sys.elastic_k * sys.side*sys.side / (8*sys.kT_ref());
	double pressure_analytical = std::sqrt(sys.elastic_k * sys.kT_ref() / math::two_pi<>) * std::exp(-a) / std::erf(std::sqrt(a));

	// initial condition
	sys.x = 1;

	unsigned n_steps = 10'000'000;
	double sump = 0, sumtemp = 0;

	for (unsigned i = 0; i < n_steps; ++i)
	{
		integ.step(sys, 1e-3);

		sump += sys.pressure();
		sumtemp += sys.temperature();
	}

	assert(std::abs((sump/n_steps - pressure_analytical) / pressure_analytical) < 1e-2);
	assert(std::abs((sumtemp/n_steps - sys.temperature_ref) / sys.temperature_ref) < 1e-2);
}

void test_mtk()
// simulate an isothermal-isobaric (PT) ensemble for the periodic harmonic oscillator using an
// integrator for Martyna-Tobias-Klein (MTK) equations.
// The test passes if:
// * the average volume of the harmonic oscillator corresponds to the one calculated analytically
//   (using the isothermal-isobaric ensemble);
// * the average temperature corresponds to the reference temperature;
// * the average pressure corresponds to the reference pressure.
{
	periodic_harmonic_oscillator sys;
	physics::mtk<periodic_harmonic_oscillator> integ(5, 5e-2, 2e-2);

	double a = 2 * sys.pressure_ref*sys.pressure_ref / (sys.elastic_k * sys.kT_ref());
	double volume_analytical
		= sys.kT_ref()/sys.pressure_ref
		- 4*sys.pressure_ref/sys.elastic_k
		+ 2*std::sqrt(2/(std::numbers::pi*sys.elastic_k*sys.kT_ref()))*std::exp(-a)/std::erfc(std::sqrt(a));

	// initial condition
	sys.x = 1;

	unsigned n_steps = 10'000'000;
	double sump = 0, sumtemp = 0, sumvol = 0;

	for (unsigned i = 0; i < n_steps; ++i)
	{
		integ.step(sys, 1e-3);

		sump += sys.pressure();
		sumtemp += sys.temperature();
		sumvol += sys.volume();
	}

	assert(std::abs((sumvol/n_steps - volume_analytical) / volume_analytical) < 1e-2);
	assert(std::abs((sumtemp/n_steps - sys.temperature_ref) / sys.temperature_ref) < 1e-2);
	assert(std::abs((sump/n_steps - sys.pressure_ref) / sys.pressure_ref) < 1e-2);
}

int main()
{
	test_leapfrog_energy();
	test_symplectic_euler_same();
	test_pefrl_order();
	test_vefrl_order();
	test_isokinetic_temperature();
	test_nose_hoover();
	test_mtk();

	std::cout << "All tests passed successfully!" << std::endl;

	return 0;
}
























