//  Particle-particle, particle-mesh method tests
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
#include <cmath> // abs

/*

Compilation:
g++ pppm.cpp -o pppm -std=c++20 -Wall -Wextra -pedantic -Ofast -pthread -fmax-errors=1

*/

#include "../../physics/physics.hpp"

using namespace physics::literals;

void test_charge_assignment_function(size_t order, double x)
// test property of charge assignment function, i.e.:
// sum_k W_k (x) = 1  for all x (charge must be conserved)
// with k = 0, ..., order-1
// `order` is the order of the charge assignment function to test.
// `x` is the argument to feed to the function to test.
{
	double sum = 0;
	for (std::size_t k = 0; k < order; ++k)
		sum += physics::charge_assignment_function(x, k, order);

	assert(std::abs(sum - 1) < 1e-15);
}

void test_force_accuracy(std::string scheme)
// Test PPPM algorithm against Ewald summation for a system of TIP3P water molecules.
// The test passes if the RMS error is not greater than 3 times the estimated RMS error.
// `scheme` is the scheme to use for PPPM (either "ik" or "ad").
{
	physics::molecular_system<double, physics::pppm> sys;
	physics::molecular_system<double, physics::ewald> sys_ref;
	physics::leapfrog<physics::molecular_system<double, physics::pppm>> integ;

	// Set parameters in order to increase accuracy
	// Also, cutoff radius must be large enough so that LJ truncation error is negligible
	sys.lrsum.cutoff_radius(20);
	sys.lrsum.charge_assignment_order(7);
	sys.lrsum.cell_multiplier(2);
	sys.lrsum.precise(true);
	// Set differentiation scheme
	sys.lrsum.set_diff_scheme(scheme);

	sys_ref.lrsum.max_n(12);
	sys_ref.lrsum.precise(true);

	double dist = std::cbrt(physics::water_mass<> / physics::water_density25<>);

	sys.primitive_cubic_lattice(8, dist, physics::water_tip3p<>);
	// simulate for 20 steps, to break the perfect-lattice structure
	integ.simulate(sys, 1_fs, 20);

	sys_ref = sys; // copy content

	sys.force();
	sys_ref.force();

	double error = std::sqrt(sys.lrsum.estimated_error * sys.lrsum.estimated_error
	                       + sys_ref.lrsum.estimated_error * sys_ref.lrsum.estimated_error);
	assert(rms(sys.f - sys_ref.f) < 3*error); // should be around 1e-7
}

void test_energy_same()
// Test that the energy calculated from ik- and ad-differentiation with the same
// parameters are the same up to rounding error.
{
	double dist = std::cbrt(physics::water_mass<> / physics::water_density25<>);

	physics::molecular_system<double, physics::pppm> sys;
	sys.primitive_cubic_lattice(3, dist, physics::water_tip3p<>);

	sys.lrsum.set_diff_scheme("ik");
	sys.force();
	double ik_energy = sys.potential;
	double param = sys.lrsum.ewald_par(); // get automatically optimized Ewald parameter

	sys.lrsum.set_diff_scheme("ad");
	sys.lrsum.ewald_par(param); // set same Ewald parameter manually

	sys.force();
	double ad_energy = sys.potential;

	assert(std::abs(ik_energy - ad_energy) < 1e-15);
}

void test_madelung_nacl()
// Test that the calculated electrostatic energy for a sodium chloride lattice is related to the
// Madelung constant M = -1.747565 with this formula:
// E = k_C e^2 N M / a
// where `k_C` is the Coulomb constant, `e` is the elementary charge, `a` is the NaCl
// lattice constant and `N` is the number of atoms. Note that in AKMA units e = 1.
{
	physics::molecular_system<double, physics::pppm> sys;

	// set dielectric to infinity to ignore dipole correction
	sys.lrsum.dielectric(std::numeric_limits<double>::infinity());
	// Set parameters in order to increase accuracy
	sys.lrsum.charge_assignment_order(7);
	sys.lrsum.cutoff_radius(20);
	sys.lrsum.cell_multiplier(2);

	sys.face_centered_cubic_lattice(2, physics::nacl_lattice<>, physics::sodium_ion<>, physics::chloride_ion<>);

	sys.force();

	double calculated_madelung = sys.lrsum.energy_coulomb*physics::nacl_lattice<>/(physics::kC<>*sys.n);
	assert(std::abs(calculated_madelung - physics::nacl_madelung<>) < 1e-5);
}

void test_madelung_cscl()
// Test that the calculated electrostatic energy for a caesium chloride lattice is related to the
// Madelung constant M = -1.762675 with this formula:
// E = k_C e^2 N M / (4 sqrt(3) a)
// where `k_C` is the Coulomb constant, `e` is the elementary charge, `a` is the CsCl
// lattice constant and `N` is the number of atoms. Note that in AKMA units e = 1.
{
	physics::molecular_system<double, physics::pppm> sys;

	// set dielectric to infinity to ignore dipole correction
	sys.lrsum.dielectric(std::numeric_limits<double>::infinity());
	// Set parameters in order to increase accuracy
	sys.lrsum.charge_assignment_order(7);
	sys.lrsum.cutoff_radius(20);
	sys.lrsum.cell_multiplier(1);

	sys.primitive_cubic_lattice(8, physics::cscl_lattice<>, physics::caesium_ion<>, physics::chloride_ion<>);

	sys.force();

	double calculated_madelung = sys.lrsum.energy_coulomb*physics::cscl_lattice<>*std::numbers::sqrt3/(physics::kC<>*sys.n);
	assert(std::abs(calculated_madelung - physics::cscl_madelung<>) < 1e-4);
}

int main()
{
	for (auto i = physics::pppm_min_order; i <= physics::pppm_max_order; ++i)
		test_charge_assignment_function(i, 0.3141592); // a random number
	test_force_accuracy("ik");
	test_force_accuracy("ad");
	test_energy_same();
	test_madelung_nacl();
	test_madelung_cscl();

	std::cout << "All tests passed successfully!" << std::endl;

	return 0;
}
























