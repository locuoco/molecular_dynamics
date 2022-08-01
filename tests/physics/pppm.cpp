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
#include "../../physics/tensor.hpp" // rms

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
// Test PPPM algorithm against Ewald summation for a lattice of TIP3P water molecules.
// The test passes if the RMS error is not greater than 3 times the estimated error.
// `scheme` is the scheme to use for PPPM (either "ik" or "ad").
{
	physics::molecular_system<double, physics::pppm> sys;
	physics::molecular_system<double, physics::ewald> sys_ref;

	// Set parameters in order to increase accuracy
	sys.lrsum.cutoff_radius(9);
	sys.lrsum.charge_assignment_order(7);
	// Set differentiation scheme
	sys.lrsum.set_diff_scheme(scheme);

	double dist = std::cbrt(physics::water_mass<> / physics::water_density25<>);

	sys.primitive_cubic_lattice(16, dist, physics::water_tip3p<>);

	sys_ref = sys;

	assert(rms(sys.force() - sys_ref.force()) < 3*sys.lrsum.estimated_error);
}

void test_madelung_nacl()
// Test that the calculated electrostatic energy for NaCl lattice is related to the
// Madelung constant M = -1.747565 with this formula:
// E = k_C e^2 M / a
// where `k_C` is the Coulomb constant, `e` is the elementary charge and `a` is the
// NaCl lattice constant. Note that in AKMA units e = 1.
{
	physics::molecular_system<double, physics::pppm> sys;

	// Set parameters in order to increase accuracy
	sys.lrsum.charge_assignment_order(7);
	sys.lrsum.cutoff_radius(20);
	sys.lrsum.cell_multiplier(2);

	int n_side = 2;
	double lattice_constant = 5.6402;

	sys.face_centered_cubic_lattice(n_side, lattice_constant, physics::sodium_ion<>, physics::chloride_ion<>);

	sys.force();
	double calculated_madelung = sys.lrsum.energy_coulomb*lattice_constant/(physics::kC<>*sys.n);

	assert(std::abs(calculated_madelung - physics::madelung_nacl<>) < 1e-5);
}

int main()
{
	for (auto i = physics::pppm_min_order; i <= physics::pppm_max_order; ++i)
		test_charge_assignment_function(i, 0.3141592); // a random number
	test_force_accuracy("ik");
	test_force_accuracy("ad");
	test_madelung_nacl();

	std::cout << "All tests passed successfully!" << std::endl;

	return 0;
}
























