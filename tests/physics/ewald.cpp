//  Ewald summation tests
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
#include <iomanip> // setprecision
#include <cassert>
#include <cmath> // abs

/*

Compilation:
g++ ewald.cpp -o ewald -std=c++20 -Wall -Wextra -pedantic -Ofast -pthread -fmax-errors=1

*/

#include "../../physics/physics.hpp"
#include "nist_test.hpp" // load_nist_sample

void test_accuracy_nist()
// test Ewald summation and dispersion calculation against reference data.
// Samples and data taken from:
// https://www.nist.gov/mml/csd/chemical-informatics-group/spce-water-reference-calculations-10a-cutoff
{
	using std::abs;
	physics::molecular_system<double, physics::ewald> sys;
	// set parameters used in reference calculation
	sys.lrsum.ewald_par(5.6/20);
	// The reference calculation used a reciprocal-space spherical cutoff (n^2 < 27), while
	// this library use a cubic cutoff. Using a cubic cutoff of 4 leads to a very similar result.
	sys.lrsum.max_n(4);

	load_nist_sample(sys);

	sys.force();

	assert(sys.side[0] == 20 && sys.side[1] == 20 && sys.side[2] == 20);
	assert(sys.n == 300); // there are 100 water molecules in the sample (i.e. 300 atoms)

	// total energy
	double tot_energy_ref = -4.88604e5;
	assert(abs((sys.potential/physics::kB<> - tot_energy_ref) / tot_energy_ref) < 1e-5);

	// reciprocal-space energy contribution
	// this value is slightly off since the reciprocal-space cutoff is done differently in this library
	// w.r.t. the one used for the reference calculation
	double reciprocal_energy_ref = 6.27009e3;
	assert(abs((sys.lrsum.energy_k/physics::kB<> - reciprocal_energy_ref) / reciprocal_energy_ref) < 1e-4);

	// self-energy correction contribution
	double self_energy_ref = -2.84469e6;
	assert(abs((sys.lrsum.energy_scd/physics::kB<> - self_energy_ref) / self_energy_ref) < 1e-6);

	// dispersion interactions contributions (both attractive and repulsive, excluding intra-molecular intractions)
	double lj_energy_ref = 9.95387e4;
	assert(abs(((sys.lrsum.energy_lj-sys.energy_intra_lj)/physics::kB<> - lj_energy_ref) / lj_energy_ref) < 1e-6);

	// long-range dispersion correction (both attractive and repulsive)
	double long_range_correction_ref = -8.23715e2;
	assert(abs((sys.energy_lrc/physics::kB<> - long_range_correction_ref) / long_range_correction_ref) < 1e-6);

	// real-space energy contribution (excluding intra-molecular interactions)
	// the error for this term is presumably dominated by the limited precision of `math::fasterfc`, but
	// it is acceptable enough. The reference calculation used the "Numerical recipes" implementation for erfcc,
	// which has similar accuracy to the implementation used for the present library.
	double real_energy_ref = -5.58889e5 + 2.80999e6;
	assert(abs(((sys.lrsum.energy_r+sys.energy_intra_coulomb)/physics::kB<> - real_energy_ref) / real_energy_ref) < 1e-5);
}

void test_madelung_nacl()
// Test that the calculated electrostatic energy for a sodium chloride lattice is related to the
// Madelung constant M = 1.747565 with this formula:
// E = -k_C e^2 N M / a
// where `k_C` is the Coulomb constant, `e` is the elementary charge, `a` is the NaCl
// lattice constant and `N` is the number of atoms. Note that in AKMA units e = 1.
{
	physics::molecular_system<double, physics::ewald> sys;

	sys.face_centered_cubic_lattice(2, physics::nacl_lattice<>, physics::sodium_ion<>, physics::chloride_ion<>);

	sys.force();

	double calculated_madelung = -sys.lrsum.energy_coulomb*physics::nacl_lattice<>/(physics::kC<>*sys.n);

	std::cout << std::setprecision(13);
	std::cout << calculated_madelung << '\n';
	std::cout << physics::nacl_madelung<> << '\n';

	assert(std::abs(calculated_madelung - physics::nacl_madelung<>) < 1e-6);
}

void test_madelung_cscl()
// Test that the calculated electrostatic energy for a caesium chloride lattice is related to the
// Madelung constant M = 1.762675 with this formula:
// E = -k_C e^2 N M / (sqrt(3) a)
// where `k_C` is the Coulomb constant, `e` is the elementary charge, `a` is the CsCl
// lattice constant and `N` is the number of atoms. Note that in AKMA units e = 1.
{
	physics::molecular_system<double, physics::ewald> sys;

	sys.primitive_cubic_lattice(2, physics::cscl_lattice<>, physics::caesium_ion<>, physics::chloride_ion<>);

	sys.force();

	double calculated_madelung = -sys.lrsum.energy_coulomb*physics::cscl_lattice<>*std::numbers::sqrt3/(physics::kC<>*sys.n);

	std::cout << std::setprecision(13);
	std::cout << calculated_madelung << '\n';
	std::cout << physics::cscl_madelung<> << '\n';

	assert(std::abs(calculated_madelung - physics::cscl_madelung<>) < 1e-5);
}

int main()
{
	test_accuracy_nist();
	test_madelung_nacl();
	test_madelung_cscl();

	std::cout << "All tests passed successfully!" << std::endl;

	return 0;
}
























