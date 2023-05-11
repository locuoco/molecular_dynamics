//  Molecule tests
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
g++ molecule.cpp -o molecule -std=c++20 -Wall -Wextra -pedantic -Ofast -pthread -fmax-errors=1

*/

#include "../../physics/molecule.hpp"
#include "../../physics/tensor.hpp" // rms

constexpr double error_threshold = 1e-15;

void test_water_mass()
// Test that molecular mass for all water models is equal to `water_mass<>` = 18.0154 u.
{
	using std::abs;
	assert(abs(mass_of(physics::water_tip3p<>) - physics::water_mass<>) < error_threshold);
	assert(abs(mass_of(physics::water_tip3p_lr<>) - physics::water_mass<>) < error_threshold);
	assert(abs(mass_of(physics::water_fba_eps<>) - physics::water_mass<>) < error_threshold);
	assert(abs(mass_of(physics::water_spc_e<>) - physics::water_mass<>) < error_threshold);
}

void test_water_charge()
// Test that the total charge for all water models is equal to 0.
{
	using std::abs;
	assert(abs(charge_of(physics::water_tip3p<>)) < error_threshold);
	assert(abs(charge_of(physics::water_tip3p_lr<>)) < error_threshold);
	assert(abs(charge_of(physics::water_fba_eps<>)) < error_threshold);
	assert(abs(charge_of(physics::water_spc_e<>)) < error_threshold);
}

void test_pc_lattice()
// Test that the simulation box length is correct after calling
// `primitive_cubic_lattice`
{
	int n_side = 5;
	double lattice_constant = 3.333;
	physics::molecular_system sys;

	sys.primitive_cubic_lattice(n_side, lattice_constant, physics::water_tip3p<>);

	assert(std::abs(sys.side[0] - n_side * lattice_constant) < error_threshold);
	assert(std::abs(sys.side[1] - n_side * lattice_constant) < error_threshold);
	assert(std::abs(sys.side[2] - n_side * lattice_constant) < error_threshold);
}

void test_fcc_lattice()
// Test that the simulation box length is correct after calling
// `face_centered_cubic_lattice`
{
	int n_side = 5;
	double lattice_constant = 3.333;
	physics::molecular_system sys;

	sys.face_centered_cubic_lattice(n_side, lattice_constant, physics::water_tip3p<>);

	assert(std::abs(sys.side[0] - n_side * lattice_constant) < error_threshold);
	assert(std::abs(sys.side[1] - n_side * lattice_constant) < error_threshold);
	assert(std::abs(sys.side[2] - n_side * lattice_constant) < error_threshold);
}

void test_cubic_lattice()
// Test that a face-centered cubic lattice with species1 == species2 degenerates
// into a primitive cubic lattice with only one species.
{
	physics::molecular_system sys1, sys2;
	auto my_molecule = physics::sodium_ion<>;

	sys1.face_centered_cubic_lattice(1, 2, my_molecule, my_molecule);
	sys2.primitive_cubic_lattice(2, 1, my_molecule);

	assert(sys1.n == sys2.n);
	assert(rms(sys1.x - sys2.x) < error_threshold);
	// note that test assumes the atom positions are put in the same order for
	// both sys1 and sys2!
}

void test_num_bonds()
// test that the number of elements of the bonds member variable of the molecules
// corresponds to the number of atoms `n` in the molecule.
{
	assert(physics::water_tip3p<>.bonds.size() == physics::water_tip3p<>.n);
	assert(physics::water_tip3p_lr<>.bonds.size() == physics::water_tip3p_lr<>.n);
	assert(physics::water_fba_eps<>.bonds.size() == physics::water_fba_eps<>.n);
	assert(physics::water_spc_e<>.bonds.size() == physics::water_spc_e<>.n);
	assert(physics::sodium_ion<>.bonds.size() == physics::sodium_ion<>.n);
	assert(physics::chloride_ion<>.bonds.size() == physics::chloride_ion<>.n);
	assert(physics::caesium_ion<>.bonds.size() == physics::caesium_ion<>.n);
}

int main()
{
	test_water_mass();
	test_water_charge();
	test_pc_lattice();
	test_fcc_lattice();
	test_cubic_lattice();
	test_num_bonds();

	std::cout << "All tests passed successfully!" << std::endl;

	return 0;
}
























