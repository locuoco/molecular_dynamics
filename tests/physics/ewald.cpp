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
#include <cassert>
#include <cmath> // abs

/*

Compilation:
g++ ewald.cpp -o ewald -std=c++20 -Wall -Wextra -pedantic -Ofast -pthread -fmax-errors=1

*/

#include "../../physics/physics.hpp"

void test_madelung_nacl()
// Test that the calculated electrostatic energy for NaCl lattice is related to the
// Madelung constant M = -1.747565 with this formula:
// E = k_C e^2 M / a
// where `k_C` is the Coulomb constant, `e` is the elementary charge and `a` is the
// NaCl lattice constant. Note that in AKMA units e = 1.
{
	physics::molecular_system<double, physics::ewald> sys;

	int n_side = 2;
	double lattice_constant = 5.6402;

	sys.face_centered_cubic_lattice(n_side, lattice_constant, physics::sodium_ion<>, physics::chloride_ion<>);

	sys.force();
	double calculated_madelung = sys.lrsum.energy_coulomb*lattice_constant/(physics::kC<>*sys.n);

	assert(std::abs(calculated_madelung - physics::madelung_nacl<>) < 1e-5);
}

int main()
{
	test_madelung_nacl();

	std::cout << "All tests passed successfully!" << std::endl;

	return 0;
}
























