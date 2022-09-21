//  Energy minimizer tests
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

/*

Compilation:
g++ energy_minimizer.cpp -o energy_minimizer -std=c++20 -Wall -Wextra -pedantic -Ofast -fmax-errors=1

*/

#include "periodic_harmonic_oscillator.hpp"
#include "../../physics/energy_minimizer.hpp"

void test_fire()
// test that FIRE algorithm converges for the harmonic oscillator in
// less than 10 iterations with force less than 1e-6 (absolute value)
{
	periodic_harmonic_oscillator sys;
	physics::fire<periodic_harmonic_oscillator> minim(1e-2);

	// initial condition
	sys.x = 1;

	assert(minim.minimize(sys, 1e-6, 10));
}

int main()
{
	test_fire();

	std::cout << "All tests passed successfully!" << std::endl;

	return 0;
}
























