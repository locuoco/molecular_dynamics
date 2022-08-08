//  Molecular dynamics example: caesium chloride lattice
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

#include <iostream>
#include <vector>
#include <limits> // infinity

/*

Compilation:
g++ caesium_chloride.cpp -o caesium_chloride -std=c++20 -I C:\Users\aless\Desktop\myLib\include -L C:\Users\aless\Desktop\myLib\lib -lopengl32 -lglu32 -lglew32.dll -lglfw3dll -lfreetype -Wall -Wextra -pedantic -Ofast -pthread -fmax-errors=1

*/

#include "../physics/physics.hpp"
#include "../gui/graphics.hpp"

int main()
{
	graphics window;

	unsigned n_side = 12;

	physics::molecular_system<double, physics::pppm, physics::nose_hoover> molsys;

	// set dielectric to infinity to ignore dipole correction
	molsys.lrsum.dielectric(std::numeric_limits<double>::infinity());

	// ik-differentiation scheme conserves momentum better
	molsys.lrsum.set_diff_scheme("ik");
	molsys.lrsum.cutoff_radius(7);
	molsys.lrsum.cell_multiplier(1);

	molsys.primitive_cubic_lattice(n_side, physics::cscl_lattice<>, physics::caesium_ion<>, physics::chloride_ion<>);

	while (!window.should_close())
	{
		molsys.step(4e-3);
		window.draw(molsys);
	}

	return 0;
}





















