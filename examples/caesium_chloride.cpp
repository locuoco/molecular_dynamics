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

#include <sstream>
#include <limits> // infinity

/*

Compilation (MinGW):
g++ caesium_chloride.cpp -o caesium_chloride -std=c++20 -I C:\Users\aless\Desktop\myLib\include -L C:\Users\aless\Desktop\myLib\lib -lopengl32 -lglu32 -lglew32.dll -lglfw3dll -lfreetype -Wall -Wextra -pedantic -Ofast -pthread -fmax-errors=1

*/

#include "../physics/physics.hpp"
#include "../gui/graphics.hpp"

using namespace physics::literals;

int main()
{
	graphics window;

	unsigned n_side = 12; // 12x12x12 = 1728 unit cubic cells = 3456 ions

	physics::molecular_system molsys;
	physics::nose_hoover integ(molsys);

	// adding more cells to improve PPPM accuracy
	molsys.lrsum.cell_multiplier(1);

	// CsCl is a primitive cubic lattice with two elements
	molsys.primitive_cubic_lattice(n_side, physics::cscl_lattice<>, physics::caesium_ion<>, physics::chloride_ion<>);

	while (!window.should_close())
	{
		integ.step(4_fs);
		std::stringstream custom_text;
		custom_text << "Ewald parameter = " << molsys.lrsum.ewald_par();
		custom_text << " A^-1\nEstimated electrostatic force RMSE = " << molsys.lrsum.estimated_error_coulomb;
		custom_text << " kcal/(mol A)\nEstimated dispersion force RMSE = " << molsys.lrsum.estimated_error_lj;
		custom_text << " kcal/(mol A)\nEstimated total force RMSE = " << molsys.lrsum.estimated_error;

		// the primitive cubic structure is stable!
		window.draw(molsys, custom_text.str());
	}

	return 0;
}





















