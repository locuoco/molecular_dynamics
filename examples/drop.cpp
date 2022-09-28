//  Molecular dynamics example: water drop
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
#include <cmath> // sqrt

/*

Compilation (MinGW):
g++ drop.cpp -o drop -std=c++20 -I C:\Users\aless\Desktop\myLib\include -L C:\Users\aless\Desktop\myLib\lib -lopengl32 -lglu32 -lglew32.dll -lglfw3dll -lfreetype -Wall -Wextra -pedantic -Ofast -pthread -fmax-errors=1

*/

#include "../physics/physics.hpp"
#include "../gui/graphics.hpp"

using namespace physics::literals;

int main()
{
	graphics window;

	double dist = std::cbrt(physics::water_mass<> / physics::water_density25<>);

	unsigned n_side = 4; // 4x4x4 = 64 water molecules

	physics::molecular_system molsys;
	physics::isokinetic_leapfrog integ(molsys);

	molsys.temperature_ref = 273.15+50;
	molsys.lrsum.cell_multiplier(1);

	// set a simple cubic lattice as initial condition
	molsys.primitive_cubic_lattice(n_side, dist, physics::water_fba_eps<>);

	molsys.side *= 3;

	while (!window.should_close())
	{
		std::stringstream custom_text;

		integ.simulate(1_fs, 10);

		custom_text << "Ewald parameter = " << molsys.lrsum.ewald_par();
		custom_text << " A^-1\nEstimated electrostatic force RMSE = " << molsys.lrsum.estimated_error_coulomb;
		custom_text << " kcal/(mol A)\nEstimated dispersion force RMSE = " << molsys.lrsum.estimated_error_lj;
		custom_text << " kcal/(mol A)\nEstimated total force RMSE = " << molsys.lrsum.estimated_error;

		// draw all steps in real time
		window.draw(molsys, custom_text.str());
	}

	return 0;
}





















