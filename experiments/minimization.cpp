//  Molecular dynamics example: minimization
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

#include <cmath>
#include <sstream>

/*

Compilation (MinGW):
g++ minimization.cpp -o minimization -std=c++20 -I C:\Users\aless\Desktop\myLib\include -L C:\Users\aless\Desktop\myLib\lib -lopengl32 -lglu32 -lglew32.dll -lglfw3dll -lfreetype -Wall -Wextra -pedantic -Ofast -pthread -fmax-errors=1

*/

#include "../physics/physics.hpp"
#include "../gui/graphics.hpp"

using namespace physics::literals;

int main()
{
	graphics window;

	// setting `dist` as the average intermolecular distance
	double dist = std::cbrt(physics::water_mass<> / physics::ice_density<>);

	unsigned n_side = 2; // 2x2x2 = 8 water molecules

	physics::molecular_system molsys;
	physics::fire integ(molsys, 0.1_fs);
	physics::thermodynamic_statistics stat_integ(integ);

	// set a simple cubic lattice as initial condition
	molsys.primitive_cubic_lattice(n_side, dist, physics::water_fba_eps<>);

	molsys.side *= 3;

	for (int i = 1; !window.should_close(); ++i)
	{
		stat_integ.step(0);

		std::wstringstream custom_text;

		custom_text << "SIMULATION (" << i << " steps)\n";

		custom_text << "Ewald parameter = " << molsys.lrsum.ewald_par();
		custom_text << L" Å^-1\nEstimated electrostatic force RMSE = " << molsys.lrsum.estimated_error_coulomb;
		custom_text << L" kcal/(mol Å)\nEstimated dispersion force RMSE = " << molsys.lrsum.estimated_error_lj;
		custom_text << L" kcal/(mol Å)\nEstimated total force RMSE = " << molsys.lrsum.estimated_error;

		// draw all steps in real time
		window.draw(molsys, custom_text.str());
	}

	return 0;
}





















