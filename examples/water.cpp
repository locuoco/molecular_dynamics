//  Molecular dynamics example: water molecules
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

/*

Compilation:
g++ water.cpp -o water -std=c++20 -I C:\Users\aless\Desktop\myLib\include -L C:\Users\aless\Desktop\myLib\lib -lopengl32 -lglu32 -lglew32.dll -lglfw3dll -lfreetype -Wall -Wextra -pedantic -Ofast -pthread -fmax-errors=1

*/

#include "../physics/physics.hpp"
#include "../gui/graphics.hpp"

int main()
{
	graphics window;

	double dist = std::cbrt(physics::water_mass<> / physics::water_density25<>);

	unsigned n_side = 16;

	physics::molecular_system<double, physics::pppm, physics::nose_hoover> molsys;

	// set a simple cubic lattice as initial condition
	molsys.primitive_cubic_lattice(n_side, dist, physics::water_tip3p_lr<>);

	while (!window.should_close())
	{
		molsys.step(1e-3);
		window.draw(molsys);
	}

	return 0;
}





















