//  Molecular dynamics
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

g++ main.cpp -o pmd -std=c++20 -I C:\Users\aless\Desktop\myLib\include -L C:\Users\aless\Desktop\myLib\lib -lopengl32 -lglu32 -lglew32.dll -lglfw3dll -lfreetype -Wall -Wextra -pedantic -Ofast -fmax-errors=1

*/

#include "physics/physics.hpp"
#include "graphics.hpp"

int main(const int, const char**)
{
	graphics window;

	//test();

	int side = 6, hside = side/2;
	//double volume = (side*side*side*18)/0.55; // 0.55 = density of ice in AKMA units
	//double dist = std::cbrt(volume)/side;
	double dist = 12;
	double init_temp = 1000, end_temp = 220;
	std::cout << "dist = " << dist << std::endl;

	physics::molecular_system<double, physics::forest_ruth, physics::leapfrog> molsys(
		dist*side,
		init_temp,
		physics::DW5<>
	);

	/*physics::molecule pert_water = physics::water_tip3p_lr<>;
	for (unsigned int i = 0; i < pert_water.n; ++i)
		pert_water.x[i][1] = -pert_water.x[i][1];*/

	for (int i = 0; i <= side; ++i)
		for (int j = 0; j <= side; ++j)
			for (int k = 0; k <= side; ++k)
				molsys.add_molecule(physics::water_fba_eps<>, {(i-hside)*dist, (j-hside)*dist, (k-hside)*dist});

	/*int side = 10, hside = side/2;
	for (int i = 0; i <= side; ++i)
		for (int j = 0; j <= side; ++j)
			molsys.add_molecule(physics::water_tip3p_lr<>, {(i-hside)*dist, (j-hside)*dist, 0});*/

	for (int i = 0; i < 100; ++i)
		molsys.step(1e-3);

	while (!window.should_close())
	{
		molsys.temperature = end_temp + (molsys.temperature - end_temp) * 0.9998;
		for (int i = 0; i < 2; ++i)
			molsys.step(2e-3, true);

		window.draw(molsys);
	}

	return 0;
}





















