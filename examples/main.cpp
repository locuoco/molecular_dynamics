//  Molecular dynamics example
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
#include <random>

/*

g++ examples/main.cpp -o examples/mold -std=c++20 -I C:\Users\aless\Desktop\myLib\include -L C:\Users\aless\Desktop\myLib\lib -lopengl32 -lglu32 -lglew32.dll -lglfw3dll -lfreetype -Wall -Wextra -pedantic -Ofast -pthread -fmax-errors=1

*/

#include "../physics/physics.hpp"
#include "../gui/graphics.hpp"

int main(const int, const char**)
{
	graphics window;

	int side = 16, hside = side/2;
	int Nmol = side*side*side;
	double volume = (Nmol*18.0154)/0.60043754468329832909213892745879;
	// 0.5522 = density of ice in AKMA units
	// 0.6022 = density of water at 4 °C
	// 0.6004 = density of water at 25 °C
	// 0.5901 = density of TIP3P water model
	double dist = std::cbrt(volume)/side;
	//dist = 7;
	double init_temp = 298.15, end_temp = 220;
	std::cout << "Nmol = " << Nmol << ", dist = " << dist << std::endl;

	physics::molecular_system<double, physics::pppm, physics::nose_hoover> molsys(
		dist*side,
		init_temp
	);

	std::mt19937_64 mersenne_twister(0);
	std::uniform_real_distribution<double> u_dist(0, 1);
	//double angles[3];
	//physics::mat3<double> rot;

	for (int i = 0; i < side; ++i)
		for (int j = 0; j < side; ++j)
			for (int k = 0; k < side; ++k)
			{
				physics::molecule w = physics::water_tip3p_lr<>;
				physics::vec3d p = {(i-hside)*dist, (j-hside)*dist, (k-hside)*dist};
				for (int l = 0; l < 3; ++l)
					p[l] += (u_dist(mersenne_twister) - 0.5) * 0.1;
				/*for (int l = 0; l < 3; ++l)
					angles[l] = u_dist(mersenne_twister)*math::two_pi<double>();
				rot = physics::rotation_yaw_pitch_roll(angles[0], angles[1], angles[2]);
				for (int l = 0; l < 3; ++l)
					w.x[l] = rot % w.x[l];*/
				molsys.add_molecule(w, p);
			}

	//molsys.rescale_temperature = true;

	while (!window.should_close())
	{
		//molsys.temperature_fixed = end_temp + (molsys.temperature_fixed - end_temp) * 0.9999;
		molsys.step(1e-3);
		window.draw(molsys);
	}

	return 0;
}





















