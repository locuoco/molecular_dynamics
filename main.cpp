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
#include <random>

/*

g++ main.cpp -o pmd -std=c++20 -I C:\Users\aless\Desktop\myLib\include -L C:\Users\aless\Desktop\myLib\lib -lopengl32 -lglu32 -lglew32.dll -lglfw3dll -lfreetype -Wall -Wextra -pedantic -Ofast -pthread -fmax-errors=1

*/

#include "physics/physics.hpp"
#include "graphics.hpp"

int main(const int, const char**)
{
	graphics window;

	//test();

	int side = 8, hside = side/2;
	int Nmol = side*side*side;
	double volume = (Nmol*18.0154)/0.602214076; // 0.5522 = density of ice in AKMA units
	double dist = std::cbrt(volume)/side;
	//dist = 7;
	double init_temp = 298.15, end_temp = 220;
	std::cout << "Nmol = " << Nmol << ", dist = " << dist << std::endl;

	physics::molecular_system<double, physics::nose_hoover> molsys(
		dist*side,
		init_temp,
		physics::DW5<>
	);

	std::mt19937_64 oliver_twist(0);
	std::uniform_real_distribution<double> u_dist(0, 1);
	double angles[3];
	physics::mat3<double> rot;

	for (int i = 0; i < side; ++i)
		for (int j = 0; j < side; ++j)
			for (int k = 0; k < side; ++k)
			{
				physics::molecule w = physics::water_tip3p<>;
				physics::point3<double> p = {(i-hside)*dist, (j-hside)*dist, (k-hside)*dist};
				for (int l = 0; l < 3; ++l)
					p[l] += (u_dist(oliver_twist) - 0.5) * (dist - 2);
				for (int l = 0; l < 3; ++l)
					angles[l] = u_dist(oliver_twist)*math::two_pi<double>();
				rot = physics::rotation_yaw_pitch_roll(angles[0], angles[1], angles[2]);
				for (int l = 0; l < 3; ++l)
					w.x[l] = rot % w.x[l];
				molsys.add_molecule(w, p);
			}

	//molsys.rescale_temperature = true;

	while (!window.should_close())
	{
		//molsys.temperature_fixed = end_temp + (molsys.temperature_fixed - end_temp) * 0.9999;
		for (int i = 0; i < 2; ++i)
			molsys.step(1e-3);
		//std::this_thread::sleep_for(std::chrono::milliseconds(100));
		window.draw(molsys);
	}

	return 0;
}





















