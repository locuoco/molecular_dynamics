

#include <iostream>

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
	double volume = (side*side*side*18)/0.55; // 0.55 = density of ice in AKMA units
	
	double dist = std::cbrt(volume)/side;
	std::cout << "dist = " << dist << std::endl;
	
	physics::molecular_system<double, physics::suzuki8<double>> molsys(dist*side, 252, physics::DW5<double>);

	physics::molecule pert_water = physics::water_fba_eps<double>;
	for (unsigned int i = 0; i < pert_water.n; ++i)
		pert_water.x[i][1] = -pert_water.x[i][1];

	for (int i = 0; i <= side; ++i)
		for (int j = 0; j <= side; ++j)
			for (int k = 0; k <= side; ++k)
				molsys.add_molecule(physics::water_fba_eps<double>, {(i-hside)*dist, (j-hside)*dist, (k-hside)*dist});

	/*int side = 10, hside = side/2;
	for (int i = 0; i <= side; ++i)
		for (int j = 0; j <= side; ++j)
			molsys.add_molecule(physics::water_fba_eps<double>, {(i-hside)*dist, (j-hside)*dist, 0});*/

	float *atomPosType = new float[molsys.max_atoms * 4];
	// x, y, z, type
	while (!window.should_close())
	{
		for (int i = 0; i < 1; ++i)
			molsys.step(1e-3, true);
		for (int i = 0; i < (int)molsys.n; ++i)
		{
			atomPosType[i*4+0] = molsys.x[i][0];
			atomPosType[i*4+1] = molsys.x[i][1];
			atomPosType[i*4+2] = molsys.x[i][2];
			atomPosType[i*4+3] = physics::atom_number[int(molsys.id[i])];
		}

		window.draw(atomPosType, molsys.n);
	}

	return 0;
}





















