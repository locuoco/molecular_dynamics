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

#include <cmath>
#include <sstream>

/*

Compilation (MinGW):
g++ water.cpp -o water -std=c++20 -I C:\Users\aless\Desktop\myLib\include -L C:\Users\aless\Desktop\myLib\lib -lopengl32 -lglu32 -lglew32.dll -lglfw3dll -lfreetype -Wall -Wextra -pedantic -Ofast -pthread -fmax-errors=1

*/

#include "../physics/physics.hpp"
#include "../gui/graphics.hpp"

using namespace physics::literals;

int main()
{
	graphics window;

	// setting `dist` as the average intermolecular distance
	double dist = std::cbrt(physics::water_mass<> / physics::water_density25<>);

	unsigned n_side = 12; // 12x12x12 = 1728 water molecules

	physics::molecular_system molsys;
	physics::mtk integ(molsys);
	physics::thermodynamic_statistics stat_integ(integ);

	// set a simple cubic lattice as initial condition
	molsys.primitive_cubic_lattice(n_side, dist, physics::water_fba_eps<>);

	// run `equilibration_steps` steps before starting to calculate statistics
	int equilibration_steps = 5000;

	for (int i = -equilibration_steps; !window.should_close(); ++i)
	{
		if (i == 0)
			continue;

		std::wstringstream custom_text;
		if (i <= 0)
		{
			integ.step(1_fs);
			custom_text << "EQUILIBRATION (other " << -i << " steps to go)\n";
		}
		else
		{
			stat_integ.step(1_fs);
			custom_text << "SIMULATION (" << i << " steps)\n";
		}

		custom_text << "Ewald parameter = " << molsys.lrsum.ewald_par();
		custom_text << L" Å^-1\nEstimated electrostatic force RMSE = " << molsys.lrsum.estimated_error_coulomb;
		custom_text << L" kcal/(mol Å)\nEstimated dispersion force RMSE = " << molsys.lrsum.estimated_error_lj;
		custom_text << L" kcal/(mol Å)\nEstimated total force RMSE = " << molsys.lrsum.estimated_error;
		custom_text << L" kcal/(mol Å)\n<P> = " << stat_integ.pressure_mean()/physics::atm<>
		            << " +/- " << stat_integ.pressure_sd()/physics::atm<>;
		// average pressure should be 1 atm, although, due to the very large oscillations, the convergence
		// is very slow. These oscillations are due to water incompressibility (very small changes in configuration
		// space result into significant changes in instant pressure). Oscillations decrease as the number of atoms
		// increase.
		custom_text << L" atm\n<ρ> = " << stat_integ.density_mean()/physics::kg_per_m3<>
		            << " +/- " << stat_integ.density_sd()/physics::kg_per_m3<>;
		// average density should be around 997 kg/m^3
		custom_text << " kg/m^3\n<T> = " << stat_integ.temperature_mean() << " +/- " << stat_integ.temperature_sd();
		// average temperature should be 298.15 K
		custom_text << " K\ntotal momentum magnitude = " << norm(molsys.p.sum());
		// total momentum should be around 0 for PPPM ik-differentiation, but increase as the square root of the
		// number of time steps (due to rounding errors). If ad-differentiation is used instead, total momentum
		// might drift randomly, since ad-differentiation does not conserve momentum.

		// draw all steps in real time
		window.draw(molsys, custom_text.str());
	}

	return 0;
}

// result from a run at 70318 steps:
//	<P> = 1.02049 +/- 2.69440 atm
//	<ρ> = 996.637 +/- 0.021 kg/m^3
// 	<T> = 298.169 +/- 0.013 K





















