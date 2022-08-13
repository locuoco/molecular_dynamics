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

#include <string>
#include <cmath> // sqrt

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

	unsigned n_side = 12; // 12x12x12 = 1728 water molecules

	physics::molecular_system<double, physics::pppm, physics::martyna_tobias_klein> molsys;

	// set a simple cubic lattice as initial condition
	molsys.primitive_cubic_lattice(n_side, dist, physics::water_fba_eps<>);

	double sumpressure = 0, sumdensity = 0, sumtemp = 0;
	double sumpressure2 = 0, sumdensity2 = 0, sumtemp2 = 0;
	int equilibration_steps = 5000;

	for (int i = -equilibration_steps; !window.should_close(); ++i)
	{
		if (i == 0)
			continue;

		std::stringstream custom_text;
		if (i <= 0)
		{
			custom_text << "EQUILIBRATION (other " << -i << " steps to go)\n";
			molsys.step(1e-3);
		}
		else
		{
			custom_text << "SIMULATION (" << i << " steps)\n";
			molsys.step(1e-3);

			// update statistics
			double pressure = molsys.pressure()/physics::atm<>;
			sumpressure += pressure;
			sumpressure2 += pressure*pressure;
			double density = molsys.density()/physics::kg_per_m3<>;
			sumdensity += density;
			sumdensity2 += density*density;
			double temp = molsys.temperature();
			sumtemp += temp;
			sumtemp2 += temp*temp;
		}

		custom_text << "Ewald parameter = " << molsys.lrsum.ewald_par();
		custom_text << " A^-1\nEstimated electrostatic force RMSE = " << molsys.lrsum.estimated_error_coulomb;
		custom_text << " kcal/(mol A)\nEstimated dispersion force RMSE = " << molsys.lrsum.estimated_error_lj;
		custom_text << " kcal/(mol A)\nEstimated total force RMSE = " << molsys.lrsum.estimated_error;
		double meanpressure = sumpressure/i;
		custom_text << " kcal/(mol A)\n<P> = " << meanpressure << " +/- " << std::sqrt((sumpressure2/i - meanpressure*meanpressure)/i);
		// average pressure should be 1 atm, although, due to the very large oscillations, the convergence
		// is very slow. These oscillations are due to water incompressibility (very small changes in configuration
		// space result into significant changes in instant pressure). Oscillations decrease as the number of atoms
		// increase.
		double meandensity = sumdensity/i;
		custom_text << " atm\n<rho> = " << meandensity << " +/- " << std::sqrt((sumdensity2/i - meandensity*meandensity)/i);
		// average density should be around 997 kg/m^3
		double meantemp = sumtemp/i;
		custom_text << " kg/m^3\n<T> = " << meantemp << " +/- " << std::sqrt((sumtemp2/i - meantemp*meantemp)/i);
		// average temperature should be 298.15 K
		custom_text << " K\ntotal momentum magnitude = " << norm(molsys.p.sum());
		// total momentum should be around 0 for PPPM ik-differentiation, but increase as the square root of the
		// number of time stpes (due to numerical errors). If ad-differentiation is used instead, total momentum
		// might drift randomly, since ad-differentiation does not conserve momentum.

		// draw all steps in real time
		window.draw(molsys, custom_text.str());
	}

	return 0;
}

// result from a run at 70318 steps:
//	<P> = 1.02049 +/- 2.69440 atm
//	<rho> = 996.637 +/- 0.021 kg/m^3
// 	<T> = 298.169 +/- 0.013 K





















