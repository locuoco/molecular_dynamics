//  Molecular dynamics example: critical point
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
#include <fstream>

/*

Compilation (MinGW):
g++ critical_point.cpp -o critical_point -std=c++20 -I C:\Users\aless\Desktop\myLib\include -L C:\Users\aless\Desktop\myLib\lib -lopengl32 -lglu32 -lglew32.dll -lglfw3dll -lfreetype -Wall -Wextra -pedantic -Ofast -pthread -fmax-errors=1

*/

#include "../../physics/physics.hpp"
#include "../../gui/graphics.hpp"

using namespace physics::literals;

int main()
{
	graphics window;

	window.camera_controls->pos = physics::vec3d(0, 0, 50);

	physics::molecular_system molsys;
	physics::leapfrog integ(molsys);

	// setting `dist` as the average intermolecular distance
	double dist = physics::lj_params<double>[int(physics::atom_type::AR)].half_Rmin * 2 / std::pow(2., 1./6);

	double dt = 16_fs;

	physics::vec3i n_side(3, 4, 3); // 3x4x3x4 = 144 argon atoms

	std::vector<double> init_temps, final_temps_mean, final_temps_sd, p_mean, p_sd, energy_mean, energy_sd;

	int num_temps = 41;
	double min_temp = 425, max_temp = 625;

	for (int i = 0; i < num_temps; ++i)
		init_temps.push_back(min_temp + i * (max_temp - min_temp) / std::max(num_temps-1, 1));

	int sim_i = 1;

	for (auto temp : init_temps)
	{
		physics::thermodynamic_statistics stat_integ(integ);

		molsys.lrsum.cell_list_multiplier({0, 2, 0});
		molsys.lrsum.cutoff_radius(36); // 36
		molsys.lrsum.limit_cutoff(false);

		molsys.temperature_ref = temp;

		// set a face-centered cubic lattice as initial condition
		molsys.face_centered_cubic_lattice(n_side, std::cbrt(4)*dist, physics::argon<>);

		molsys.side[1] *= 3;

		// run `equilibration_steps` steps before starting to calculate statistics
		int equilibration_steps = 4000;
		int simulation_steps = 30000;

		for (int i = -equilibration_steps; i < simulation_steps; ++i)
		{
			if (i == 0)
				continue;

			std::wstringstream custom_text;
			if (i <= 0)
			{
				integ.step(dt);
				custom_text << "EQUILIBRATION " << sim_i << " (other " << -i << " steps to go)\n";
			}
			else
			{
				stat_integ.step(dt);
				custom_text << "SIMULATION " << sim_i << " (" << i << " steps)\n";
			}

			custom_text << "Ewald parameter = " << molsys.lrsum.ewald_par();
			custom_text << L" Å^-1\nEstimated electrostatic force RMSE = " << molsys.lrsum.estimated_error_coulomb;
			custom_text << L" kcal/(mol Å)\nEstimated dispersion force RMSE = " << molsys.lrsum.estimated_error_lj;
			custom_text << L" kcal/(mol Å)\nEstimated total force RMSE = " << molsys.lrsum.estimated_error;
			custom_text << L" kcal/(mol Å)\n<P> = " << stat_integ.pressure_mean()/physics::atm<>
						<< " +/- " << stat_integ.pressure_sd()/physics::atm<>;

			custom_text << L" atm\n<ρ> = " << stat_integ.density_mean()/physics::kg_per_m3<>
						<< " +/- " << stat_integ.density_sd()/physics::kg_per_m3<>;

			custom_text << " kg/m^3\n<T> = " << stat_integ.temperature_mean() << " +/- " << stat_integ.temperature_sd();

			custom_text << " K\ntotal momentum magnitude = " << norm(molsys.p.sum());
			// total momentum should be around 0 for PPPM ik-differentiation, but increase as the square root of the
			// number of time steps (due to rounding errors). If ad-differentiation is used instead, total momentum
			// might drift randomly, since ad-differentiation does not conserve momentum.

			// draw all steps in "real time"
			window.draw(molsys, custom_text.str());
		}

		final_temps_mean.push_back(stat_integ.temperature_mean());
		final_temps_sd.push_back(stat_integ.temperature_sd());
		p_mean.push_back(stat_integ.pressure_mean());
		p_sd.push_back(stat_integ.pressure_sd());
		energy_mean.push_back(stat_integ.energy_mean());
		energy_sd.push_back(stat_integ.energy_sd());

		std::cout << sim_i << ": " << init_temps[sim_i-1] << ' ' << final_temps_mean[sim_i-1] << ' ' << final_temps_sd[sim_i-1] << '\n';

		++sim_i;
	}

	std::ofstream file("examples/phase_transitions/critical_point.txt");

	file << "T0 T dT P dP E dE\n";
	for (std::size_t i = 0; i < init_temps.size(); ++i)
		file << init_temps[i] << ' ' << final_temps_mean[i] << ' ' << final_temps_sd[i]
	                          << ' ' << p_mean[i] << ' ' << p_sd[i]
	                          << ' ' << energy_mean[i] << ' ' << energy_sd[i] << '\n';

	return 0;
}





















