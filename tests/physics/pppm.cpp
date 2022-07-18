//  Particle-particle, particle-mesh method tests
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

#include <iostream> // cout, endl
#include <cassert>
#include <cmath> // fabs
#include <random> // mt19937, uniform_real_distribution

/*

g++ pppm.cpp -o pppm -std=c++20 -Wall -Wextra -pedantic -Ofast -pthread -fmax-errors=1

*/

#include "../../physics/physics.hpp"

std::mt19937 mersenne_twister;

template <typename T>
std::uniform_real_distribution<T> dist(0, 1);

void test_charge_assignment_function()
// test property of charge assignment function, i.e.:
// sum_p W_p (x) = 1 for all x (charge must be conserved)
{
	using std::size_t;
	double x = dist<double>(mersenne_twister) - 0.5, sum;

	for (size_t order = physics::pppm_min_order; order <= physics::pppm_max_order; ++order)
	{
		sum = 0;
		for (size_t k = 0; k < order; ++k)
			sum += physics::charge_assignment_function(x, k, order);

		std::cout << "order " << order << ": "; 
		assert(std::fabs(sum - 1) < 1.e-14);
		std::cout << "ok\n";
	}
}

void test_force_accuracy()
{
	using std::sqrt;

	int side = 15, hside = side/2;
	int Nmol = side*side*side;
	double volume = (Nmol*18.0154)/0.602214076;
	double dist = std::cbrt(volume)/side;

	physics::molecular_system<double, physics::pppm, physics::leapfrog> sys(dist*side);
	physics::molecular_system<double, physics::ewald, physics::leapfrog> sys_ref;

	std::mt19937_64 mersenne_twister(0);
	std::uniform_real_distribution<double> u_dist(0, 1);

	sys.lrsum.cutoff_radius(9);
	sys.lrsum.charge_assignment_order(7);
	sys.lrsum.set_diff_scheme("ik");
	sys.lrsum.cell_multiplier(0);

	for (int i = 0; i < side; ++i)
		for (int j = 0; j < side; ++j)
			for (int k = 0; k < side; ++k)
			{
				physics::molecule w = physics::water_tip3p<>;
				physics::point3<double> p {(i-hside)*dist, (j-hside)*dist, (k-hside)*dist};
				for (int l = 0; l < 3; ++l)
					p[l] += (u_dist(mersenne_twister) - 0.5) * 0.1;
				sys.add_molecule(w, p);
			}

	sys_ref = sys;

	sys.step(0);
	sys_ref.step(0);

	auto sqr = [](physics::point3d x) { return x*x; };
	auto sum_coords = [](physics::point3d x) { return x[0]+x[1]+x[2]; };

	double rmse = sqrt(sum_coords((sys.f - sys_ref.f).apply(sqr).sum()) / sys.n);
	double rms_elec = sqrt(sum_coords(sys_ref.lrsum.f.apply(sqr).sum()) / sys.n);
	double rms_tot = sqrt(sum_coords(sys_ref.f.apply(sqr).sum()) / sys.n);

	std::cout << "Calculated total RMS error = " << rmse << '\n';
	std::cout << "Force electrostatic RMS = " << rms_elec << '\n';
	std::cout << "Total force RMS (excluded intra-molecular pairs) = " << rms_tot << '\n';
}

int main()
{
	test_charge_assignment_function();
	test_force_accuracy();

	return 0;
}
























