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

#include "../../physics/pppm.hpp"

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
		assert(std::fabs(sum - 1) < 1.e-12);
		std::cout << "ok\n";
	}
}

int main()
{
	test_charge_assignment_function();

	return 0;
}
























