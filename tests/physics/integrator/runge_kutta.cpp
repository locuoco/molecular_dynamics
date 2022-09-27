//  Runge-Kutta methods tests
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
#include <cmath> // abs

/*

Compilation:
g++ runge_kutta.cpp -o runge_kutta -std=c++20 -Wall -Wextra -pedantic -Ofast -fmax-errors=1

*/

#include "../periodic_harmonic_oscillator.hpp"
#include "../../../physics/integrator/runge_kutta.hpp"

template <template<typename> typename ERK>
void test_explicit_runge_kutta()
// Test that the values of the parameters of an explicit Runge-Kutta method are correct.
// `ERK` is the template class corresponding to the explicit Runge-Kutta method.
// The test passed if the sum of the space parameters of each row corresponds
// to the time parameter, and if the sum of parameters in the last line corresponds
// to 1 (up to rounding errors).
{
	using erk = ERK<periodic_harmonic_oscillator>;

	periodic_harmonic_oscillator sys;
	erk integ(sys);
	typename erk::scalar_type sumpars;

	for (std::size_t i = 0; i < erk::num_stages-1; ++i)
	{
		sumpars = 0;
		for (std::size_t j = 1; j < erk::num_stages; ++j)
			sumpars += integ.pars[i*erk::num_stages + j];
		assert(abs(integ.pars[i*erk::num_stages] - sumpars) < 1e-15);
	}
	sumpars = 0;
	for (std::size_t j = 0; j < erk::num_stages; ++j)
		sumpars += integ.pars[erk::num_stages*(erk::num_stages-1) + j];
	assert(abs(sumpars - 1) < 1e-15);
}

int main()
{
	test_explicit_runge_kutta<physics::euler>();
	test_explicit_runge_kutta<physics::midpoint>();
	test_explicit_runge_kutta<physics::heun2>();
	test_explicit_runge_kutta<physics::ralston2>();
	test_explicit_runge_kutta<physics::rk4>();
	test_explicit_runge_kutta<physics::rk4_3_8>();
	test_explicit_runge_kutta<physics::ralston4>();
	test_explicit_runge_kutta<physics::butcher6>();
	test_explicit_runge_kutta<physics::verner8>();

	std::cout << "All tests passed successfully!" << std::endl;

	return 0;
}
























