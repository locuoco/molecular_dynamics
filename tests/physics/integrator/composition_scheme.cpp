//  Composition schemes tests
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
g++ composition_scheme.cpp -o composition_scheme -std=c++20 -Wall -Wextra -pedantic -Ofast -fmax-errors=1

*/

#include "../periodic_harmonic_oscillator.hpp"
#include "../../../physics/integrator/composition_scheme.hpp"

template <template<typename> typename CSInteg>
void test_composition_scheme()
// Test that the values of the parameters of a composition scheme are correct.
// `CSInteg` is the template class corresponding to the composition scheme.
// The test passed if the sum of the parameters is 1.
{
	using cs = CSInteg<periodic_harmonic_oscillator>;

	cs integ;
	typename cs::scalar_type sumpars = 0;
	for (std::size_t i = 0; i < cs::num_stages; ++i)
		sumpars += integ.pars[std::min(i, cs::num_stages-i-1)];

	assert(abs(sumpars - 1) < 1e-15);
}

int main()
{
	test_composition_scheme<physics::forest_ruth>();
	test_composition_scheme<physics::suzuki4>();
	test_composition_scheme<physics::kahan_li4a>();
	test_composition_scheme<physics::kahan_li4b>();
	test_composition_scheme<physics::yoshida6>();
	test_composition_scheme<physics::kahan_li6a>();
	test_composition_scheme<physics::kahan_li6b>();
	test_composition_scheme<physics::suzuki8>();
	test_composition_scheme<physics::kahan_li8a>();
	test_composition_scheme<physics::kahan_li8b>();
	test_composition_scheme<physics::kahan_li10a>();
	test_composition_scheme<physics::kahan_li10b>();
	test_composition_scheme<physics::kahan_li10c>();
	test_composition_scheme<physics::kahan_li10d>();

	std::cout << "All tests passed successfully!" << std::endl;

	return 0;
}
























