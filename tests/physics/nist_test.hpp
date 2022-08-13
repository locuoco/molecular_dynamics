//  NIST test sample loader
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

#ifndef TESTS_PHYSICS_NIST_TEST
#define TESTS_PHYSICS_NIST_TEST

#include <fstream>

#include "../../physics/molecule.hpp"

template <typename System>
void load_nist_sample(System& sys)
// Load NIST sample for Ewald summation reference calculation
// Samples and data taken from:
// https://www.nist.gov/mml/csd/chemical-informatics-group/spce-water-reference-calculations-10a-cutoff
// `s` is the physical system where the data is loaded to.
{
	sys.clear();

	std::ifstream input("spce_sample_config_periodic1.txt");

	double dummy;
	char dummy_char;

	input >> sys.side;
	input >> dummy;
	input >> dummy;
	input >> dummy;
	input >> dummy;
	while (input)
	{
		physics::molecule w = physics::water_spc_e<typename System::scalar_type>;
		input >> w.x[1][0];
		input >> w.x[1][1];
		input >> w.x[1][2];
		input >> dummy_char;
		input >> dummy;
		input >> w.x[0][0];
		input >> w.x[0][1];
		input >> w.x[0][2];
		input >> dummy_char;
		input >> dummy;
		input >> w.x[2][0];
		input >> w.x[2][1];
		input >> w.x[2][2];
		input >> dummy_char;
		input >> dummy;
		w.x[0] = w.x[1] + remainder(w.x[0] - w.x[1], sys.side);
		w.x[2] = w.x[1] + remainder(w.x[2] - w.x[1], sys.side);
		sys.add_molecule(w);
	}
	sys.fetch();
}

#endif // TESTS_PHYSICS_NIST_TEST








