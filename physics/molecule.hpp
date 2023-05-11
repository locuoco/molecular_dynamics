//  Molecules
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

#ifndef PHYSICS_MOLECULE_H
#define PHYSICS_MOLECULE_H

#include <concepts> // floating_point
#include <valarray>
#include <vector>
#include <map>
#include <tuple>
#include <algorithm> // min, max, generate
#include <cmath> // sqrt, atan2
#include <random> // mt19937_64, normal_distribution
#include <stdexcept> // runtime_error
#include <numbers> // numbers::pi_v

#include "../math/helper.hpp" // deg2rad, two1_6, is_even

#include "tensor.hpp"
#include "pppm.hpp"
#include "../utils/thread_pool.hpp"

namespace physics
{
	// Coulomb constant in AKMA units
	template <std::floating_point T = double>
	inline constexpr T kC = 332.06371330062682040436798930743L;

	// square root of the Coulomb constant in AKMA units
	template <std::floating_point T = double>
	inline constexpr T sqrtkC = 18.222615435239443706926657644703L;

	// Boltzmann constant in AKMA units (kcal/(mol K))
	template <std::floating_point T = double>
	inline constexpr T kB = 0.00198720425864083173996175908222L;

	// Avogadro number (mol^-1)
	template <std::floating_point T = double>
	inline constexpr T NA = 6.02214076e23L;

	// 1 atmosphere in AKMA units
	template <std::floating_point T = double>
	inline constexpr T atm = 1.458397257425908222e-5L;

	// 1 kg/m^3 in AKMA units
	template <std::floating_point T = double>
	inline constexpr T kg_per_m3 = 6.02214076e-4L;

	// 1 g/cm^3 in AKMA units
	template <std::floating_point T = double>
	inline constexpr T g_per_cm3 = 6.02214076e-1L;

	// 1 picosecond (1 ps = 10^-12 s) in AKMA units
	template <std::floating_point T = double>
	inline constexpr T picosecond = 20.4548282844073286665866518779L;

	// 1 femtosecond (1 fs = 10^-15 s) in AKMA units
	template <std::floating_point T = double>
	inline constexpr T femtosecond = 0.0204548282844073286665866518779L;

	// water self-diffusion constant at 5 °C in AKMA units
	template <std::floating_point T = double>
	inline constexpr T DW5 = 0.006164803646683544007L;

	// water self-diffusion constant at 25 °C in AKMA units
	template <std::floating_point T = double>
	inline constexpr T DW25 = 0.01123940014569822972L;

	// density of ice in AKMA units (a.m.u. / A^3)
	template <std::floating_point T = double>
	inline constexpr T ice_density = 917*kg_per_m3<long double>;

	// density of water at 4 °C in AKMA units (a.m.u. / A^3)
	template <std::floating_point T = double>
	inline constexpr T water_density4 = 1000*kg_per_m3<long double>;

	// density of water at 25 °C in AKMA units (a.m.u. / A^3)
	template <std::floating_point T = double>
	inline constexpr T water_density25 = 997.13L*kg_per_m3<long double>;

	// water molecular mass (in atomic mass units)
	template <std::floating_point T = double>
	inline constexpr T water_mass = 18.0154L;

	// kcal to kJ conversion factor
	template <std::floating_point T = double>
	inline constexpr T kcal2kj = 4.184L;

	// NaCl Madelung constant
	template <std::floating_point T = double>
	inline constexpr T nacl_madelung = 1.74756459463318L;

	// CsCl Madelung constant
	template <std::floating_point T = double>
	inline constexpr T cscl_madelung = 1.76267477307098L;

	// ZnS Madelung constant
	template <std::floating_point T = double>
	inline constexpr T zns_madelung = 1.6380550533L;

	// NaCl lattice constant
	template <std::floating_point T = double>
	inline constexpr T nacl_lattice = 5.6402L;

	// CsCl lattice constant
	template <std::floating_point T = double>
	inline constexpr T cscl_lattice = 4.119L;

	namespace literals
	{
		// literals suffixes to convert to AKMA units

		// picoseconds to AKMA units
		constexpr long double operator""_ps(long double x)
		{
			return x * picosecond<long double>;
		}
		constexpr long double operator""_ps(unsigned long long x)
		{
			return x * picosecond<long double>;
		}

		// femtoseconds to AKMA units
		constexpr long double operator""_fs(long double x)
		{
			return x * femtosecond<long double>;
		}
		constexpr long double operator""_fs(unsigned long long x)
		{
			return x * femtosecond<long double>;
		}

		// atmospheres to AKMA units
		constexpr long double operator""_atm(long double x)
		{
			return x * atm<long double>;
		}
		constexpr long double operator""_atm(unsigned long long x)
		{
			return x * atm<long double>;
		}
	}

	enum class atom_type : unsigned char
	{
		HT,  // TIP3P water hydrogen atom
		HTL, // TIP3P-Ewald water hydrogen atom
		HF,  // FBA/eps water hydrogen atom
		HSE, // SPC/E water hydrogen atom
		OT,  // TIP3P water oxygen atom
		OTL, // TIP3P-Ewald water oxygen atom
		OF,  // FBA/eps water oxygen atom
		OSE, // SPC/E water oxygen atom
		SOD, // sodium (Na+) cation
		CLA, // chloride (Cl-) anion
		AR,  // argon (Ar) atom
		CES, // caesium (Cs+) cation
		N    // number of atom types
	};
	inline constexpr std::size_t num_atom_types = std::size_t(atom_type::N);

	// maximum bond coordination number assumed (not more than `max_bonds` bonds per atoms are allowed)
	inline constexpr std::size_t max_bonds = 4;

	// maximum dihedral multiplicity
	inline constexpr std::size_t max_dihedrals = 6;

	// maximum number of atoms
	inline constexpr std::size_t max_atoms = 2'000'000;

	inline constexpr unsigned char atom_number[num_atom_types] {1, 1, 1, 1, 8, 8, 8, 8, 11, 17, 18, 55};

	// atom masses in atomic units
	// note that for ions, the electron mass is not added/subtracted, as it is not usually significant
	template <std::floating_point T>
	inline constexpr T atom_mass[num_atom_types]
	{
		  1.008L,   // HT
		  1.008L,   // HTL
		  1.008L,   // HF
		  1.008L,   // HSE
		 15.9994L,  // OT
		 15.9994L,  // OTL
		 15.9994L,  // OF
		 15.9994L,  // OSE
		 22.98977L, // SOD
		 35.45L,    // CLA
		 39.948L,   // AR
		132.90545L, // CES
	};

	template <std::floating_point T>
	struct bond_struct { T Kb, b0; };

	template <std::floating_point T>
	struct angle_struct { T Ktheta, theta0, Kub, S0; };

	template <std::floating_point T>
	struct dihedral_struct { T Kchi[max_dihedrals], delta[max_dihedrals]; };

	template <std::floating_point T>
	struct improper_struct { T Kpsi, psi0; };

	template <std::floating_point T>
	struct lj_struct { T epsilon, half_Rmin; };

	template <std::floating_point T>
	inline std::map<std::pair<atom_type, atom_type>, bond_struct<T>> bond_params =
	// 1-2 bonds coefficients
	// (Kb, b0) parameters associated to a pair of adjacent atoms
	{
		{ {atom_type::HT,  atom_type::OT }, {450,       0.9572L} },
		{ {atom_type::HTL, atom_type::OTL}, {450,       0.9572L} },
		{ {atom_type::HF,  atom_type::OF }, {358.5086L, 1.027L } },
		{ {atom_type::HSE, atom_type::OSE}, {0,         1      } },
	};

	template <std::floating_point T>
	inline std::map<std::tuple<atom_type, atom_type, atom_type>, angle_struct<T>> angle_params =
	// 1-3 bonds coefficients
	// (Ktheta, theta0, Kub, S0) parameters associated to a triplet of adjacent atoms
	// Kub, S0 are the Urey-Bradley parameters
	{
		{ {atom_type::HT,  atom_type::OT,  atom_type::HT }, {55,       math::deg2rad(104.52L), 0, 0} },
		{ {atom_type::HTL, atom_type::OTL, atom_type::HTL}, {55,       math::deg2rad(104.52L), 0, 0} },
		{ {atom_type::HF,  atom_type::OF,  atom_type::HF }, {45.7696L, math::deg2rad(114.70L), 0, 0} },
		{ {atom_type::HSE, atom_type::OSE, atom_type::HSE}, {0,        math::deg2rad(109.47L), 0, 0} },
	};

	template <std::floating_point T>
	inline lj_struct<T> lj_params[num_atom_types] =
	// Lennard-Jones parameters
	// (epsilon, R_min/2) associated to each atom types
	// taken from CHARMM force field
	// R_min is given by sigma * 2^(1/6)
	// for ions, see also Beglov, Roux, 1994 and Stote, Karplus, 1995
	{
		{0.046L,         0.2245L   }, // HT
		{0,              0         }, // HTL
		{0,              0         }, // HF
		{0,              0         }, // HSE
		{0.1521L,        1.7682L   }, // OT
		{0.102L,         1.7892L   }, // OTL
		{0.1894L,        1.78337L  }, // OF
		{0.15539426812L, 1.7766093L}, // OSE
		{0.0469L,        1.41075L  }, // SOD
		{0.150L,         2.27L     }, // CLA
		{0.22165623256L, 2.03334L  }, // AR
		{0.190L,         2.100L    }, // CES
	};

	template <std::floating_point T = double>
	struct molecule
	{
		state<T, 3>                          x{{0}};      // positions of the atoms
		std::vector<atom_type>               id;          // identities of the atoms
		std::vector<T>                       q{0};        // partial charges
		std::vector<fixed_list<max_bonds>>   bonds{{}};   // bonds list as an adjacency list for each atom
		std::vector<std::array<unsigned, 4>> impropers{}; // list of impropers in the whole molecule
		unsigned                             n = 1;       // number of atoms in the molecule
	};

	template <typename T>
	T mass_of(const molecule<T>& species)
	// Return the total mass of a chemical species.
	{
		T mass = 0;
		for (auto id : species.id)
			mass += atom_mass<T>[ int(id) ];
		return mass;
	}

	template <typename T>
	T charge_of(const molecule<T>& species)
	// Return the total charge of a chemical species.
	{
		T charge = 0;
		for (auto q : species.q)
			charge += q;
		return charge;
	}

	template <std::floating_point T = double>
	inline const molecule<T> empty_molecule
	// Empty molecule
	{
		.x     = {},
		.id    = {},
		.q     = {},
		.bonds = {},
		.n     = 0
	};

	template <std::floating_point T = double>
	inline const molecule<T> water_tip3p
	// TIP3P water model
	// see also W. L. Jorgensen, J. Chandrasekhar, J. D. Madura, R. W. Impey,
	// M. L. Klein, "Comparison of simple potential functions for
	// simulating liquid water", J. Chem. Phys. 79 926-935 (1983).
	// This variant uses flexible bonds, which tend to overestimate the density.
	// This is the model used in CHARMM.
	{
		.x     = { {-0.75695L, 0.58588L}, {0}, {0.75695L, 0.58588L} }, // oxygen at the origin
		.id    = {atom_type::HT, atom_type::OT, atom_type::HT},
		.q     = {0.417L, -0.834L, 0.417L},
		.bonds = { {1}, {0, 2}, {1} },
		.n     = 3
	};

	template <std::floating_point T = double>
	inline const molecule<T> water_tip3p_lr
	// TIP3P water model optimized for long-range interactions
	// see also D. J. Price, C. L. Brooks, "A modified TIP3P water potential for simulation with Ewald summation", 2004
	// This variant uses flexible bonds, which tend to overestimate the density.
	{
		.x     = { {-0.75695L, 0.58588L}, {0}, {0.75695L, 0.58588L} }, // oxygen at the origin
		.id    = {atom_type::HTL, atom_type::OTL, atom_type::HTL},
		.q     = {0.415L, -0.830L, 0.415L},
		.bonds = { {1}, {0, 2}, {1} },
		.n     = 3
	};

	template <std::floating_point T = double>
	inline const molecule<T> water_fba_eps
	// FBA/epsilon water model
	// see R. Fuentes-Azcatl, M. Barbosa, "Flexible Bond and Angle, FBA/epsilon model of water", 2018
	// This model gives an accurate estimate of water density at 25 °C using Ewald summation.
	// FBA/epsilon model is the most accurate flexible water model currently implemented in this library.
	{
		.x     = { {-0.75695L, 0.58588L}, {0}, {0.75695L, 0.58588L} }, // oxygen at the origin
		.id    = {atom_type::HF, atom_type::OF, atom_type::HF},
		.q     = {0.4225L, -0.8450L, 0.4225L},
		.bonds = { {1}, {0, 2}, {1} },
		.n     = 3
	};

	template <std::floating_point T = double>
	inline const molecule<T> water_spc_e
	// SPC/E water model
	// only used for testing Ewald summation, do not use for simulations! (fixed bonds not supported)
	{
		.x     = { {-0.75695L, 0.58588L}, {0}, {0.75695L, 0.58588L} }, // oxygen at the origin
		.id    = {atom_type::HSE, atom_type::OSE, atom_type::HSE},
		.q     = {0.4238L, -0.8476L, 0.4238L},
		.bonds = { {1}, {0, 2}, {1} },
		.n     = 3
	};

	template <std::floating_point T = double>
	inline const molecule<T> sodium_ion
	// Sodium ion
	{
		.id = {atom_type::SOD},
		.q  = {1},
	};

	template <std::floating_point T = double>
	inline const molecule<T> chloride_ion
	// Chloride ion
	{
		.id = {atom_type::CLA},
		.q  = {-1},
	};

	template <std::floating_point T = double>
	inline const molecule<T> caesium_ion
	// Caesium ion
	{
		.id = {atom_type::CES},
		.q  = {1},
	};

	template <std::floating_point T = double>
	inline const molecule<T> argon
	// Argon atom
	{
		.id = {atom_type::AR},
	};

	template <
		std::floating_point T = double,
		template <typename, typename> typename LRSummationT = pppm
	>
	class molecular_system : public physical_system_base<T, state<T, 3>>
	{
		using LRSum = LRSummationT<T, state<T, 3>>;

		std::normal_distribution<T> n_dist = std::normal_distribution<T>(0, 1);
		T kin; // kinetic energy

		public:

		state<T, 3> x, p, v, f; // position, momentum, velocity, force
		state<T, 3> noise;      // gaussian noise
		std::valarray<T> m, gamma; // mass, damping factor
		std::vector<vec3<T>> x_tmp, p_tmp, v_tmp, f_tmp; // used for initialization
		std::vector<T> m_tmp; // used for initialization
		std::vector<T> z, lj_sqrteps, lj_halfR; // normalized partial charges + Lennard-Jones parameters
		std::vector<atom_type> id; // identity of the atom
		std::vector<fixed_list<max_bonds>> bonds; // bonds as an adjacency list
		std::vector<std::array<unsigned, 4>> impropers; // list of impropers in the whole system
		utils::thread_pool tp; // thread pool
		LRSum lrsum; // long-range summation algorithm
		std::mt19937_64 mersenne_twister = std::mt19937_64(0);
		vec3<T> side; // box sides
		T M = 0, Z = 0, Z2 = 0, tracedisp6 = 0, sumdisp6 = 0, sumdisp62 = 0, tracedisp12 = 0, sumdisp12 = 0;
			// total mass, total charge, sum of charges^2, trace of dispersion matrix, sum (of squares) of all dispersion terms coefficients
		unsigned n = 0, dof = 0; // total number of atoms, degrees of freedom
		T temperature_ref, pressure_ref; // reference temperature (for constant-T and initialization), reference pressure (for constant-P)
		T energy_lrc, energy_intra_coulomb, energy_intra_lj; // LJ long-range correction, intra-molecular energy corrections
		T t = 0, potential, virial, D; // time, potential energy, virial, diffusion coefficient
		bool rescale_temperature = false, first_step = true, dispersion_correction = true;

		molecular_system(T temp = 298.15, T p = 1*atm<T>, vec3<T> side = 50, T D = DW25<T>, LRSum lrsum = LRSum())
		// constructor:
		//	`temp` is the temperature to give to the initial configuration and reference temperature (in Kelvin)
		//	`p` is the reference pressure (in AKMA units)
		//	`side` is the array of sides that characterize the dimensions of the box of the system to simulate (in angstrom)
		//	`D` is the diffusion coefficient used in stochastic integrators (in AKMA units)
		//	`lrsum` is the long-range summation algorithm to be used
			: lrsum(lrsum), side(side), temperature_ref(temp), pressure_ref(p), D(D)
		{}

		void add_molecule(const molecule<T>& mol, const vec3<T>& pos = 0)
		// add a molecule to the system.
		// `mol` is the molecule to add.
		// `pos` is the position of the molecule.
		// Throw a `std::runtime_error` exception if the number of atoms is going
		// to exceed the maxium allowed (i.e. `max_atoms`).
		{
			using std::size_t;
			using std::sqrt;

			if (mol.n == 0)
				return;

			if (!first_step)
				copy_to_vector();

			if (n + mol.n > max_atoms)
				throw std::runtime_error("Error: Exceeded maximum number of atoms");

			x_tmp.resize(n + mol.n); p_tmp.resize(n + mol.n);
			v_tmp.resize(n + mol.n); f_tmp.resize(n + mol.n);
			m_tmp.resize(n + mol.n);

			for (size_t i = 0; i < mol.n; ++i)
				x_tmp[n + i] = mol.x[i] + pos;
			for (auto atom_id : mol.id)
				id.push_back(atom_id);
			for (auto q : mol.q)
				z.push_back(q*sqrtkC<long double>);

			for (auto atom_id : mol.id)
				lj_sqrteps.push_back(sqrt(lj_params<long double>[int(atom_id)].epsilon));
			for (auto atom_id : mol.id)
				lj_halfR.push_back(lj_params<T>[int(atom_id)].half_Rmin);
			for (size_t i = 0; i < mol.n; ++i)
				m_tmp[n + i] = atom_mass<T>[ int(id[n + i]) ];
			for (size_t i = 0; i < mol.n; ++i)
				M += m_tmp[n + i];
			// momenta are initialized by sampling from a gaussian distribution
			// with variance proportional to mass and temperature. With a proper
			// configuration, this would correspond to a Maxwell-Boltzmann distribution.
			for (size_t i = 0; i < mol.n; ++i)
				p_tmp[n + i] = gen_gaussian() * sqrt(m_tmp[n + i] * kT_ref());

			for (auto atom_bonds : mol.bonds)
			{
				for (size_t j = 0; j < atom_bonds.n; ++j)
					atom_bonds[j] += n;
				bonds.push_back(atom_bonds);
			}

			for (size_t i = 0; i < mol.impropers.size(); ++i)
				impropers.push_back({mol.impropers[i][0] + n,
									 mol.impropers[i][1] + n,
									 mol.impropers[i][2] + n,
									 mol.impropers[i][3] + n});

			// update the sum and the sum square of the charges
			for (size_t i = 0; i < mol.n; ++i)
			{
				Z += z[n + i];
				Z2 += z[n + i]*z[n + i];
			}
			// update the trace and the sum of the dispersion coefficients matrix.
			// Dispersion coefficients are calculated using the Lorentz-Berthelot mixing rule.
			// These coefficients are used for long-range dispersion energy correction. If not
			// needed, set `dispersion_correction` to false for faster initialization.
			if (dispersion_correction)
			{
				for (size_t i = 0; i < mol.n; ++i)
				{
					T eps = lj_sqrteps[n + i]*lj_sqrteps[n + i];
					T C6 = 2*lj_halfR[n + i];
					C6 *= C6;
					C6 = C6 * C6 * C6;
					T C12 = eps * C6 * C6;
					C6 = 2 * eps * C6;
					tracedisp6 += C6;
					sumdisp6 += C6;
					sumdisp62 += C6 * C6;
					tracedisp12 += C12;
					sumdisp12 += C12;
				}
				for (size_t i = n; i < n+mol.n; ++i)
					for (size_t j = 0; j < i; ++j)
					{
						T epsij = lj_sqrteps[i]*lj_sqrteps[j];
						T C6 = lj_halfR[i]+lj_halfR[j];
						C6 *= C6;
						C6 = C6 * C6 * C6;
						T C12 = epsij * C6 * C6;
						C6 = 2 * epsij * C6;
						sumdisp6 += 2*C6;
						sumdisp62 += 2*C6*C6;
						sumdisp12 += 2*C12;
					}
			}
			else
				tracedisp6 = sumdisp6 = sumdisp62 = tracedisp12 = sumdisp12 = 0;

			n += mol.n;
			dof = 3*n;
			first_step = true;
		}

		void primitive_cubic_lattice(
			const vec3i&       n_side,
			T                  lattice_const,
			const molecule<T>& species1,
			const molecule<T>& species2 = empty_molecule<T>
		)
		// create a primitive cubic lattice (overwriting the current system). If two chemical
		// species are given, create an interpenetrating (caesium-chloride-like) structure.
		// Examples: Caesium chloride (CsCl), caesium bromide (CsBr), caesium iodide (CsI), many
		// binary metallic alloys.
		// `n_side` is the number of cubic unit cells along each dimension. The total number of lattice
		// points for species 1 (or 2) is thus n_side[0] * n_side[1] * n_side[2].
		// `lattice_const` is the distance between nearest cubic unit cells (i.e. the linear size
		// of a cubic unit cell).
		// `species1` is the first chemical species to insert in the lattice.
		// `species2` is the second chemical species to insert in the lattice. By default it
		// is a empty molecule (no atoms are placed).
		// If species1 == species2 the result is a body-centered cubic lattice for a single
		// species.
		// After this method is called, calling `fetch` is not needed.
		{
			clear();
			side = vec3<T>(n_side)*lattice_const;
			vec3<T> hside = vec3<T>(n_side)/2 - T(0.25);
			for (int i = 0; i < n_side[0]; ++i)
				for (int j = 0; j < n_side[1]; ++j)
					for (int k = 0; k < n_side[2]; ++k)
						add_molecule(species1, {(i-hside[0])*lattice_const, (j-hside[1])*lattice_const, (k-hside[2])*lattice_const});
			hside -= T(0.5);
			for (int i = 0; i < n_side[0]; ++i)
				for (int j = 0; j < n_side[1]; ++j)
					for (int k = 0; k < n_side[2]; ++k)
						add_molecule(species2, {(i-hside[0])*lattice_const, (j-hside[1])*lattice_const, (k-hside[2])*lattice_const});
			fetch();
		}

		void face_centered_cubic_lattice(
			const vec3i&       n_side,
			T                  lattice_const,
			const molecule<T>& species1,
			const molecule<T>& species2 = empty_molecule<T>
		)
		// create a face-centered cubic lattice (overwriting the current system). If two chemical
		// species are given, create a rock-salt-like structure.
		// Examples: Rock salt (NaCl), lithium chloride (LiCl), lithium floride (LiF), many
		// salts composed with other alkali metals (expect caesium, and probably francium).
		// `n_side` is the number of cubic unit cells along each dimension. The total number of lattice
		// points for species 1 (or 2) is thus 4 * n_side[0] * n_side[1] * n_side[2].
		// `lattice_const` is the distance between nearest cubic unit cells (i.e. the linear size
		// of a cubic unit cell).
		// `species1` is the first chemical species to insert in the lattice.
		// `species2` is the second chemical species to insert in the lattice. By default it
		// is a empty molecule (no atoms are placed).
		// If species1 == species2 the result is a primitive cubic lattice for a single species
		// with 8x the number of cubic unit cells (1/2x the linear cell size).
		// After this method is called, calling `fetch` is not needed.
		{
			clear();
			vec3i dside = n_side*2;
			side = vec3<T>(n_side)*lattice_const;
			T dist = lattice_const/2;
			for (int i = 0; i < dside[0]; ++i)
				for (int j = 0; j < dside[1]; ++j)
					for (int k = 0; k < dside[2]; ++k)
					{
						physics::vec3d p = {(i-n_side[0] + T(0.25))*dist, (j-n_side[1] + T(0.25))*dist, (k-n_side[2] + T(0.25))*dist};
						if (math::is_even(i+j+k))
							add_molecule(species1, p);
						else
							add_molecule(species2, p);
					}
			fetch();
		}

		T kinetic_energy()
		// kinetic energy of the system in kcal/mol
		{
			fetch();

			kin = 0;
			for (std::size_t i = 0; i < n; ++i)
				kin += dot(p[i], p[i]) / m[i];
			return kin/2;
		}

		T internal_energy()
		// total energy of the system in kcal/mol.
		// `force` must be called before this method so that the potential is updated.
		{
			return kinetic_energy() + potential;
		}

		T temperature()
		// instantaneous temperature of the system, in K
		{
			return 2*kinetic_energy() / (dof*kB<T>);
		}

		T kT()
		// instantaneous temperature multiplied by Boltzmann constant, in kcal/mol
		{
			return kB<T> * temperature();
		}

		T kT_ref() const noexcept
		// reference temperature multiplied by Boltzmann constant, in kcal/mol
		{
			return kB<T> * temperature_ref;
		}

		T volume() const noexcept
		// return the volume of the simulation box
		{
			return side[0]*side[1]*side[2]; // in angstrom^3 = 10^-30 m^3
		}

		T density() const noexcept
		// return the average mass density in the simulation box
		{
			return M / volume(); // in amu/angstrom^3
		}

		T number_density() const noexcept
		// return the average number density in the simulation box
		{
			return n / volume(); // in angstrom^-3
		}

		T pressure()
		// instantaneous pressure calculated with instantaneous temperature in kcal/(mol angstrom^3)
		{
			return (2*kinetic_energy() + virial)/(3*volume());
		}

		T pressure_T_ref() const noexcept
		// instantaneous pressure calculated with reference temperature in kcal/(mol angstrom^3)
		{
			return (n*kT_ref() + virial/3) / volume();
		}

		T external_pressure()
		// instantaneous pressure calculated from external virial in kcal/(mol angstrom^3);
		// not valid if the system is periodic
		{
			fetch();
			return (virial - dot(x, f)) / (3*volume());
		}

		const state<T, 3>& force(bool eval = true) override
		// calculate all forces and update energies;
		// if `eval` is false, reuse old values
		{
			fetch();
			if (eval)
			{
				f = 0;
				potential = 0;

				force_dihedrals();
				force_impropers();

				force_angles();
				force_bonds();

				virial = dot(x, f);

				force_nonbonded();

				if constexpr (std::is_same_v<LRSum, direct<T, state<T, 3>>>)
					diff_box_confining(2);
			}
			return f;
		}

		const state<T, 3>& force_long(bool eval = true)
		// calculate long-range forces only;
		// if `eval` is false, reuse old values
		// required by multi_timestep_leapfrog integrator
		{
			fetch();
			if (eval)
			{
				f = 0;

				force_nonbonded();
			}
			return f;
		}

		const state<T, 3>& force_short(bool eval = true)
		// calculate short-range forces only;
		// if `eval` is false, reuse old values
		// required by multi_timestep_leapfrog integrator
		{
			fetch();
			if (eval)
			{
				f = 0;
				potential = 0;

				force_dihedrals();
				force_impropers();

				force_angles();
				force_bonds();

				virial = dot(x, f);

				if constexpr (std::is_same_v<LRSum, direct<T, state<T, 3>>>)
					diff_box_confining(2);
			}
			return f;
		}

		const state<T, 3>& vel(bool eval = true) override
		// calculate velocities from momenta;
		// if `eval` is false, reuse old values
		{
			fetch();
			if (eval)
				for (std::size_t i = 0; i < n; ++i)
					v[i] = p[i] / m[i];
			return v;
		}

		void rand()
		// generate random numbers for use in stochastic integrators (for Langevin's equations)
		// required by `stochastic_leapfrog`
		{
			using std::sqrt;
			fetch();

			std::ranges::generate(noise, gen_gaussian);
			for (size_t i = 0; i < n; ++i)
				noise[i] *= sqrt(m[i] * kT_ref());
			gamma = (kT_ref() / D) / m;
		}

		void rescale_temp()
		// rescale the temperature of the system to a temperature fixed by the member variable: temperature_ref
		{
			using std::sqrt;

			fetch();

			T temp = temperature();
			if (temp > 0)
				p *= sqrt(temperature_ref / temp);
		}

		void fetch()
		// during simulation, a std::valarray is used because its implementation supports efficient
		// expression templates through operator overloading, making it a more useful choice
		// with respect to std::vector. During initialization, std::vector is preferred because it allows
		// an efficient implementation of `resize`.
		// Important: this function needs to be explicitly called if x,p,m need to be accessed
		// right after `add_molecule` is called (i.e. before any of `step`, `force`, `force_short`
		// `force_long`, `vel`, `rand`, `kinetic_energy`, `total_energy`, `temperature`, `kT`,
		// `pressure`, `external_pressure`, `primitive_cubic_lattice`, `face_centered_cubic_lattice`
		// is called).
		{
			// if valarray have already the correct size, just return
			if (x.size() == n)
				return;
			x.resize(n); p.resize(n);
			v.resize(n); f.resize(n);
			m.resize(n);
			noise.resize(n); gamma.resize(n);
			copy(begin(x_tmp), end(x_tmp), begin(x));
			copy(begin(p_tmp), end(p_tmp), begin(p));
			copy(begin(v_tmp), end(v_tmp), begin(v));
			copy(begin(f_tmp), end(f_tmp), begin(f));
			copy(begin(m_tmp), end(m_tmp), begin(m));

			zero_total_momentum();
			rescale_temp();
		}

		void zero_total_momentum()
		// set total momentum to zero
		{
			fetch();
			p -= p.sum()/n;
		}

		void random_momenta()
		// randomize atoms momenta based on Maxwell-Boltzmann distribution
		{
			fetch();

			for (size_t i = 0; i < n; ++i)
				p[i] = gen_gaussian() * sqrt(m[i] * kT_ref());

			zero_total_momentum();
			rescale_temp();
		}

		private:

		void copy_to_vector()
		// when many molecules are added, a std::vector is used, since std::valarray<T>::resize
		// reallocates and zeroes all memory at each call, which would lead to poor performance
		// during initialization
		{
			x_tmp.resize(n); p_tmp.resize(n);
			v_tmp.resize(n); f_tmp.resize(n);
			m_tmp.resize(n);
			copy(begin(x), end(x), begin(x_tmp));
			copy(begin(p), end(p), begin(p_tmp));
			copy(begin(v), end(v), begin(v_tmp));
			copy(begin(f), end(f), begin(f_tmp));
			copy(begin(m), end(m), begin(m_tmp));
		}

		vec3<T> gen_gaussian()
		// generate a vec3 with three random numbers sampled from a normal distribution
		// with 0 mean and 1 std.dev.
		{
			return vec3<T>(
				n_dist(mersenne_twister),
				n_dist(mersenne_twister),
				n_dist(mersenne_twister)
			);
		}

		void force_bonds()
		// calculate 1-2 bonded forces and energies for all atoms
		{
			using std::size_t;

			for (size_t i = 0; i < n; ++i)
				for (size_t k = 0; k < bonds[i].n; ++k)
				{
					const size_t j = bonds[i][k];
					if (!(i < j))
						continue;
					atom_type idi = id[i], idj = id[j];
					if (idi > idj)
						std::swap(idi, idj);
					const auto& par = bond_params<T>[{ idi, idj }];

					vec3<T> r = x[i] - x[j];
					T d = norm(r);
					T d_ = 1/d;
					T diff = d - par.b0;

					r = (2 * par.Kb * diff * d_) * r;
					f[i] -= r;
					f[j] += r;
					potential += par.Kb * diff*diff;
				}
		}

		void force_angles()
		// calculate 1-3 bonded forces and energies for all atoms
		{
			using std::sqrt;
			using std::atan2;
			using std::size_t;

			for (size_t i = 0; i < n; ++i)
				for (size_t mj = 0; mj < bonds[i].n; ++mj)
				{
					size_t j = bonds[i][mj];
					atom_type idj = id[j];

					vec3<T> rij = x[i] - x[j];
					T rij2 = dot(rij, rij);

					for (size_t mk = 0; mk < bonds[j].n; ++mk)
					{
						size_t k = bonds[j][mk];
						if (!(i < k))
							continue;
						atom_type idi = id[i], idk = id[k];
						if (idi > idk)
							std::swap(idi, idk);
						const auto& par = angle_params<T>[{ idi, idj, idk }];

						vec3<T> rkj = x[k] - x[j];
						T rkj2 = dot(rkj, rkj);

						T num = dot(rij, rkj), den = sqrt(rij2 * rkj2);
						T costheta = num / den;
						T sintheta = sqrt((1 + costheta) * (1 - costheta));
						T theta = atan2(sintheta, costheta);
						T diff = theta - par.theta0;
						T Kthetasinden = 2 * par.Ktheta * diff / (sintheta * den);
						vec3<T> rijt = rij * (num / rij2), rkjt = rkj * (num / rkj2);
						// Urey-Bradley potential
						vec3<T> rik = x[i] - x[k];
						T dik = norm(rik);
						T dik_ = 1/dik;
						T diffUB = dik - par.S0;

						rik = (2 * par.Kub * diffUB * dik_) * rik;

						f[i] += (rkj - rijt) * Kthetasinden - rik;
						f[j] += (rijt + rkjt - rij - rkj) * Kthetasinden;
						f[k] += (rij - rkjt) * Kthetasinden + rik;
						potential += par.Ktheta * diff*diff + par.Kub * diffUB*diffUB;
					}
				}
		}

		void force_dihedrals()
		// (TODO) calculate dihedral angles forces and energies for all atoms
		{}

		void force_impropers()
		// (TODO) calculate improper angles forces and energies for all atoms
		{}

		void eval_nonbonded(std::size_t i, std::size_t j, T cutoff2) noexcept
		// evaluate non-bonded forces and energies for a pair of atoms with indices `i`, `j`.
		// `cutoff2` is the square of the Lennard-Jones cutoff radius (0 for no cutoff).
		// This method is called inside `correct_nonbonded`.
		{
			vec3<T> r = x[i] - x[j];
			T r2 = dot(r, r), r2_ = 1/r2;

			T d_ = sqrt(r2_);
			T epsijs = 0, s = 0;
			if (cutoff2 == 0 || r2 <= cutoff2)
			{
				T Rij = lj_halfR[i] + lj_halfR[j];
				s = Rij * Rij * r2_;
				s = s * s * s;
				epsijs = lj_sqrteps[i] * lj_sqrteps[j] * s;
			}
			T coulomb = z[i] * z[j] * d_;
			T vir = -(12 * epsijs * (s - 1) + coulomb);
			vec3<T> fij = (vir * r2_) * r;
			f[i] += fij;
			f[j] -= fij;
			energy_intra_coulomb -= coulomb;
			energy_intra_lj -= epsijs * (s - 2);
			virial += vir;
		}

		void correct_nonbonded(T cutoff2)
		// the long-range forces must be excluded for 1-2 and 1-3 bonded atoms.
		// This correction is applied so that the force calculated through a long-range
		// summation algorithm vanish for these atom pairs.
		// `cutoff2` is the square of the Lennard-Jones cutoff radius (0 for no cutoff).
		{
			using std::size_t;
			energy_intra_coulomb = 0;
			energy_intra_lj = 0;
			for (size_t i = 0; i < n; ++i)
				for (size_t mj = 0; mj < bonds[i].n; ++mj)
				{
					const size_t j = bonds[i][mj];
					if (i < j)
						eval_nonbonded(i, j, cutoff2);
					for (size_t mk = 0; mk < bonds[j].n; ++mk)
					{
						const size_t k = bonds[j][mk];
						if (i < k)
							eval_nonbonded(i, k, cutoff2);
					}
				}
			potential += energy_intra_coulomb + energy_intra_lj;
		}

		void force_nonbonded()
		// calculate non-bonded forces and energies (both electrostatic and Lennard-Jones)
		// for all atoms.
		{
			// call long-range summation algorithm operator()
			lrsum(*this, tp);

			// subtract non-bonded forces for bonded atoms
			T cutoff = lrsum.cutoff_radius();
			correct_nonbonded(cutoff*cutoff);

			// the following applies a long-range correction to the energy due to dispersion interactions
			// since they are truncated at a cutoff radius (see Allen et al., 2017, pp. 79-81)
			if (cutoff)
			{
				T cutoff3 = cutoff*cutoff*cutoff;
				T cutoff9 = cutoff3*cutoff3*cutoff3;

				T common_factor = 2 * std::numbers::pi_v<T> / volume();
				T term6 = common_factor * sumdisp6 / (3*cutoff3);
				T term12 = common_factor * sumdisp12 / (9*cutoff9);
				energy_lrc = term12 - term6;
				potential += energy_lrc;
				virial += 12*term12 - 6*term6;
			}
		}

		void diff_box_confining(T steepness)
		// a potential well or box may be approximated by a potential which gives a constant
		// force outside a cube and no force inside it
		// `steepness` is the steepness of the potential outside the cube, which is equal
		// to the force applied
		{
			vec3<T> hside = side / 2;
			for (std::size_t i = 0; i < n; ++i)
			{
				for (std::size_t j = 0; j < 3; ++j)
				{
					if (x[i][j] > hside[j])
					{
						f[i][j] -= steepness;
						potential += (x[i][j] - hside[j])*steepness;
					}

					if (x[i][j] < -hside[j])
					{
						f[i][j] += steepness;
						potential += (-hside[j] - x[i][j])*steepness;
					}
				}
			}
		}

		public:

		// constructors and utilities that occupy lots of space and thus they are put at the end

		void clear()
		// clear all data so that the object may be reused for a new simulation
		{
			x.resize(0); p.resize(0);
			v.resize(0); f.resize(0);
			m.resize(0);
			noise.resize(0); gamma.resize(0);
			x_tmp.clear(); p_tmp.clear();
			v_tmp.clear(); f_tmp.clear();
			m_tmp.clear();
			z.clear();
			lj_sqrteps.clear();
			lj_halfR.clear();
			id.clear();
			bonds.clear();
			impropers.clear();
			M = Z = Z2 = 0;
			tracedisp6 = sumdisp6 = 0;
			tracedisp12 = sumdisp12 = 0;
			n = dof = 0;
			t = 0;
			first_step = true;
		}

		molecular_system(const molecular_system& other)
		// copy constructor
			: x(other.x), p(other.p), v(other.v), f(other.f), noise(other.noise), m(other.m), gamma(other.gamma),
			x_tmp(other.x_tmp), p_tmp(other.p_tmp), v_tmp(other.v_tmp), f_tmp(other.f_tmp), m_tmp(other.m_tmp),
			z(other.z), lj_sqrteps(other.lj_sqrteps), lj_halfR(other.lj_halfR), id(other.id), bonds(other.bonds), impropers(other.impropers),
			mersenne_twister(other.mersenne_twister),
			side(other.side), M(other.M), Z(other.Z), Z2(other.Z2), tracedisp6(other.tracedisp6), sumdisp6(other.sumdisp6),
			sumdisp62(other.sumdisp62), tracedisp12(other.tracedisp12), sumdisp12(other.sumdisp12), n(other.n), dof(other.dof),
			temperature_ref(other.temperature_ref),
			energy_lrc(other.energy_lrc), energy_intra_coulomb(other.energy_intra_coulomb), energy_intra_lj(other.energy_intra_lj),
			t(other.t), potential(other.potential), virial(other.virial), D(other.D),
			first_step(other.first_step), dispersion_correction(other.dispersion_correction)
		{}

		template <template <typename, typename> typename LRSummationU>
		molecular_system(const molecular_system<T, LRSummationU>& other)
		// template copy constructor
			: x(other.x), p(other.p), v(other.v), f(other.f), noise(other.noise), m(other.m), gamma(other.gamma),
			x_tmp(other.x_tmp), p_tmp(other.p_tmp), v_tmp(other.v_tmp), f_tmp(other.f_tmp), m_tmp(other.m_tmp),
			z(other.z), lj_sqrteps(other.lj_sqrteps), lj_halfR(other.lj_halfR), id(other.id), bonds(other.bonds), impropers(other.impropers),
			mersenne_twister(other.mersenne_twister),
			side(other.side), M(other.M), Z(other.Z), Z2(other.Z2), tracedisp6(other.tracedisp6), sumdisp6(other.sumdisp6),
			sumdisp62(other.sumdisp62), tracedisp12(other.tracedisp12), sumdisp12(other.sumdisp12), n(other.n), dof(other.dof),
			temperature_ref(other.temperature_ref),
			energy_lrc(other.energy_lrc), energy_intra_coulomb(other.energy_intra_coulomb), energy_intra_lj(other.energy_intra_lj),
			t(other.t), potential(other.potential), virial(other.virial), D(other.D),
			first_step(other.first_step), dispersion_correction(other.dispersion_correction)
		{}

		molecular_system& operator=(const molecular_system& other)
		// copy-assignment operator
		{
			x = other.x; p = other.p;
			v = other.v; f = other.f;
			noise = other.noise;
			m = other.m;
			gamma = other.gamma;
			x_tmp = other.x_tmp; p_tmp = other.p_tmp;
			v_tmp = other.v_tmp; f_tmp = other.f_tmp;
			m_tmp = other.m_tmp;
			z = other.z;
			lj_sqrteps = other.lj_sqrteps;
			lj_halfR = other.lj_halfR;
			id = other.id;
			bonds = other.bonds;
			impropers = other.impropers;
			mersenne_twister = other.mersenne_twister;
			side = other.side;
			M = other.M;
			Z = other.Z;
			Z2 = other.Z2;
			tracedisp6 = other.tracedisp6;
			sumdisp6 = other.sumdisp6;
			sumdisp62 = other.sumdisp62;
			tracedisp12 = other.tracedisp12;
			sumdisp12 = other.sumdisp12;
			n = other.n;
			dof = other.dof;
			temperature_ref = other.temperature_ref;
			energy_lrc = other.energy_lrc;
			energy_intra_coulomb = other.energy_intra_coulomb;
			energy_intra_lj = other.energy_intra_lj;
			t = other.t;
			potential = other.potential;
			virial = other.virial;
			D = other.D;
			first_step = other.first_step;
			dispersion_correction = other.dispersion_correction;
			return *this;
		}

		template <template <typename, typename> typename LRSummationU>
		molecular_system& operator=(const molecular_system<T, LRSummationU>& other)
		// template copy-assignment operator
		{
			x = other.x; p = other.p;
			v = other.v; f = other.f;
			noise = other.noise;
			m = other.m;
			gamma = other.gamma;
			x_tmp = other.x_tmp; p_tmp = other.p_tmp;
			v_tmp = other.v_tmp; f_tmp = other.f_tmp;
			m_tmp = other.m_tmp;
			z = other.z;
			lj_sqrteps = other.lj_sqrteps;
			lj_halfR = other.lj_halfR;
			id = other.id;
			bonds = other.bonds;
			impropers = other.impropers;
			mersenne_twister = other.mersenne_twister;
			side = other.side;
			M = other.M;
			Z = other.Z;
			Z2 = other.Z2;
			tracedisp6 = other.tracedisp6;
			sumdisp6 = other.sumdisp6;
			sumdisp62 = other.sumdisp62;
			tracedisp12 = other.tracedisp12;
			sumdisp12 = other.sumdisp12;
			n = other.n;
			dof = other.dof;
			temperature_ref = other.temperature_ref;
			energy_lrc = other.energy_lrc;
			energy_intra_coulomb = other.energy_intra_coulomb;
			energy_intra_lj = other.energy_intra_lj;
			t = other.t;
			potential = other.potential;
			virial = other.virial;
			D = other.D;
			first_step = other.first_step;
			dispersion_correction = other.dispersion_correction;
			return *this;
		}
	};

} // namespace physics

#endif // PHYSICS_MOLECULE_H
































