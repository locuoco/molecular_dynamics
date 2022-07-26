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
#include <algorithm> // min, max, ranges::generate
#include <cmath> // sqrt, atan2
#include <random> // mt19937_64, normal_distribution
#include <stdexcept> // runtime_error
#include <numbers> // numbers::pi_v

#include "../math/helper.hpp" // deg2rad, two1_6

#include "tensor.hpp"
#include "integrator/integrator.hpp"
#include "pppm.hpp"
#include "../utils/thread_pool.hpp"

namespace physics
{
	// Coulomb constant in AKMA units
	template <std::floating_point T>
	constexpr T kC = 332.06371330062682040436798930743L;

	// square root of the Coulomb constant in AKMA units
	template <std::floating_point T>
	constexpr T sqrtkC = 18.222615435239443706926657644703L;

	// Boltzmann constant in AKMA units (kcal/(mol K))
	template <std::floating_point T>
	constexpr T kB = 0.00198720425864083173996175908222L;

	// Avogadro number (mol^-1)
	template <std::floating_point T>
	constexpr T NA = 6.02214076e23L;

	// 1 atmosphere in AKMA units
	template <std::floating_point T>
	constexpr T atm = 1.4583972574259082217973231357553e-5L;

	// 1 kg/m^3 in AKMA units
	template <std::floating_point T>
	constexpr T kg_per_m3 = 6.02214076e-4L;

	// water self-diffusion constant at 5 °C in AKMA units
	template <std::floating_point T = double>
	constexpr T DW5 = 0.00616480364668354400695529510607L;

	// water self-diffusion constant at 25 °C in AKMA units
	template <std::floating_point T = double>
	constexpr T DW25 = 0.01123940014569822971609058164065L;

	enum class atom_type : unsigned char { HT, OT, HTL, OTL, HF, OF, N };
	constexpr std::size_t num_atom_types = std::size_t(atom_type::N);

	// maximum bond coordination number assumed (not more than `max_bonds` bonds per atoms are allowed)
	constexpr std::size_t max_bonds = 4;

	constexpr unsigned char atom_number[num_atom_types] = {1, 8, 1, 8, 1, 8};

	// atom masses in atomic units
	template <std::floating_point T>
	constexpr T atom_mass[num_atom_types] = {1.008L, 15.9994L, 1.008L, 15.9994L, 1.008L, 15.9994L};

	template <std::floating_point T>
	struct lj_struct { T epsilon, half_Rmin; };

	template <std::floating_point T>
	struct bond_struct { T Kr, r0; };

	template <std::floating_point T>
	struct angle_struct { T Kangle, angle0, KUB, rUB; };

	template <std::floating_point T>
	lj_struct<T> lj_params[num_atom_types] =
	// Lennard-Jones parameters
	// (epsilon, R_min/2) associated to each atom types
	{
		{0.046L, 0.4L*math::two1_6<long double>/2}, // HT
		{0.1521L, 3.1507L*math::two1_6<long double>/2}, // OT
		{0, 0}, // HTL
		{0.102L, 3.188L*math::two1_6<long double>/2}, // OTL
		{0, 0}, // HF
		{0.18936998087954110898661567877629L, 3.1776L*math::two1_6<long double>/2}, // OF
	};

	template <std::floating_point T>
	std::map<std::pair<atom_type, atom_type>, bond_struct<T>> bond_params =
	// 1-2 bonds coefficients
	// (Kr, r0) parameters associated to a pair of adjacent atoms
	{
		{ {atom_type::HT, atom_type::OT}, {450, 0.9572L} },
		{ {atom_type::HTL, atom_type::OTL}, {450, 0.9572L} },
		{ {atom_type::HF, atom_type::OF}, {358.50860420650095602294455066922L, 1.027L} }
	};

	template <std::floating_point T>
	std::map<std::tuple<atom_type, atom_type, atom_type>, angle_struct<T>> angle_params =
	// 1-3 bonds coefficients
	// (Kangle, angle0, KUB, rUB) parameters associated to a triplet of adjacent atoms
	{
		{ {atom_type::HT, atom_type::OT, atom_type::HT}, {55, math::deg2rad(104.52L), 0, 0} },
		{ {atom_type::HTL, atom_type::OTL, atom_type::HTL}, {55, math::deg2rad(104.52L), 0, 0} },
		{ {atom_type::HF, atom_type::OF, atom_type::HF}, {45.769598470363288718929254302103L, math::deg2rad(114.70L), 0, 0} }
	};

	template <std::floating_point T = double>
	struct molecule
	{
		state<T, 3> x; // position
		std::vector<atom_type> id; // identity of the atom
		std::vector<T> part_q; // partial charges
		std::vector<fixed_list<max_bonds>> bonds; // bonds list as an adjacency list for each atom
		std::vector<std::array<unsigned int, 4>> impropers; // list of impropers in the whole molecule
		unsigned int n; // number of atoms in the molecule
	};

	template <std::floating_point T = double>
	const molecule<T> water_tip3p
	// TIP3P water model
	{
		.x = { {-0.75695032726366116174L, 0.58588227661829493656L}, {0}, {0.75695032726366116174L, 0.58588227661829493656L} },
			// oxygen at the origin
		.id = {atom_type::HT, atom_type::OT, atom_type::HT},
		.part_q = {0.417L, -0.834L, 0.417L},
		.bonds = { {1}, {0, 2}, {1} },
		.impropers = {}, // no impropers
		.n = 3
	};

	template <std::floating_point T = double>
	const molecule<T> water_tip3p_lr
	// TIP3P water model optimized for long-range interactions
	// see D. J. Price, C. L. Brooks, "A modified TIP3P water potential for simulation with Ewald summation"
	{
		.x = { {-0.75695032726366116174L, 0.58588227661829493656L}, {0}, {0.75695032726366116174L, 0.58588227661829493656L} },
			// oxygen at the origin
		.id = {atom_type::HTL, atom_type::OTL, atom_type::HTL},
		.part_q = {0.415L, -0.830L, 0.415L},
		.bonds = { {1}, {0, 2}, {1} },
		.impropers = {}, // no impropers
		.n = 3
	};

	template <std::floating_point T = double>
	const molecule<T> water_fba_eps
	// FBA/eps water model
	// see R. Fuentes-Azcatl, M. Barbosa, "Flexible Bond and Angle, FBA/epsilon model of water", 2018
	{
		.x = { {-0.75695032726366116174L, 0.58588227661829493656L}, {0}, {0.75695032726366116174L, 0.58588227661829493656L} },
			// oxygen at the origin
		.id = {atom_type::HF, atom_type::OF, atom_type::HF},
		.part_q = {0.4225L, -0.8450L, 0.4225L},
		.bonds = { {1}, {0, 2}, {1} },
		.impropers = {}, // no impropers
		.n = 3
	};

	template <
		std::floating_point T = double,
		template <typename, typename> typename LRSummationT = pppm,
		template <typename, typename, template <typename...> typename...> typename IntegT = leapfrog,
		template <typename...> typename ... IntegPars
	>
	requires integrator<IntegT<T, state<T, 3>, IntegPars...>, T, state<T, 3>>
	class molecular_system : public physical_system_base<T, state<T, 3>>
	{
		using LRSum = LRSummationT<T, state<T, 3>>;
		using Integ = IntegT<T, state<T, 3>, IntegPars...>;

		std::normal_distribution<T> n_dist = std::normal_distribution<T>(0, 1);
		mutable T kin; // kinetic energy
		mutable bool kin_updated = false;

		public:

		using scalar_type = T;

		state<T, 3> x, p, v, f; // position, momentum, velocity, force
		state<T, 3> noise; // gaussian noise
		std::valarray<T> m, gamma; // mass, damping factor
		std::vector<vec3<T>> x_tmp, p_tmp, v_tmp, f_tmp; // used for initialization
		std::vector<T> m_tmp; // used for initialization
		std::vector<T> z, lj_sqrteps, lj_halfR; // normalized partial charges + Lennard-Jones parameters
		std::vector<atom_type> id; // identity of the atom
		std::vector<fixed_list<max_bonds>> bonds; // bonds as an adjacency list
		std::vector<std::array<unsigned int, 4>> impropers; // list of impropers in the whole system
		T M = 0, Z = 0, Z2 = 0, tracedisp = 0, sumdisp = 0;
			// total mass, total charge, sum of charges^2, trace of dispersion matrix, sum of all dispersion terms coefficients
		unsigned int n = 0, dof = 0; // total number of atoms, degrees of freedom
		utils::thread_pool tp; // thread pool
		Integ integ; // integrator
		LRSum lrsum; // long-range summation algorithm
		std::mt19937_64 mersenne_twister = std::mt19937_64(0);
		T side, temperature_ref;
		T t = 0, potential, virial, D; // time, potential energy, virial, diffusion coefficient
		bool rescale_temperature = false, first_step = true;

		static const std::size_t max_atoms = 2'000'000;

		molecular_system(T side = 50, T temp = 298.15, T D = DW25<T>, Integ integ = Integ(), LRSum lrsum = LRSum())
		// constructor:
		// `side` is the side of the cubic box of the system to simulate (in angstrom)
		// `temp` is the temperature to give to the initial configuration and reference temperature (in Kelvin)
		// `D` is the diffusion coefficient used in stochastic integrators (in AKMA units)
		// `integ` is the integrator to be used
		// `lrsum` is the long-range summation algorithm to be used
			: integ(integ), lrsum(lrsum), side(side), temperature_ref(temp), D(D)
		{}

		void add_molecule(const molecule<T>& mol, const vec3<T>& pos = 0)
		{
			using std::size_t;
			using std::sqrt;

			if (!first_step)
				copy_to_vector();

			if (n + mol.n > max_atoms)
				throw std::runtime_error("Error: Exceeded maximum number of atoms");
			x_tmp.resize(n + mol.n); p_tmp.resize(n + mol.n);
			v_tmp.resize(n + mol.n); f_tmp.resize(n + mol.n);
			m_tmp.resize(n + mol.n);

			for (size_t i = 0; i < mol.n; ++i)
				x_tmp[n + i] = mol.x[i] + pos;
			for (size_t i = 0; i < mol.n; ++i)
				id.push_back(mol.id[i]);
			for (size_t i = 0; i < mol.n; ++i)
				z.push_back(mol.part_q[i]*sqrtkC<long double>);

			for (size_t i = 0; i < mol.n; ++i)
				lj_sqrteps.push_back(sqrt(lj_params<long double>[int(mol.id[i])].epsilon));
			for (size_t i = 0; i < mol.n; ++i)
				lj_halfR.push_back(lj_params<T>[int(mol.id[i])].half_Rmin);
			for (size_t i = 0; i < mol.n; ++i)
				m_tmp[n + i] = atom_mass<T>[ int(id[n + i]) ];
			for (size_t i = 0; i < mol.n; ++i)
				M += m_tmp[n + i];
			// momenta are initialized by sampling from a gaussian distribution
			// with variance proportional to mass and temperature
			for (size_t i = 0; i < mol.n; ++i)
				p_tmp[n + i] = gen_gaussian() * sqrt(m_tmp[n + i] * kT_ref());

			for (size_t i = 0; i < mol.n; ++i)
			{
				fixed_list<max_bonds> ls(mol.bonds[i]);
				for (size_t j = 0; j < ls.n; ++j)
					ls[j] += n;
				bonds.push_back(ls);
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
			// update the trace and the sum of the dispersion coefficients matrix
			// dispersion coefficients are calculated from the Lorentz-Berthelot rule
			for (size_t i = 0; i < mol.n; ++i)
			{
				T C6 = 2*lj_halfR[n + i];
				C6 *= C6;
				C6 = C6 * C6 * C6;
				C6 = 2 * lj_sqrteps[n + i]*lj_sqrteps[n + i] * C6;
				tracedisp += C6;
				sumdisp += C6;
			}
			for (size_t i = n; i < n+mol.n; ++i)
				for (size_t j = 0; j < i; ++j)
				{
					T C6 = lj_halfR[i]+lj_halfR[j];
					C6 *= C6;
					C6 = C6 * C6 * C6;
					C6 = 2 * lj_sqrteps[i]*lj_sqrteps[j] * C6;
					sumdisp += 2*C6;
				}
			n += mol.n;
			dof = 3*n;
			first_step = true;
			kin_updated = false;
		}

		T kinetic_energy() const
		{
			if (!kin_updated)
			{
				kin = 0;
				for (std::size_t i = 0; i < n; ++i)
					kin += dot(p[i], p[i]) / m[i];
				kin /= 2;

				kin_updated = true;
			}
			return kin;
		}

		T total_energy() const
		{
			return kinetic_energy() + potential; // in kcal/mol
		}

		T temperature() const
		// instantaneous temperature of the system, in K
		{
			return 2*kinetic_energy() / (dof*kB<T>);
		}

		T kT() const
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
		{
			return side*side*side; // in angstrom^3 = 10^-30 m^3
		}

		T density() const noexcept
		{
			return M / volume(); // in amu/angstrom^3
		}

		T number_density() const noexcept
		{
			return n / volume(); // in angstrom^-3
		}

		T pressure() const
		// instantaneous pressure calculated with instantaneous temperature in kcal/(mol angstrom^3)
		{
			return (2*kinetic_energy() + virial)/(3*volume());
		}

		T pressure_T_ref() const noexcept
		// instantaneous pressure calculated with reference temperature in kcal/(mol angstrom^3)
		{
			return (n*kT_ref() + virial/3) / volume();
		}

		T external_pressure() const
		// instantaneous pressure calculated from external virial in kcal/(mol angstrom^3)
		// not valid if the system is periodic
		{
			return (virial - dot(x, f)) / (3*volume());
		}

		void step(T dt = 1e-3L)
		// perform an integration step of `dt` picoseconds
		{
			if (first_step)
				copy_to_valarray();
			integ.step(*this, T(dt * 20.4548282844073286665866518779L), first_step); // picoseconds to AKMA time unit
			first_step = false;

			// rescale temperatures if asked for
			rescale_temp();
		}

		const state<T, 3>& force(bool eval = true) override
		// calculate all forces
		// if `eval` is false, reuse old values
		{
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

				kin_updated = false;
			}
			return f;
		}

		const state<T, 3>& force_long(bool eval = true)
		// calculate long-range forces only
		// if `eval` is false, reuse old values
		// required by multi_timestep_leapfrog integrator
		{
			if (eval)
			{
				f = 0;

				force_nonbonded();

				kin_updated = false;
			}
			return f;
		}

		const state<T, 3>& force_short(bool eval = true)
		// calculate short-range forces only
		// if `eval` is false, reuse old values
		// required by multi_timestep_leapfrog integrator
		{
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

				kin_updated = false;
			}
			return f;
		}

		const state<T, 3>& vel(bool eval = true) override
		// calculate velocities from momenta
		// if `eval` is false, reuse old values
		{
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

			std::ranges::generate(noise, gen_gaussian);
			for (size_t i = 0; i < n; ++i)
				noise[i] *= sqrt(m[i] * kT_ref());
			gamma = (kT_ref() / D) / m;
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

		void copy_to_valarray()
		// during simulation, a std::valarray because its implementation supports efficient
		// expression templates through operator overloading, making it a more useful choice
		// with respect to std::vector
		{
			x.resize(n); p.resize(n);
			v.resize(n); f.resize(n);
			m.resize(n);
			noise.resize(n); gamma.resize(n);
			copy(begin(x_tmp), end(x_tmp), begin(x));
			copy(begin(p_tmp), end(p_tmp), begin(p));
			copy(begin(v_tmp), end(v_tmp), begin(v));
			copy(begin(f_tmp), end(f_tmp), begin(f));
			copy(begin(m_tmp), end(m_tmp), begin(m));
		}

		void rescale_temp()
		// rescale the temperature of the system to a temperature fixed by the member
		// variable: temperature_ref
		{
			using std::sqrt;

			if (rescale_temperature)
			{
				T temp = temperature();
				if (temp > 0)
					p *= sqrt(temperature_ref / temp);
			}
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
					T diff = d - par.r0;

					r = (2 * par.Kr * diff * d_) * r;
					f[i] -= r;
					f[j] += r;
					potential += par.Kr * diff*diff;
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
						T cosangle = num / den;
						T sinangle = sqrt((1 + cosangle) * (1 - cosangle));
						T angle = atan2(sinangle, cosangle);
						T diff = angle - par.angle0;
						T Kanglesinden = 2 * par.Kangle * diff / (sinangle * den);
						vec3<T> rijt = rij * (num / rij2), rkjt = rkj * (num / rkj2);
						// Urey-Bradley potential
						vec3<T> rik = x[i] - x[k];
						T dik = norm(rik);
						T dik_ = 1/dik;
						T diffUB = dik - par.rUB;

						rik = (2 * par.KUB * diffUB * dik_) * rik;

						f[i] += (rkj - rijt) * Kanglesinden - rik;
						f[j] += (rijt + rkjt - rij - rkj) * Kanglesinden;
						f[k] += (rij - rkjt) * Kanglesinden + rik;
						potential += par.Kangle * diff*diff + par.KUB * diffUB*diffUB;
					}
				}
		}

		void force_dihedrals()
		// (TODO) calculate dihedral angles forces and energies for all atoms
		{

		}

		void force_impropers()
		// (TODO) calculate improper angles forces and energies for all atoms
		{

		}

		void eval_nonbonded(std::size_t i, std::size_t j, T fact = 1) noexcept
		// evaluate non-bonded forces and energies for a pair of atoms with indices `i`, `j`
		// `fact` is a factor which forces and energies are multiplied by
		// this method is called inside `correct_nonbonded`
		{
			T Rij = lj_halfR[i] + lj_halfR[j];

			vec3<T> r = x[i] - x[j];
			T r2_ = 1/dot(r, r);
			T d_ = sqrt(r2_);
			T s = Rij * Rij * r2_;
			s = s * s * s;
			T epsijs = lj_sqrteps[i] * lj_sqrteps[j] * s;
			T coulomb = z[i] * z[j] * d_;
			T vir = (12 * epsijs * (s - 1) + coulomb) * fact;
			vec3<T> fij = (vir * r2_) * r;
			f[i] += fij;
			f[j] -= fij;
			potential += (epsijs * (s - 2) + coulomb) * fact;
			virial += vir;
		}

		void correct_nonbonded()
		// the long-range forces must be excluded for 1-2, 1-3 and 1-4 bonded atoms
		// a correction is applied so that the force calculated through a long-range
		// summation algorithm vanish for these atom pairs
		{
			using std::size_t;
			for (size_t i = 0; i < n; ++i)
				for (size_t mj = 0; mj < bonds[i].n; ++mj)
				{
					const size_t j = bonds[i][mj];
					if (i < j)
						eval_nonbonded(i, j, -1);
					for (size_t mk = 0; mk < bonds[j].n; ++mk)
					{
						const size_t k = bonds[j][mk];
						if (i < k)
							eval_nonbonded(i, k, -1);
					}
				}
		}

		void force_nonbonded()
		// calculate non-bonded forces and energies (both electrostatic and Lennard-Jones)
		// for all atoms
		{
			correct_nonbonded();

			// call long-range summation algorithm operator()
			lrsum(*this, tp);

			// the following applies a long-range correction to the energy due to dispersion interactions
			// since they are truncated at a cutoff radius (see Allen et al., 2017, pp. 79-81)
			T cutoff = lrsum.cutoff_radius();
			if (cutoff)
			{
				T cutoff3 = cutoff*cutoff*cutoff;
				T disp_average = 0;
				if (n > 1)
					disp_average = (sumdisp - tracedisp) / (n * (n-1));
				T factor = 2 * std::numbers::pi_v<T> * n * number_density() * disp_average / cutoff3;
				potential -= factor/3;
				virial -= 2*factor;
			}
		}

		void elastic_confining(vec3<T> k)
		// an elastic potential may be used to confine the atoms (if the system is not periodic)
		{
			using std::size_t;
			f -= k * x;
			for (size_t i = 0; i < n; ++i)
				potential += T(.5L) * dot(k, x[i]*x[i]);
		}

		void diff_box_confining(T steepness)
		// a potential well or box may be approximated by a potential which gives a constant
		// force outside a cube and no force inside it
		// `steepness` is the steepness of the potential outside the cube, which is equal
		// to the force applied
		{
			T hside = side / 2;
			for (std::size_t i = 0; i < n; ++i)
			{
				if (x[i][0] > hside)
				{
					f[i][0] -= steepness;
					potential += (x[i][0] - hside)*steepness;
				}
				if (x[i][1] > hside)
				{
					f[i][1] -= steepness;
					potential += (x[i][1] - hside)*steepness;
				}
				if (x[i][2] > hside)
				{
					f[i][2] -= steepness;
					potential += (x[i][2] - hside)*steepness;
				}

				if (x[i][0] < -hside)
				{
					f[i][0] += steepness;
					potential += (-hside - x[i][0])*steepness;
				}
				if (x[i][1] < -hside)
				{
					f[i][1] += steepness;
					potential += (-hside - x[i][1])*steepness;
				}
				if (x[i][2] < -hside)
				{
					f[i][2] += steepness;
					potential += (-hside - x[i][2])*steepness;
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
			tracedisp = sumdisp = 0;
			n = dof = 0;
			t = 0;
			first_step = true;
			kin_updated = false;
		}

		molecular_system(const molecular_system& other)
		// copy constructor
			: x(other.x), p(other.p), v(other.v), f(other.f), noise(other.noise), m(other.m), gamma(other.gamma),
			x_tmp(other.x_tmp), p_tmp(other.p_tmp), v_tmp(other.v_tmp), f_tmp(other.f_tmp), m_tmp(other.m_tmp),
			z(other.z), lj_sqrteps(other.lj_sqrteps), lj_halfR(other.lj_halfR), id(other.id), bonds(other.bonds), impropers(other.impropers),
			M(other.M), Z(other.Z), Z2(other.Z2), tracedisp(other.tracedisp), sumdisp(other.sumdisp), n(other.n), dof(other.dof),
			mersenne_twister(other.mersenne_twister), side(other.side), temperature_ref(other.temperature_ref),
			t(other.t), potential(other.potential), virial(other.virial), D(other.D), first_step(other.first_step)
		{}

		template <
			template <typename, typename, template <typename...> typename...> typename IntegU,
			template <typename, typename> typename LRSummationU,
			template <typename...> typename ... IntegPars2
		>
		molecular_system(const molecular_system<T, IntegU, LRSummationU, IntegPars2...>& other)
		// template copy constructor
			: x(other.x), p(other.p), v(other.v), f(other.f), noise(other.noise), m(other.m), gamma(other.gamma),
			x_tmp(other.x_tmp), p_tmp(other.p_tmp), v_tmp(other.v_tmp), f_tmp(other.f_tmp), m_tmp(other.m_tmp),
			z(other.z), lj_sqrteps(other.lj_sqrteps), lj_halfR(other.lj_halfR), id(other.id), bonds(other.bonds), impropers(other.impropers),
			M(other.M), Z(other.Z), Z2(other.Z2), tracedisp(other.tracedisp), sumdisp(other.sumdisp), n(other.n), dof(other.dof),
			mersenne_twister(other.mersenne_twister), side(other.side), temperature_ref(other.temperature_ref),
			t(other.t), potential(other.potential), virial(other.virial), D(other.D), first_step(other.first_step)
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
			M = other.M;
			Z = other.Z;
			Z2 = other.Z2;
			tracedisp = other.tracedisp;
			sumdisp = other.sumdisp;
			n = other.n;
			dof = other.dof;
			mersenne_twister = other.mersenne_twister;
			side = other.side;
			temperature_ref = other.temperature_ref;
			t = other.t;
			potential = other.potential;
			virial = other.virial;
			D = other.D;
			first_step = other.first_step;
			return *this;
		}

		template <
			template <typename, typename, template <typename...> typename...> typename IntegU,
			template <typename, typename> typename LRSummationU,
			template <typename...> typename ... IntegPars2
		>
		molecular_system& operator=(const molecular_system<T, IntegU, LRSummationU, IntegPars2...>& other)
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
			M = other.M;
			Z = other.Z;
			Z2 = other.Z2;
			tracedisp = other.tracedisp;
			sumdisp = other.sumdisp;
			n = other.n;
			dof = other.dof;
			mersenne_twister = other.mersenne_twister;
			side = other.side;
			temperature_ref = other.temperature_ref;
			t = other.t;
			potential = other.potential;
			virial = other.virial;
			D = other.D;
			first_step = other.first_step;
			return *this;
		}
	};

} // namespace physics

#endif // PHYSICS_MOLECULE_H
































