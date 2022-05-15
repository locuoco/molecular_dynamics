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
#include <cmath> // sqrt, atan2, pow
#include <random> // mt19937_64, normal_distribution
#include <thread>

#include "../math.hpp" // deg2rad

#include "point.hpp"
#include "integrator.hpp"

namespace physics
{
	enum class atom_type : unsigned char { HT, OT, HTL, OTL, HF, OF, N };

	constexpr std::size_t max_bonds = 4; // max coordination number assumed

	constexpr unsigned char atom_number[] = {1, 8, 1, 8, 1, 8};

	template <std::floating_point T>
	constexpr T kC = 332.06371330062682040436798930743L; // Coulomb constant in AKMA units

	template <std::floating_point T>
	constexpr T kB = 0.00198720425864083173996175908222L; // Boltzmann constant in AKMA units

	template <std::floating_point T = double>
	constexpr T DW5 = 0.00616480364668354400695529510607L; // water self-diffusion constant in AKMA units at 5 °C

	template <std::floating_point T = double>
	constexpr T DW25 = 0.01123940014569822971609058164065L; // water self-diffusion constant in AKMA units at 25 °C

	template <std::floating_point T>
	constexpr T atom_mass[] = {1.008L, 15.9994L, 1.008L, 15.9994L, 1.008L, 15.9994L};

	template <std::floating_point T>
	std::map<const std::pair<atom_type, atom_type>, const std::pair<T, T>> bond_params = // k, r0
	{
		{ {atom_type::HT, atom_type::OT}, {450, 0.9572L} },
		{ {atom_type::HTL, atom_type::OTL}, {450, 0.9572L} },
		{ {atom_type::HF, atom_type::OF}, {358.50860420650095602294455066922L, 1.027L} }
	};

	template <std::floating_point T>
	std::map<const std::tuple<atom_type, atom_type, atom_type>, const std::pair<T, T>> angle_params = // k, theta0 // 383 kJ/mol
	{
		{ {atom_type::HT, atom_type::OT, atom_type::HT}, {55, math::deg2rad(104.52L)} },
		{ {atom_type::HTL, atom_type::OTL, atom_type::HTL}, {55, math::deg2rad(104.52L)} },
		{ {atom_type::HF, atom_type::OF, atom_type::HF}, {45.769598470363288718929254302103L, math::deg2rad(114.70L)} }
	};

	long double two1_6 = 1.1224620483093729814335330496792L;

	template <std::floating_point T> // epsilon, sigma
	std::pair<T, T> LJ_params[][int(atom_type::N)] =
	{
		{ {0.046L, 0.4L*two1_6}, {0.0836L, 1.7753L*two1_6}, {0, 0}, {0, 0}, {0, 0}, {0, 0} }, // HT
		{ {0.0836L, 1.7753L*two1_6}, {0.1521L, 3.1507L*two1_6}, {0, 0}, {0, 0}, {0, 0}, {0, 0} }, // OT
		{ {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0} }, // HTL
		{ {0, 0}, {0, 0}, {0, 0}, {0.102L, 3.188L*two1_6}, {0, 0}, {0, 0} }, // OTL
		{ {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0} }, // HF
		{ {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0.18936998087954110898661567877629L, 3.1776L*two1_6} } // OF
	};

	template <std::floating_point T = double>
	struct molecule
	{
		state<T, 3> x; // position
		std::vector<atom_type> id; // identity of the atom
		std::vector<T> part_q; // partial charges
		std::vector<fixed_list<max_bonds>> bonds; // as an adjacency list for each atom
		std::vector<std::array<unsigned int, 4>> impropers; // list of impropers in the molecule
		unsigned int n; // number of atoms in the molecule
	};

	template <std::floating_point T = double>
	const molecule<T> water_tip3p // TIP3P model
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
	const molecule<T> water_tip3p_lr // TIP3P model long-range
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
	const molecule<T> water_fba_eps // FBA/eps model
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
		template <typename, typename, template <typename...> typename...> typename IntegT = leapfrog,
		template <typename...> typename ... TPars
	>
	requires integrator<IntegT<T, state<T, 3>, TPars...>, T, state<T, 3>>
	class molecular_system : public physical_system_base<T, state<T, 3>>
	{
		using Integ = IntegT<T, state<T, 3>, TPars...>;

		std::normal_distribution<T> n_dist;
		mutable T kin;
		mutable bool kin_updated;

		public:

		using scalar_type = T;

		state<T, 3> x, p, v, f; // position, momentum, velocity, force
		state<T, 3> noise; // gaussian noise
		std::valarray<T> m, gamma; // mass, damping factor
		std::vector<T> part_q; // partial charges
		std::vector<atom_type> id; // identity of the atom
		std::vector<fixed_list<max_bonds>> bonds; // bonds as an adjacency list
		std::vector<std::array<unsigned int, 4>> impropers; // list of impropers in the whole system
		unsigned int n, dof; // total number of atoms, degrees of freedom
		Integ integ; // integrator
		std::mt19937_64 mersenne_twister;
		T side, temperature_fixed, D; // D = diffusion coefficient
		// temperature_fixed is a fixed number independent of the simulation
		T t, U; // time, kinetic energy, potential energy
		bool rescale_temperature, first_step;

		static const std::size_t max_atoms = 2'000'000;

		molecular_system(T side = 100, T temp = 298.15, T D = DW25<T>, Integ integ = Integ())
			: n_dist(0, 1), kin_updated(false), n(0), integ(integ), mersenne_twister(0), side(side), temperature_fixed(temp), D(D), t(0),
			  rescale_temperature(false), first_step(true)
		{}

		void add_molecule(const molecule<T>& mol, const point3<T>& pos = 0)
		{
			using std::size_t;
			using std::sqrt;

			if (n + mol.n > max_atoms)
				throw("Error: Exceeded maximum number of atoms");
			state<T, 3> x_tmp = x, p_tmp = p, v_tmp = v, f_tmp = f;
			x.resize(n + mol.n);
			p.resize(n + mol.n);
			v.resize(n + mol.n);
			f.resize(n + mol.n);
			m.resize(n + mol.n);
			noise.resize(n + mol.n);
			gamma.resize(n + mol.n);
			for (size_t i = 0; i < n; ++i)
				x[i] = x_tmp[i];
			for (size_t i = 0; i < n; ++i)
				p[i] = p_tmp[i];
			for (size_t i = 0; i < n; ++i)
				v[i] = v_tmp[i];
			for (size_t i = 0; i < n; ++i)
				f[i] = f_tmp[i];

			for (size_t i = 0; i < mol.n; ++i)
				x[n + i] = mol.x[i] + pos;
			for (size_t i = 0; i < mol.n; ++i)
				id.push_back(mol.id[i]);
			for (size_t i = 0; i < mol.n; ++i)
				part_q.push_back(mol.part_q[i]);
			for (size_t i = 0; i < n + mol.n; ++i)
				m[i] = atom_mass<T>[ int(id[i]) ];
			for (size_t i = 0; i < mol.n; ++i)
				p[n + i] = gen_gaussian() * sqrt(m[n + i] * kB<T> * temperature_fixed);

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
			n += mol.n;
			dof = 3*n;
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
			return kinetic_energy() + U;
		}

		T temperature() const
		{
			return 2*kinetic_energy() / (3*n*kB<T>);
		}

		T kT_fixed() const
		{
			return kB<T> * temperature_fixed;
		}

		void step(T dt = 1e-3L) // picoseconds
		{
			integ.step(*this, T(dt * 20.4548282844073286665866518779L), first_step); // picoseconds to AKMA time unit
			first_step = false;

			rescale_temp();
		}

		const state<T, 3>& force(bool eval = true)
		{
			if (eval)
			{
				f = 0;
				U = 0;

				force_nonbonded(8);

				force_dihedrals();
				force_impropers();

				force_angles();
				force_bonds();

				diff_box_confining(2);
			}
			kin_updated = false;
			return f;
		}

		const state<T, 3>& force_long(bool eval = true)
		{
			if (eval)
			{
				f = 0;

				force_nonbonded(8);
			}
			kin_updated = false;
			return f;
		}

		const state<T, 3>& force_short(bool eval = true)
		{
			if (eval)
			{
				f = 0;
				U = 0;

				force_dihedrals();
				force_impropers();

				force_angles();
				force_bonds();

				diff_box_confining(10);
			}
			kin_updated = false;
			return f;
		}

		const state<T, 3>& vel(bool eval = true)
		{
			if (eval)
				for (std::size_t i = 0; i < n; ++i)
					v[i] = p[i] / m[i];
			return v;
		}

		void rand()
		{
			using std::sqrt;

			std::ranges::generate(noise, gen_gaussian);
			for (size_t i = 0; i < n; ++i)
				noise[i] *= sqrt(m[i] * kB<T> * temperature_fixed);
			gamma = (kB<T> * temperature_fixed / D) / m;
			/*for (size_t i = 0; i < n; ++i)
			{
				T minx = x[i][0], maxx = x[i][0];
				minx = std::min(minx, x[i][1]);
				maxx = std::max(maxx, x[i][1]);
				minx = std::min(minx, x[i][2]);
				maxx = std::max(maxx, x[i][2]);
				if (minx >= -side/2 && maxx <= side/2)
				{
					noise[i] = 0;
					gamma[i] = 0;
				}
			}*/
		}

		private:

		void rescale_temp()
		{
			using std::sqrt;

			if (rescale_temperature)
			{
				T temp = temperature();
				if (temp > 0)
					p *= sqrt(temperature_fixed / temp);
			}
		}

		point3<T> gen_gaussian()
		{
			return point3<T>(
				n_dist(mersenne_twister),
				n_dist(mersenne_twister),
				n_dist(mersenne_twister)
			);
		}

		void force_bonds()
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

					point3<T> r = x[i] - x[j];
					T d = norm(r);
					T diff = d - par.second;
					r = (2 * par.first * diff / d) * r;
					f[i] -= r;
					f[j] += r;
					U += par.first * diff*diff;
				}
		}

		void force_angles()
		{
			using std::sqrt;
			using std::atan2;
			using std::size_t;

			for (size_t i = 0; i < n; ++i)
				for (size_t l = 0; l < bonds[i].n; ++l)
				{
					size_t j = bonds[i][l];
					atom_type idj = id[j];

					point3<T> rij = x[i] - x[j];
					T rij2 = dot(rij, rij);
					for (size_t m = 0; m < bonds[j].n; ++m)
					{
						size_t k = bonds[j][m];
						if (!(i < k))
							continue;
						atom_type idi = id[i], idk = id[k];
						if (idi > idk)
							std::swap(idi, idk);
						const auto& par = angle_params<T>[{ idi, idj, idk }];

						point3<T> rkj = x[k] - x[j];
						T rkj2 = dot(rkj, rkj);
						T num = dot(rij, rkj), den = sqrt(rij2 * rkj2);
						T cosangle = num / den;
						T sinangle = sqrt((1 + cosangle) * (1 - cosangle));
						T angle = atan2(sinangle, cosangle);
						T diff = angle - par.second;
						T Kanglesinden = 2 * par.first * diff / (sinangle * den);
						point3<T> rijt = rij * (num / rij2), rkjt = rkj * (num / rkj2);
						/* // for Urey-Bradley
						point3<T> rik = x[i] - x[k];
						T dik = norm(rik);
						rik /= dik;
						rik = (get<2>(par) * (dik - get<3>(par))) * rik;
						*/
						f[i] += (rkj - rijt) * Kanglesinden; // - rik
						f[j] += (rijt + rkjt - rij - rkj) * Kanglesinden;
						f[k] += (rij - rkjt) * Kanglesinden; // + rik
						U += par.first * diff*diff;
					}
				}
		}

		void force_dihedrals()
		{
			
		}

		void force_impropers()
		{
			using std::size_t;

			/*for (size_t i = 0; i < impropers.size(); ++i)
			{
				atom_type idi = impropers[i][0], idj = impropers[i][1],
						  idk = impropers[i][2], idl = impropers[i][3];
				const auto& par = improper_params<T>[{ idi, idj, idk, idl }];
				// TODO
			}*/
		}

		bool check_vicinity(std::size_t i, std::size_t i0)
		{
			using std::size_t;

			if (i0 == i)
				return true;
			//atom_type idi = id[i];
			for (size_t mj = 0; mj < bonds[i].n; ++mj)
			{
				size_t j = bonds[i][mj];
				//atom_type idj = id[j];
				if (i0 == j)
					return true;
				else
					for (size_t mk = 0; mk < bonds[j].n; ++mk)
					{
						size_t k = bonds[j][mk];
						//atom_type idk = id[k];
						if (i0 == k)
							return true;
						/*else
							for (size_t ml = 0; ml < bonds[k].n; ++ml)
							{
								size_t l = bonds[k][ml];
								atom_type idl = id[l];
								if (i0 == l)
									if (get<0>(dihedral_params<T>[{ idi, idj, idk, idl }]) != 0)
										return true;
									else
										return false;
							}*/
					}
			}
			return false;
		}

		void force_nonbonded(size_t num_threads = 1)
		{
			using std::size_t;

			if (num_threads == 1)
			{
				for (size_t i = 0; i < n; ++i)
					for (size_t j = i+1; j < n; ++j)
						if (!check_vicinity(i, j))
						{
							atom_type idi = id[i], idj = id[j];
							const auto& par = LJ_params<T>[int(idi)][int(idj)];

							point3<T> r = x[i] - x[j];
							T r2_ = 1/dot(r, r);
							T d_ = sqrt(r2_);
							T s = par.second * par.second * r2_;
							s = s * s * s;
							T coulomb = kC<T> * part_q[i] * part_q[j] * d_;
							r = ((12 * par.first * s * (s - 1) + coulomb) * r2_) * r;
							f[i] += r;
							f[j] -= r;
							U += par.first * s * (s - 2) + coulomb;
						}
			}
			else
			{
				std::vector<T> partU(num_threads);
				auto eval = [this, num_threads, &partU](size_t idx)
				{
					T u = 0;
					size_t block = (n-1)/num_threads + 1;
					for (size_t i = idx*block; (i < n) && (i < (idx+1)*block); ++i)
						for (size_t j = 0; j < n; ++j)
							if (!check_vicinity(i, j))
							{
								atom_type idi = id[i], idj = id[j];
								const auto& par = LJ_params<T>[int(idi)][int(idj)];

								point3<T> r = x[i] - x[j];
								T r2_ = 1/dot(r, r);
								T d_ = sqrt(r2_);
								T s = par.second * par.second * r2_;
								s = s * s * s;
								T coulomb = kC<T> * part_q[i] * part_q[j] * d_;
								r = ((12 * par.first * s * (s - 1) + coulomb) * r2_) * r;
								f[i] += r;
								u += par.first * s * (s - 2) + coulomb;
							}
					partU[idx] = u/2;
				};
				std::vector<std::thread> threads;
				for (size_t i = 0; i < num_threads; ++i)
					threads.push_back(std::move(std::thread(eval, i)));
				for (size_t i = 0; i < num_threads; ++i)
					threads[i].join();
				for (size_t i = 0; i < num_threads; ++i)
					U += partU[i];
			}
		}

		void elastic_confining(point3<T> k)
		{
			using std::size_t;
			f -= k * x;
			for (size_t i = 0; i < n; ++i)
				U += T(.5L) * dot(k, x[i]*x[i]);
		}

		void diff_box_confining(T steepness)
		{
			T hside = side / 2;
			for (std::size_t i = 0; i < n; ++i)
			{
				if (x[i][0] > hside)
				{
					f[i][0] -= steepness;
					U += (x[i][0] - hside)*steepness;
				}
				if (x[i][1] > hside)
				{
					f[i][1] -= steepness;
					U += (x[i][1] - hside)*steepness;
				}
				if (x[i][2] > hside)
				{
					f[i][2] -= steepness;
					U += (x[i][2] - hside)*steepness;
				}

				if (x[i][0] < -hside)
				{
					f[i][0] += steepness;
					U += (-hside - x[i][0])*steepness;
				}
				if (x[i][1] < -hside)
				{
					f[i][1] += steepness;
					U += (-hside - x[i][1])*steepness;
				}
				if (x[i][2] < -hside)
				{
					f[i][2] += steepness;
					U += (-hside - x[i][2])*steepness;
				}
			}
		}
	};

} // namespace physics

#endif // PHYSICS_MOLECULE_H
































