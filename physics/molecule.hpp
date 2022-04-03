#ifndef PHYSICS_MOLECULE_H
#define PHYSICS_MOLECULE_H

#include "../math.hpp" // deg2rad

#include "point.hpp"
#include "integrator.hpp"

#include <concepts> // floating_point
#include <valarray>
#include <vector>
#include <map>
#include <tuple>
#include <algorithm> // min, max, ranges::generate
#include <cmath> // sqrt, atan2, pow
#include <random> // mt19937_64, normal_distribution

namespace physics
{
	enum class atom_type : unsigned char { HW, OW, N };

	constexpr std::size_t max_bonds = 4;

	constexpr unsigned char atom_number[] = {1, 8};

	template <std::floating_point T>
	constexpr T kC = 332.06371330062682040436798930743L; // Coulomb constant in AKMA units

	template <std::floating_point T>
	constexpr T kB = 0.00198720425864083173996175908222L; // Boltzmann constant in AKMA units

	template <std::floating_point T>
	constexpr T DW5 = 0.00616480364668354400695529510607L; // water self-diffusion constant in AKMA units at 5 °C

	template <std::floating_point T>
	constexpr T DW25 = 0.01123940014569822971609058164065L; // water self-diffusion constant in AKMA units at 25 °C

	template <std::floating_point T>
	constexpr T atom_mass[] = {1.008L, 15.9994L};

	/*template <std::floating_point T>
	std::map<const std::pair<atom_type, atom_type>, const std::pair<T, T>> bond_params = { { {atom_type::HW, atom_type::OW}, {450, 0.9572L} } }; // k, r0

	template <std::floating_point T>
	std::map<const std::tuple<atom_type, atom_type, atom_type>, const std::pair<T, T>> angle_params =
		{ { {atom_type::HW, atom_type::OW, atom_type::HW}, {55, math::deg2rad(104.52L)} } }; // k, theta0

	template <std::floating_point T>
	std::pair<T, T> LJ_params[][int(atom_type::N)] = { { {0.046L, 0.4L}, {0.0836L, 1.7753L} },
													   { {0.0836L, 1.7753L}, {0.1521L, 3.1507L} } }; // epsilon, sigma*/

	template <std::floating_point T>
	std::map<const std::pair<atom_type, atom_type>, const std::pair<T, T>> bond_params = { { {atom_type::HW, atom_type::OW}, {358.50860420650095602294455066922L, 1.027L} } }; // k, r0 // 3000 kJ / mol A^2

	template <std::floating_point T>
	std::map<const std::tuple<atom_type, atom_type, atom_type>, const std::pair<T, T>> angle_params =
		{ { {atom_type::HW, atom_type::OW, atom_type::HW}, {45.769598470363288718929254302103L, math::deg2rad(114.70L)} } }; // k, theta0 // 383 kJ/mol

	template <std::floating_point T>
	std::pair<T, T> LJ_params[][int(atom_type::N)] = { { {0, 0}, {0.0836L, 1.7753L} },
													   { {0.0836L, 1.7753L}, {0.18936998087954110898661567877629L, 3.1776L} } }; // epsilon, sigma

	template <std::floating_point T>
	std::pair<T, T> LJ_params2[][int(atom_type::N)] = { { {0, 0}, {0, 0} },
														{ {0, 0}, {0.18936998087954110898661567877629L, 3.1776L} } }; // epsilon, sigma

	template <std::floating_point T>
	struct molecule
	{
		state<T, 3> x; // position
		std::vector<atom_type> id; // identity of the atom
		std::vector<T> part_q; // partial charges
		std::vector<fixed_list<max_bonds>> bonds; // as an adjacency list for each atom (max 4 bonds per atom)
		std::vector<std::array<unsigned int, 4>> impropers; // list of impropers in the molecule
		unsigned int n; // number of atoms in the molecule
	};

	template <std::floating_point T>
	const molecule<T> water_tip3p // TIP3P model
	{
		.x = { {-0.75695032726366116174L, 0.58588227661829493656L}, {0}, {0.75695032726366116174L, 0.58588227661829493656L} },
			// oxygen at the origin
		.id = {atom_type::HW, atom_type::OW, atom_type::HW},
		.part_q = {0.417L, -0.834L, 0.417L},
		.bonds = { {1}, {0, 2}, {1} },
		.impropers = {}, // no impropers
		.n = 3
	};

	template <std::floating_point T>
	const molecule<T> water_fba_eps // FBA/eps model
	{
		.x = { {-0.75695032726366116174L, 0.58588227661829493656L}, {0}, {0.75695032726366116174L, 0.58588227661829493656L} },
			// oxygen at the origin
		.id = {atom_type::HW, atom_type::OW, atom_type::HW},
		.part_q = {0.4225L, -0.8450L, 0.4225L},
		.bonds = { {1}, {0, 2}, {1} },
		.impropers = {}, // no impropers
		.n = 3
	};

	template <std::floating_point T, integrator Integ = leapfrog<T>>
	class molecular_system : public physical_system_base<T>
	{
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
					r = (par.first * (d - par.second) / d) * r;
					a[i] -= r;
					a[j] += r;
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
						T Kanglesinden = par.first * (angle - par.second) / (sinangle * den);
						point3<T> rijt = rij * (num / rij2), rkjt = rkj * (num / rkj2);
						/* // for Urey-Bradley
						point3<T> rik = x[i] - x[k];
						T dik = norm(rik);
						rik /= dik;
						rik = (get<2>(par) * (dik - get<3>(par))) * rik;
						*/
						a[i] += (rkj - rijt) * Kanglesinden; // - rik
						a[j] += (rijt + rkjt - rij - rkj) * Kanglesinden;
						a[k] += (rij - rkjt) * Kanglesinden; // + rik
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
							for (size_t ml = 0; ml < f_bonds[k].n; ++ml)
							{
								size_t l = f_bonds[k][ml];
								atom_type idl = id[l];
								if (i0 == l)
									if (get<0>(dihedral_params<T>[{ idi, idj, idj, idl }]) != 0)
										return true;
									else
										return false;
							}*/
					}
			}
			return false;
		}

		void force_nonbonded()
		{
			using std::size_t;

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
						r = ((par.first * s * (s - 1) + kC<T> * part_q[i] * part_q[j] * d_) * r2_) * r;
						a[i] += r;
						a[j] -= r;
					}
		}

		void elastic_confining(point3<T> k)
		{
			a -= k * x;
		}

		void diff_box_confining(T steepness)
		{
			T hside = side / 2;
			for (std::size_t i = 0; i < n; ++i)
			{
				if (x[i][0] > hside)
					a[i][0] -= steepness;
				if (x[i][1] > hside)
					a[i][1] -= steepness;
				if (x[i][2] > hside)
					a[i][2] -= steepness;

				if (x[i][0] < -hside)
					a[i][0] += steepness;
				if (x[i][1] < -hside)
					a[i][1] += steepness;
				if (x[i][2] < -hside)
					a[i][2] += steepness;
			}
		}

		T calculate_temperature()
		{
			T two_kin = 0;
			for (std::size_t i = 0; i < n; ++i)
				two_kin += atom_mass<T>[ int(id[i]) ] * dot(v[i], v[i]);
			return two_kin / (3*n*kB<T>);
		}

		public:

		state<T, 3> x, v, a; // position, velocity, acceleration
		state<T, 3> noise;
		std::valarray<T> gamma;
		std::vector<atom_type> id; // identity of the atom
		std::vector<T> part_q; // partial charges
		std::vector<fixed_list<max_bonds>> bonds; // bonds as an adjacency list
		std::vector<std::array<unsigned int, 4>> impropers; // list of impropers in the whole system
		unsigned int n; // total number of atoms
		Integ integ; // integrator
		std::mt19937_64 mersenne_twister;
		std::normal_distribution<T> n_dist;
		T side, temperature, D; // default = 25 °C
		T t; // time

		static const std::size_t max_atoms = 2'000'000;

		molecular_system(T side = 40, T temp = 298.15, T D = DW25<T>, Integ integ = Integ())
			: n(0), integ(integ), mersenne_twister(0), n_dist(0, 1), side(side), temperature(temp), D(D), t(0)
		{
			using std::size_t;
			for (size_t i = 0; i < size_t(atom_type::N); ++i)
				for (size_t j = 0; j < size_t(atom_type::N); ++j)
				{
					LJ_params<T>[i][j].first *= 12;
					LJ_params<T>[i][j].second *= std::pow(.5L, 1.L/6);
				}
		}

		void add_molecule(const molecule<T>& mol, const point3<T>& pos = 0)
		{
			using std::size_t;

			if (n + mol.n > max_atoms)
				throw("Error: Exceeded maximum number of atoms");
			state<T, 3> x_tmp = x, v_tmp = v, a_tmp = a;
			x.resize(n + mol.n);
			v.resize(n + mol.n);
			a.resize(n + mol.n);
			noise.resize(n + mol.n);
			gamma.resize(n + mol.n);
			for (size_t i = 0; i < n; ++i)
				x[i] = x_tmp[i];
			for (size_t i = 0; i < n; ++i)
				v[i] = v_tmp[i];
			for (size_t i = 0; i < n; ++i)
				a[i] = a_tmp[i];

			for (size_t i = 0; i < mol.n; ++i)
				x[n + i] = mol.x[i] + pos;
			for (size_t i = 0; i < mol.n; ++i)
				v[n + i] = 0;
			for (size_t i = 0; i < mol.n; ++i)
				id.push_back(mol.id[i]);
			for (size_t i = 0; i < mol.n; ++i)
				part_q.push_back(mol.part_q[i]);

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
		}

		void step(T dt = 1e-3L, bool temp_rescale = false) // picoseconds
		{
			using std::sqrt;

			integ.step(*this, T(dt * 20.4548282844073286665866518779L)); // picoseconds to AKMA time unit

			if (temp_rescale)
			{
				T measured_temp = calculate_temperature();
				if (measured_temp > 0)
					v *= sqrt(temperature / measured_temp);
			}
		}

		void accel()
		{
			using std::size_t;

			a = 0;

			force_bonds();
			force_angles();
			
			a *= 2;
			// the 2 factor is due to the fact that the K constant parameters are two times the relative elastic constants

			force_nonbonded();

			diff_box_confining(10);

			for (size_t i = 0; i < n; ++i)
				a[i] /= atom_mass<T>[ int(id[i]) ];
		}

		void rand()
		{
			using std::sqrt;

			auto gen = [this]()
			{
				return point3<T>(n_dist(mersenne_twister),
								 n_dist(mersenne_twister),
								 n_dist(mersenne_twister));
			};
			std::ranges::generate(noise, gen);
			for (size_t i = 0; i < n; ++i)
				noise[i] *= sqrt(kB<T> * temperature / atom_mass<T>[ int(id[i]) ]);
			for (size_t i = 0; i < n; ++i)
				gamma[i] = kB<T> * temperature / (atom_mass<T>[ int(id[i]) ] * D);
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
	};

} // namespace physics

#endif // PHYSICS_MOLECULE_H
































