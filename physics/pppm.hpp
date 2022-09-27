//  Particle-particle, particle-mesh method
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

#ifndef PHYSICS_PPPM_H
#define PHYSICS_PPPM_H

#include <iostream> // clog
#include <cmath> // round, ceil, log2, sqrt, sin, cos, abs, remainder
#include <numeric> // iota
#include <numbers> // numbers::sqrt2_v, numbers::sqrt3_v, numbers::pi_v
#include <functional> // function
#include <algorithm> // fill, copy
#include <complex>
#include <stdexcept> // invalid_argument
#include <random> // mt19937_64, uniform_real_distribution

#include "ewald.hpp" // from math/helper.hpp: fasterfc, fastexp, mod, two1_6, is_even
	// eval_ewald_r, eval_ewald_scd, 
#include "../math/dft.hpp"
#include "../utils/parallel_sort.hpp"

namespace physics
{
	static constexpr std::size_t pppm_min_order = 1, pppm_max_order = 7;

	template <typename T>
	constexpr T charge_assignment_function_coeff[pppm_max_order][pppm_max_order][pppm_max_order] =
	// Coefficients used inside `charge_assignment_function`
	// See Appendix E of Deserno et al., I, 1998
		{
			// order = 1
			{{1}},
			// order = 2
			{
				{1.L/2, -1},
				{1.L/2, +1},
			},
			// order = 3
			{
				{1.L/8, -1.L/2, 1.L/2},
				{3.L/4,  0,    -1    },
				{1.L/8, +1.L/2, 1.L/2},
			},
			// order = 4
			{
				{ 1.L/48, -1.L/8,  1.L/4, -1.L/6},
				{23.L/48, -5.L/8, -1.L/4, +1.L/2},
				{23.L/48, +5.L/8, -1.L/4, -1.L/2},
				{ 1.L/48, +1.L/8,  1.L/4, +1.L/6},
			},
			// order = 5
			{
				{  1.L/384, -1.L/48, 1.L/16, -1.L/12, 1.L/24},
				{ 19.L/96, -11.L/24, 1.L/4,  +1.L/6, -1.L/6 },
				{115.L/192,  0,     -5.L/8,   0,      1.L/4 },
				{ 19.L/96, +11.L/24, 1.L/4,  -1.L/6, -1.L/6 },
				{  1.L/384, +1.L/48, 1.L/16, +1.L/12, 1.L/24},
			},
			// order = 6
			{
				{  1.L/3840,  -1.L/384,   1.L/96, -1.L/48,  1.L/48, -1.L/120},
				{237.L/3840, -25.L/128,  21.L/96, -1.L/16, -1.L/16, +1.L/24 },
				{841.L/1920, -77.L/192, -11.L/48, +7.L/24,  1.L/24, -1.L/12 },
				{841.L/1920, +77.L/192, -11.L/48, -7.L/24,  1.L/24, +1.L/12 },
				{237.L/3840, +25.L/128,  21.L/96, +1.L/16, -1.L/16, -1.L/24 },
				{  1.L/3840,  +1.L/384,   1.L/96, +1.L/48,  1.L/48, +1.L/120},
			},
			// order = 7
			{
				{    1.L/46080,   -1.L/3840, 1.L/768,  -1.L/288,   1.L/192, -1.L/240, 1.L/720},
				{  361.L/23040,  -59.L/960, 37.L/384,  -5.L/72,    1.L/96,  +1.L/60, -1.L/120},
				{10543.L/46080, -289.L/768, 79.L/768, +43.L/288, -17.L/192, -1.L/48,  1.L/48 },
				{ 5887.L/11520,    0,      -77.L/192,   0,        21.L/144,  0,      -1.L/36 },
				{10543.L/46080, +289.L/768, 79.L/768, -43.L/288, -17.L/192, +1.L/48,  1.L/48 },
				{  361.L/23040,  +59.L/960, 37.L/384,  +5.L/72,    1.L/96,  -1.L/60, -1.L/120},
				{    1.L/46080,   +1.L/3840, 1.L/768,  +1.L/288,   1.L/192, +1.L/240, 1.L/720},
			},
		};

	template <typename T>
	constexpr T charge_assignment_function(T x, std::size_t pos, std::size_t order) noexcept
	// Calculate the charge assignment (interpolation) function.
	// Let "xbar" be the nearest mesh point if `order` is odd or w.r.t. the midpoint between
	// the two nearest mesh points if `order` is even.
	// `x` is the relative position of the charge w.r.t. "xbar", divided by the cell size.
	// `pos` is the relative displacement of the mesh interpolation point w.r.t.
	// "xbar" - floor((order-1)/2), so that 0 <= pos < order.
	// `order` is the interpolation order, between 1 and 7 included.
	{
		T res = 0;
		for (std::ptrdiff_t i = order-1; i >= 0; --i)
			res = res * x + charge_assignment_function_coeff<T>[order-1][pos][i];
		return res;
	}

	template <typename T>
	constexpr T charge_assignment_function_der(T x, std::size_t pos, std::size_t order) noexcept
	// Calculate the derivative of the charge assignment (interpolation) function.
	// The arguments have the same meaning as in `charge_assignment_function`.
	{
		T res = 0;
		for (std::ptrdiff_t i = order-1; i >= 1; --i)
			res = res * x + i*charge_assignment_function_coeff<T>[order-1][pos][i];
		return res;
	}

	template <typename T>
	const std::function<T(T)> influence_function_den[pppm_max_order] =
	// Calculate the denominator of the influence function (or lattice Green function) for orders between 1 and 7.
		{
			// order = 1
			[](T) { return T(1); },
			// order = 2
			[](T s2) { return 1 - T(2.L/3)*s2; },
			// order = 3
			[](T s2) { return 1 - s2 + T(2.L/15)*s2*s2; },
			// order = 4
			[](T s2) { T s4(s2*s2); return 1 - T(4.L/3)*s2 + T(2.L/5)*s4 - T(4.L/315)*s4*s2; },
			// order = 5
			[](T s2) { T s4(s2*s2); return 1 - T(5.L/3)*s2 + T(7.L/9)*s4 - T(17.L/189)*s4*s2 + T(2.L/2835)*s4*s4; },
			// order = 6
			[](T s2)
			{
				T s4(s2*s2), s8(s4*s4);
				return 1 - 2*s2 + T(19.L/15)*s4 - T(256.L/945)*s4*s2 + T(62.L/4725)*s8 - T(4.L/155'925)*s8*s2;
			},
			// order = 7
			[](T s2)
			{
				T s4(s2*s2), s8(s4*s4);
				return 1 - T(7.L/3)*s2 + T(28.L/15)*s4 - T(16.L/27)*s4*s2 + T(26.L/405)*s8 - T(2.L/1485)*s8*s2 + T(4.L/6'081'075)*s8*s4;
			},
		};

	template <typename T>
	constexpr T error_expansion_coeff[pppm_max_order][pppm_max_order] =
	// Expansion coefficients for the calculation of force RMS error estimation
	// See Appendix of Deserno et al., II, 1998
		{
			// order = 1
			{2.L/3},
			// order = 2
			{1.L/50, 5.L/294},
			// order = 3
			{1.L/588, 7.L/1440, 21.L/3872},
			// order = 4
			{1.L/4320, 3.L/1936, 7601.L/2'271'360, 143.L/28800},
			// order = 5
			{1.L/23232, 7601.L/13'628'160, 143.L/69120, 517'231.L/106'536'960, 106'640'677.L/11'737'571'328},
			// order = 6
			{
				691.L/68'140'800, 13.L/57600, 47021.L/35'512'320, 9'694'607.L/2'095'994'880, 733'191'589.L/59'609'088'000,
				326'190'917.L/11'700'633'600
			},
			// order = 7
			{
				1.L/345'600, 3617.L/35'512'320, 745'739.L/838'397'952, 56'399'353.L/12'773'376'000, 25'091'609.L/1'560'084'480,
				1'755'948'832'039.L/36'229'939'200'000, 4'887'769'399.L/37'838'389'248
			},
		};

	template <typename T, typename State>
	struct pppm
	// Particle-particle, particle-mesh method (also known as PPPM or P^3M or P3M method)
	// Based on Ewald summation, it has O(N log N) asymptotic complexity thanks to FFT.
	// See Deserno et al., 1998.
	{
		void cutoff_radius(T cutoff)
		// Set the cutoff radius (both for PPPM and Lennard-Jones potential truncation).
		// Throw a `std::invalid_argument` if `cutoff` is not strictly positive.
		{
			if (this->cutoff != cutoff)
			{
				if (cutoff <= 0)
					throw std::invalid_argument("Error: cutoff radius must be positive");
				this->cutoff = cutoff;
				update = true;
			}
		}
		T cutoff_radius() const noexcept
		// Return the cutoff radius (both for PPPM and Lennard-Jones potential truncation).
		{
			return cutoff;
		}

		T dielectric() const noexcept
		// Return the external relative dielectric constant (1 by default)
		{
			return dielec;
		}
		void dielectric(T dielec) noexcept
		// Set the external relative dielectric constant (1 by default)
		{
			this->dielec = dielec;
		}

		void charge_assignment_order(std::size_t order)
		// Set the charge assignment interpolation order.
		// Throw a `std::invalid_argument` if `order` is smaller than `pppm_min_order` or
		// if it is greater than `pppm_max_order`.
		{
			if (this->order != order)
			{
				if (order < pppm_min_order)
					throw std::invalid_argument(std::string("Error: minimum PPPM order is ") + std::to_string(pppm_min_order));
				if (order > pppm_max_order)
					throw std::invalid_argument(std::string("Error: maximum PPPM order is ") + std::to_string(pppm_max_order));
				this->order = order;
				update = true;
			}
		}
		T charge_assignment_order() const noexcept
		// Return the charge assignment interpolation order.
		{
			return order;
		}

		void cell_multiplier(std::ptrdiff_t ch_num_cells) noexcept
		// Set the "cell multiplier". It can be used to increase or decrease the refinement of
		// the mesh, following a logarithmic scale. Default value is 0.
		{
			this->ch_num_cells = ch_num_cells;
		}
		std::ptrdiff_t cell_multiplier() const noexcept
		// Return the "cell multiplier".
		{
			return ch_num_cells;
		}

		void set_diff_scheme(const std::string& scheme)
		// Set differentiation scheme. Two differentiation schemes are supported:
		// "ik": ik-differentiation (gradient calculated in Fourier space).
		// "ad": analytic differentiation (gradient calculated in real space).
		// Throw a `std::invalid_argument` if `scheme` is neither "ik" nor "ad".
		// Note that std::string::operator== is case sensitive.
		{
			if (scheme == std::string("ik"))
			{
				if (!use_ik)
				{
					use_ik = true;
					update = true;
				}
			}
			else if (scheme == std::string("ad"))
			{
				if (use_ik)
				{
					use_ik = false;
					update = true;
				}
			}
			else
				throw std::invalid_argument(std::string("Error: Unknown differentiation scheme: ") + scheme);
		}

		void ewald_par(T kappa)
		// Set the Ewald parameter manually. If this method is not called, the (approximate) optimal Ewald
		// parameter is determined automatically by minimizing the force RMS estimation.
		// Throw a `std::invalid_argument` if `kappa` is not strictly positive.
		{
			if (kappa <= 0)
				throw std::invalid_argument("Error: Ewald parameter must be positive");
			this->kappa = kappa;
			manual = true;
		}
		T ewald_par() const noexcept
		// Return the Ewald parameter.
		{
			return kappa;
		}

		void precise(bool flag) noexcept
		// If `flag` is set to true, the calculations are carried on with the highest
		// precision. By default is set to false, enabling faster calculations.
		{
			fast = !flag;
		}

		void update_ewald() noexcept
		// set update flag to true (calculates optimal Ewald parameter again)
		{
			update = true;
		}

		template <coulomb_and_lj_periodic System>
		void operator()(System& s, utils::thread_pool& tp)
		// Perform summation using the PPPM method with multi-threading.
		// `s` is the physical system.
		// `tp` is a thread pool.
		{
			using std::size_t;
			energy_r = energy_k = energy_scd = energy_lj = 0;
			if (s.n < 2)
			{
				estimated_error = 0;
				energy_coulomb = 0;
				return;
			}

			init_members(s.n, tp.size(), s.side);

			if (update)
				optimize_ewald_par(s);

			// calculate the array part2mesh_inds, which maps the particle indices
			// to their nearest cell (flattened) indices.
			for (size_t i = 0; i < s.n; ++i)
			{
				vec3<size_t> a(nearestcell(s.x[i]));
				part2mesh_inds[i] = (a[0]*num_cells + a[1])*num_cells + a[2];
			}

			init_cell_list(s, tp);

			// calculate real-space contribution to force/energy
			for (size_t i = 0; i < num_threads; ++i)
				tp.enqueue([this, i, &s] { eval_r(s, i); });
			// calculate mesh charges at the same time
			mesh_charges(s);
			tp.wait(); // wait for enqueued tasks to finish
			for (size_t i = 0; i < num_threads; ++i)
			{
				energy_r += partpot_coulomb[i];
				energy_lj += partpot_lj[i];
				s.virial += partvir_lj[i];
			}
			// update influence function if some parameters/conditions have changed
			if (update)
			{
				for (size_t i = 0; i < num_threads; ++i)
					tp.enqueue([this, i] { eval_G(i); });
				tp.wait();
			}
			// calculate reciprocal-space contribution to energy
			forward_transform(tp);
			// calculate reciprocal-space contribution to forces
			backward_transform(tp);
			for (size_t i = 0; i < num_threads; ++i)
				tp.enqueue([this, i, &s] { eval_f_k(s, i); });
			tp.wait();
			// corrections
			eval_ewald_scd(s, energy_scd, kappa, side*side*side, dielec);
			energy_coulomb = energy_r + energy_k + energy_scd;
			s.potential += energy_coulomb + energy_lj;
			s.virial += energy_coulomb;
			update = false;
		}

		T estimated_error, estimated_error_coulomb, estimated_error_lj; // estimated force RMS error
		T energy_coulomb, energy_lj; // energy for electrostatic and LJ potentials
		T energy_r, energy_k, energy_scd; // real, reciprocal and corrections parts of electrostatic potential

		private:

			std::vector<vec3<T>> W, dW, cell_rel_r;
			std::vector<std::complex<T>> zMG, electric_field[3];
			std::vector<T> G, G_energy, partpot_coulomb, partpot_lj, partvir_lj;
			std::vector<unsigned> part2mesh_inds, sorted_inds, first_inds, n_part_cell;
			math::dft<T> dft;
			std::mt19937_64 mersenne_twister = std::mt19937_64(0);
			std::uniform_real_distribution<T> u_dist = std::uniform_real_distribution<T>(0, 1);
			T *zMG_real, *electric_field_real[3];
			T kappa = 0.4, cutoff = 6, dielec = std::numeric_limits<T>::infinity(), cell_size, side;
			std::size_t order = 5, num_cells, m, m_real, num_threads;
			std::ptrdiff_t ch_num_cells = 0, maxn_ = 2; // maxn_ is the max n for the influence function summations
			unsigned update_counter, update_max_count = 100;
			bool update = true, use_ik = true, manual = false, fast = true, verbose = false;

			void init_members(std::size_t n, std::size_t num_threads, T side)
			// initialize member variables
			// `n` is the number of particles (s.n).
			// `num_threads` is the number of threads (tp.size())
			// `side` is the side of the simulation box (s.side)
			{
				using std::size_t;
				int expo = std::ceil(std::log2(n)/3);
				expo = math::clamp(expo, int(2), int(7));
				expo += ch_num_cells;
				expo = math::clamp(expo, int(2), int(7));
				// `num_cells` is the number of cells along one dimension
				num_cells = 1 << expo;
				// `m` is the total number of cells
				m = num_cells*num_cells*num_cells;
				m_real = num_cells*num_cells*(num_cells/2 + 1);

				// update Ewald paremeter every once in a while for constant-P simulations
				if (update_counter % update_max_count == (update_max_count-1))
				{
					update = true;
					++update_counter;
				}

				if (!update)
					if (side != this->side)
					{
						// rescaling Ewald parameter and influence functions
						kappa *= this->side / side;
						T scal = (side * side) / (this->side * this->side);
						for (auto& Gelem : G)
							Gelem *= scal;
						for (auto& Gelem : G_energy)
							Gelem *= scal;

						++update_counter;
					}
				this->side = side;
				this->num_threads = num_threads;
				// `cell_size` is the length of a cell along one dimension
				cell_size = side / num_cells;
				if (cutoff > side/2)
				{
					cutoff = side/2;
					std::clog << "Set cutoff radius to " << cutoff << '\n';
				}

				partpot_coulomb.resize(num_threads);
				partpot_lj.resize(num_threads);
				partvir_lj.resize(num_threads);
				part2mesh_inds.resize(n);
				sorted_inds.resize(n);
				first_inds.resize(m);
				n_part_cell.resize(m);
				zMG.resize(m_real);
				// reinterpreting complex<T> as T[2] is well-defined as per C++ standard
				zMG_real = reinterpret_cast<T*>(zMG.data());
				G.resize(m_real);
				G_energy.resize(m_real);
				W.resize(num_threads*order);
				dW.resize(num_threads*order);
				cell_rel_r.resize(n);
				for (size_t i = 0; i < 3; ++i)
					electric_field[i].resize(m_real);
				for (size_t i = 0; i < 3; ++i)
					electric_field_real[i] = reinterpret_cast<T*>(electric_field[i].data());
			}

			template <typename System>
			void optimize_ewald_par(System& s)
			// Calculate (approximate) optimal Ewald parameter with Newton's method.
			// See Deserno et al., II, 1998
			// `s` is the physical system.
			{
				using std::abs;
				using std::sqrt;
				T cutoff2 = cutoff*cutoff;
				T next_kappa = kappa, errF, corr = 1;
				if (!use_ik)
					// empirical correction (ad-differentiation is 3-4x less accurate)
					corr = 10;
				unsigned counter = 0, max_counter = 100;
				do
				{
					if (!manual)
						// randomize kappa if not converging in 20 iterations
						if ((counter % 20 == 0) || kappa <= 0)
						{
							do
								kappa = u_dist(mersenne_twister);
							while (kappa <= 0);
						}
					kappa = next_kappa;
					T k2 = kappa*kappa;
					// real-space force error
					T errFr = 4/cutoff * math::fastexp(-2*cutoff2*k2);
					T errFr1 = -4*kappa*cutoff2*errFr; // 1st derivative
					T errFr2 = 4*cutoff2*errFr*(4*k2*cutoff2 - 1); // 2nd derivative
					// calculate reciprocal-space force error
					T hk = cell_size*kappa, hk2 = hk*hk;
					T expans = 0, expans1 = 0, expans2 = 0;
					for (ptrdiff_t i = order-1; i >= 0; --i)
						expans = expans * hk2 + error_expansion_coeff<T>[order-1][i];
					for (ptrdiff_t i = order-1; i >= 1; --i)
						expans1 = expans1 * hk2 + i*error_expansion_coeff<T>[order-1][i];
					expans1 *= hk2;
					for (ptrdiff_t i = order-1; i >= 1; --i)
						expans2 = expans2 * hk2 + (i*i)*error_expansion_coeff<T>[order-1][i];
					expans2 *= hk2;
					T hk2P = 1;
					for (unsigned i = 0; i < order; ++i)
						hk2P *= hk2;
					T factor = hk2P * kappa * corr * std::numbers::sqrt2_v<T> / std::numbers::inv_sqrtpi_v<T>;
					// reciprocal-space force error
					T errFk = factor * expans;
					T errFk1 = ((2*order+1)*errFk + 2*factor*expans1)/kappa; // 1st derivative
					T errFk2 = (2*order+1)*(errFk1/kappa - errFk/k2) + (order*expans1 + expans2)*4*factor/k2; // 2nd derivative
					// total force error
					errF = errFr + errFk;
					T errF1 = errFr1 + errFk1; // 1st derivative
					T errF2 = errFr2 + errFk2; // 2nd derivative

					next_kappa -= errF1/errF2;
					++counter;
				}
				while (abs((next_kappa-kappa) / kappa) > 1e-3 && counter < max_counter && !manual);

				T volume = side*side*side;
				estimated_error_coulomb = s.Z2 * sqrt(errF / (s.n * volume));
				estimated_error_lj = 2 * s.tracedisp6 / (3 * cutoff2 * cutoff2) * sqrt(std::numbers::pi_v<T> / (cutoff * s.n * volume));
				estimated_error = sqrt(estimated_error_coulomb*estimated_error_coulomb + estimated_error_lj*estimated_error_lj);

				if (verbose)
				{
					std::clog << "===== PPPM METHOD LOG =====\n";
					std::clog << "Ewald parameter (A^-1): " << kappa << '\n';
					std::clog << "Estimated electrostatic force RMS error (kcal/(mol A)): " << estimated_error_coulomb << '\n';
					std::clog << "estimated total force RMS error (kcal/(mol A)): " << estimated_error << '\n';
					std::clog << "N_M: " << num_cells << '\n';
					std::clog << "P: " << order << '\n';
					std::clog << "differentiation scheme: " << (use_ik ? "ik" : "ad") << '\n';
					std::clog << "h (A): " << cell_size << '\n';
					std::clog << "r_max (A): " << cutoff << "\n\n";
				}
				if (counter == max_counter)
					std::clog << "Warning: Ewald parameter optimization may have not converged!\n";
				if (estimated_error < 1e-5 && fast)
					std::clog << "Warning: estimated accuracy may not be reliable if `precise` is not set to true.\n";
			}

			template <typename Vec>
			vec3<int> nearestcell(const Vec& x)
			// calculate the mesh point (cell) nearest to the position `x`.
			{
				return mod(vec3<int>(round(x*num_cells/side)), int(num_cells));
			}

			vec3<int> index2vec(std::size_t index)
			// calculate the cell vector (i.e. an index for each dimension) from its flattened index.
			// May be used for lattice functions in real space.
			{
				vec3<int> v;
				v[2] = index % num_cells;
				index /= num_cells;
				v[1] = index % num_cells;
				index /= num_cells;
				v[0] = index;
				return v;
			}

			vec3<int> index2vec_complex(std::size_t index)
			// calculate the cell vector (i.e. an index for each dimension) from its flattened index.
			// May be used for lattice functions in reciprocal space, in which case half of the elements
			// is known by symmetry).
			{
				vec3<int> v;
				v[2] = index % (num_cells/2 + 1);
				index /= num_cells/2 + 1;
				v[1] = index % num_cells;
				index /= num_cells;
				v[0] = index;
				return v;
			}

			int remap(int ind)
			// used in `vec2index` and `vec2index_real`.
			{
				return math::mod(ind, int(num_cells));
			}

			std::size_t vec2index(const vec3<int>& v)
			// calculate the cell flattened index from its vector indices.
			// May be used for lattice functions in real space.
			{
				return (remap(v[0])*num_cells + remap(v[1]))*num_cells + remap(v[2]);
			}

			std::size_t vec2index_real(const vec3<int>& v)
			// calculate the cell flattened index from its vector indices.
			// May be used for lattice functions in real space (reinterpreted from complex).
			{
				return (remap(v[0])*num_cells + remap(v[1]))*(num_cells + 2) + remap(v[2]);
			}

			template <typename System>
			void init_cell_list(System& s, utils::thread_pool& tp)
			// initialize the cell list, which will be used in `eval_r`.
			// `s` is the physical system.
			// `tp` is a thread pool.
			{
				using std::size_t;
				auto comp = [this](size_t i, size_t j)
					{
						return part2mesh_inds[i] < part2mesh_inds[j];
					};
				std::iota(sorted_inds.begin(), sorted_inds.end(), 0);
				utils::parallel_sort(sorted_inds.begin(), sorted_inds.end(), tp, comp);
				auto ki = [this](size_t j)
					{
						return part2mesh_inds[sorted_inds[j]];
					};
				// calculate the first sorted index in every cell
				for (size_t k = 0; k <= ki(0); ++k)
					first_inds[k] = 0;
				for (size_t j = 1; j < s.n; ++j)
					for (size_t k = ki(j-1)+1; k <= ki(j); ++k)
						first_inds[k] = j;
				for (size_t k = ki(s.n-1)+1; k < m; ++k)
					first_inds[k] = s.n;
				// calculate number of particles in every cell
				for (size_t k = 0; k < m-1; ++k)
					n_part_cell[k] = first_inds[k+1] - first_inds[k];
				n_part_cell[m-1] = s.n - first_inds[m-1];

				// relative positions of particles w.r.t. the nearest cell
				for (size_t i = 0; i < s.n; ++i)
				{
					vec3<T> cell_r(index2vec(part2mesh_inds[i]));
					cell_r *= cell_size;
					cell_rel_r[i] = remainder(s.x[i] - cell_r, cell_size);
				}
			}

			template <typename System>
			void eval_r(System& s, std::size_t idx)
			// calculate real-space contribution to force/energy.
			// `s` is the physical system.
			// `idx` is the index of the thread.
			// the potential and the virial are accumulated inside `partpot_*[idx]`
			// and `partvir_*[idx]` respectively.
			{
				using std::size_t;
				using std::ceil;
				T uC = 0, uLJ = 0, vLJ = 0;
				T cutoff2 = cutoff*cutoff;
				T cell_cutoff = cutoff + std::numbers::sqrt3_v<T>*cell_size;
				// taking into account particles on the border of a cell
				T cell_cutoff2 = cell_cutoff*cell_cutoff;
				int cell_dist = ceil(cutoff/cell_size);
				size_t block = (s.n-1)/num_threads + 1;
				for (size_t i = idx*block; (i < s.n) && (i < (idx+1)*block); ++i)
				{
					vec3<int> cell(index2vec(part2mesh_inds[i])), kkk;
					// iterating through neighbor cells
					for (kkk[0] = -cell_dist; kkk[0] <= cell_dist; ++kkk[0])
						for (kkk[1] = -cell_dist; kkk[1] <= cell_dist; ++kkk[1])
							for (kkk[2] = -cell_dist; kkk[2] <= cell_dist; ++kkk[2])
							{
								size_t other_k = vec2index(cell + kkk);
								if (n_part_cell[other_k] == 0)
									continue; // there are no particles in this cell
								vec3<T> cell_displ(kkk); // cell displacement vector
								cell_displ *= -cell_size;
								if (dot(cell_displ, cell_displ) > cell_cutoff2)
									continue; // all particles inside the cell are too far away
								size_t first_ind = first_inds[other_k], last_ind = first_ind+n_part_cell[other_k];
								cell_displ += cell_rel_r[i]; // distance between particle i and cell other_k
								// iterating through all particles in the cell other_k
								for (size_t sj = first_ind; sj < last_ind; ++sj)
								{
									size_t j = sorted_inds[sj];
									if (i != j)
									{
										vec3<T> r = cell_displ - cell_rel_r[j];
										T r2 = dot(r, r);
										if (r2 <= cutoff2)
										{
											if (fast)
												eval_ewald_r<true>(s, uC, uLJ, vLJ, kappa, i, j, r, r2);
											else
												eval_ewald_r<false>(s, uC, uLJ, vLJ, kappa, i, j, r, r2);
										}
									}
								}
							}
				}
				partpot_coulomb[idx] = uC/2;
				partpot_lj[idx] = uLJ/2;
				partvir_lj[idx] = vLJ/2;
			}

			template <typename Vec>
			vec3<T> charge_assignment_center(const Vec& x, std::size_t p)
			// Calculate the charge assignment center.
			// `x` is the particle position.
			// `p` is the nearest cell flattened index (which may be given by part2mesh_inds[i],
			// where `i` is the particle index).
			// The charge assignment center corresponds to the nearest cell if order is odd,
			// otherwise it corresponds to the midpoint between the two nearest cells along each
			// dimension (which, in other words, may be viewed as the nearest "vertex" if cells
			// are seen as small cubes).
			// Called by `mesh_charges` and `eval_f_k`.
			{
				using std::remainder;
				// calculate nearest cell vector indices
				vec3<T> nearest(index2vec(p));
				if (math::is_even(order))
				{
					// calculate nearest vertex
					for (std::size_t j = 0; j < 3; ++j)
						if (remainder(x[j] - nearest[j] * cell_size, side) < 0)
							nearest[j] -= 1;
					return nearest + T(0.5);
				}
				else
					return nearest;
			}

			template <typename System>
			void mesh_charges(System& s)
			// calculate mesh charges, interpolating the particle charges onto a lattice.
			// `s` is the physical system.
			// The mesh charges will be stored in the `zMG` container (interpreted as real).
			{
				using std::size_t;
				std::ranges::fill(zMG, 0);
				for (size_t i = 0; i < s.n; ++i)
				{
					vec3<T> center(charge_assignment_center(s.x[i], part2mesh_inds[i]));
					vec3<T> x_prime(0);
					if (order > 1)
						x_prime = remainder(s.x[i]/cell_size - center, T(num_cells));
					for (size_t k = 0; k < order; ++k)
						for (size_t j = 0; j < 3; ++j)
							W[k][j] = charge_assignment_function(x_prime[j], k, order);
					vec3<int> kkk, cell(mod(center, T(num_cells)));
					cell -= int(order-1)/2;
					for (kkk[0] = 0; kkk[0] < int(order); ++kkk[0])
						for (kkk[1] = 0; kkk[1] < int(order); ++kkk[1])
							for (kkk[2] = 0; kkk[2] < int(order); ++kkk[2])
							{
								size_t index = vec2index_real(cell + kkk);
								zMG_real[index] += s.z[i] * W[kkk[0]][0] * W[kkk[1]][1] * W[kkk[2]][2];
							}
				}
			}

			void influence_function_ik(const vec3<T>& kvec, T den, std::size_t i)
			// calculate influence function, also known as the lattice Green function,
			// using ik-differentiation, at index `i`.
			// `kvec` is the wavevector associated to index `i`.
			// `den` is the influence function denominator (not squared).
			// The optimal influence function for forces will be `G` and the optimal
			// influence function for energy will be `G_energy`.
			{
				const T m1_4kappa2 = -1/(4*kappa*kappa);
				const T pi_h = std::numbers::pi_v<T>/cell_size, twopi_h = 2*pi_h;
				vec3<int> nnn_;
				T k2 = dot(kvec, kvec), num = 0, num_energy = 0;
				T s0, s1, s2;
				for (nnn_[0] = -maxn_; nnn_[0] <= maxn_; ++nnn_[0])
				{
					vec3<T> k_;
					k_[0] = kvec[0] + twopi_h * nnn_[0];
					s0 = math::sinc(k_[0]*cell_size/2);
					for (nnn_[1] = -maxn_; nnn_[1] <= maxn_; ++nnn_[1])
					{
						k_[1] = kvec[1] + twopi_h * nnn_[1];
						s1 = s0*math::sinc(k_[1]*cell_size/2);
						for (nnn_[2] = -maxn_; nnn_[2] <= maxn_; ++nnn_[2])
						{
							k_[2] = kvec[2] + twopi_h * nnn_[2];
							s2 = s1*math::sinc(k_[2]*cell_size/2);
							T k_2 = dot(k_, k_);
							T U2 = 1;
							for (unsigned k = 0; k < order; ++k)
								U2 *= s2;
							U2 *= U2;
							T factor = 0;
							if (fast)
								factor = U2 / k_2 * math::fastexp(k_2*m1_4kappa2);
							else
								factor = U2 / k_2 * std::exp(k_2*m1_4kappa2);
							num += dot(kvec, k_) * factor;
							num_energy += factor;
						}
					}
				}
				G[i] = num * 2 * math::two_pi<T> / (k2*den*den);
				G_energy[i] = num_energy * 2 * math::two_pi<T> / (den*den);
			}

			void influence_function_ad(const vec3<T>& kvec, T den, std::size_t i)
			// calculate influence function, also known as the lattice Green function,
			// using analytic differentiation in real space, at index `i`.
			// `kvec` is the wavevector associated to index `i`.
			// `den` is the influence function part of the denominator which can be computed with closed formulas.
			// The optimal influence function for forces will be `G` and the optimal
			// influence function for energy will be `G_energy`.
			{
				const T m1_4kappa2 = -1/(4*kappa*kappa);
				const T pi_h = std::numbers::pi_v<T>/cell_size, twopi_h = 2*pi_h;
				vec3<int> nnn_;
				T num = 0, num_energy = 0, den_ = 0;
				T s0, s1, s2;
				for (nnn_[0] = -maxn_; nnn_[0] <= maxn_; ++nnn_[0])
				{
					vec3<T> k_;
					k_[0] = kvec[0] + twopi_h * nnn_[0];
					s0 = math::sinc(k_[0]*cell_size/2);
					for (nnn_[1] = -maxn_; nnn_[1] <= maxn_; ++nnn_[1])
					{
						k_[1] = kvec[1] + twopi_h * nnn_[1];
						s1 = s0*math::sinc(k_[1]*cell_size/2);
						for (nnn_[2] = -maxn_; nnn_[2] <= maxn_; ++nnn_[2])
						{
							k_[2] = kvec[2] + twopi_h * nnn_[2];
							s2 = s1*math::sinc(k_[2]*cell_size/2);
							T k_2 = dot(k_, k_);
							T U2 = 1;
							for (unsigned k = 0; k < order; ++k)
								U2 *= s2;
							U2 *= U2;
							T factor = 0;
							if (fast)
								factor = U2 * math::fastexp(k_2*m1_4kappa2);
							else
								factor = U2 * std::exp(k_2*m1_4kappa2);
							num += factor;
							num_energy += factor / k_2;
							den_ += U2 * k_2;
						}
					}
				}
				G[i] = num * 2 * math::two_pi<T> / (den*den_);
				G_energy[i] = num_energy * 2 * math::two_pi<T> / (den*den);
			}
			
			void eval_G(std::size_t idx)
			// calculate influence function, also known as the lattice Green function.
			// `idx` is the index of the thread.
			{
				using std::size_t;
				using std::sin;
				const T twopi_L = 2*std::numbers::pi_v<T>/side;
				size_t block = (m_real-1)/num_threads + 1;
				for (size_t i = idx*block; (i < m_real) && (i < (idx+1)*block); ++i)
				{
					if (i == 0)
					{
						G[i] = G_energy[i] = 0;
						continue;
					}
					auto nvec = index2vec_complex(i);
					nvec -= int(num_cells)*(2*nvec/int(num_cells));
					vec3<T> kvec = twopi_L * vec3<T>(nvec), s = kvec*cell_size/2;
					for (size_t j = 0; j < 3; ++j)
						s[j] = sin(s[j]);
					T den = 1;
					for (size_t j = 0; j < 3; ++j)
						den *= influence_function_den<T>[order-1](s[j]*s[j]);
					if (use_ik)
						influence_function_ik(kvec, den, i);
					else
						influence_function_ad(kvec, den, i);
				}
			}

			void forward_transform(utils::thread_pool& tp)
			// calculate DFT of the mesh charges and accumulate
			// the reciprocal-space contribution to energy.
			// The DFT of the meshes is put inside `zMG` (in-place).
			// `tp` is a thread pool.
			{
				using std::size_t;
				std::array<size_t, 3> nlist{num_cells, num_cells, num_cells/2 + 1};
				dft.template rfftn<3>(zMG, nlist, tp);

				T volume = side*side*side;
				T u = 0;
				for (size_t i = 0; i < m_real; ++i)
					u += norm(zMG[i]) * G_energy[i];
				energy_k = u / volume;
			}

			void backward_transform(utils::thread_pool& tp)
			// Calculate the electric potential in reciprocal space, and then:
			// * if ik-differentiation is chosen, calculate the gradient in reciprocal
			//   space and perform three inverse DFT to obtain the three components of the
			//   electric field on the mesh points (contained inside `electric_field`).
			// * if ad-differentiation is chosen, perform an inverse DFT to obtain the
			//   electric potential on the mesh points (contained inside `zMG`).
			// `tp` is a thread pool.
			{
				using std::size_t;
				std::array<size_t, 3> nlist{num_cells, num_cells, num_cells/2 + 1};
				T cell_size3 = cell_size*cell_size*cell_size;
				for (size_t i = 0; i < m_real; ++i)
					zMG[i] *= G[i] / cell_size3;
				if (use_ik)
				{
					T twopi_L = 2*std::numbers::pi_v<T>/side;
					vec3<size_t> nnn;
					vec3<T> kvec;
					for (nnn[0] = 0; nnn[0] < num_cells; ++nnn[0])
						for (nnn[1] = 0; nnn[1] < num_cells; ++nnn[1])
							for (nnn[2] = 0; nnn[2] < num_cells/2 + 1; ++nnn[2])
							{
								size_t index = (nnn[0]*num_cells + nnn[1])*(num_cells/2 + 1) + nnn[2];
								vec3<int> nvec(nnn);
								nvec -= int(num_cells)*(2*nvec/int(num_cells));
								kvec = twopi_L * vec3<T>(nvec);
								std::complex<T> zGi = zMG[index] * std::complex<T>(0, 1);
								for (size_t j = 0; j < 3; ++j)
									electric_field[j][index] = kvec[j] * zGi;
							}
					for (size_t j = 0; j < 3; ++j)
						dft.template irfftn<3>(electric_field[j], nlist, tp);
				}
				else
					dft.template irfftn<3>(zMG, nlist, tp);
			}

			template <typename System>
			void eval_f_k(System& s, std::size_t idx)
			// compute reciprocal-space contribution to forces.
			// `s` is the physical system.
			// `idx` is the index of the thread.
			{
				using std::size_t;
				size_t block = (s.n-1)/num_threads + 1;
				for (size_t i = idx*block; (i < s.n) && (i < (idx+1)*block); ++i)
				{
					vec3<T> center(charge_assignment_center(s.x[i], part2mesh_inds[i]));
					vec3<T> x_prime(0);
					if (order > 1)
						x_prime = remainder(s.x[i]/cell_size - center, T(num_cells));
					for (size_t k = 0; k < order; ++k)
						for (size_t j = 0; j < 3; ++j)
							W[order*idx + k][j] = charge_assignment_function(x_prime[j], k, order);
					// compute charge assignment function derivatives for ad-differentiation
					if (!use_ik)
						for (size_t k = 0; k < order; ++k)
							for (size_t j = 0; j < 3; ++j)
								dW[order*idx + k][j] = charge_assignment_function_der(x_prime[j], k, order)/cell_size;
					vec3<T> ei(0);
					vec3<int> kkk, cell(mod(center, T(num_cells)));
					cell -= int(order-1)/2;
					for (kkk[0] = 0; kkk[0] < int(order); ++kkk[0])
						for (kkk[1] = 0; kkk[1] < int(order); ++kkk[1])
							for (kkk[2] = 0; kkk[2] < int(order); ++kkk[2])
							{
								size_t index = vec2index_real(cell + kkk);
								if (use_ik)
								{
									T Wxyz = W[order*idx + kkk[0]][0] * W[order*idx + kkk[1]][1] * W[order*idx + kkk[2]][2];
									for (size_t j = 0; j < 3; ++j)
										ei[j] -= electric_field_real[j][index] * Wxyz;
								}
								else
								{
									vec3<T> gradW(
										dW[order*idx + kkk[0]][0] * W[order*idx + kkk[1]][1] * W[order*idx + kkk[2]][2],
										W[order*idx + kkk[0]][0] * dW[order*idx + kkk[1]][1] * W[order*idx + kkk[2]][2],
										W[order*idx + kkk[0]][0] * W[order*idx + kkk[1]][1] * dW[order*idx + kkk[2]][2]
									);
									ei -= zMG_real[index] * gradW;
								}
							}
					s.f[i] += s.z[i] * ei;
				}
			}
	};

} // namespace physics

#endif // PHYSICS_PPPM_H
































