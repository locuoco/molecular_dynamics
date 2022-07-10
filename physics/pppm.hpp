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
#include <cmath> // round, ceil, log2, sqrt, sin, cos, fabs, remainder
#include <numeric> // iota
#include <numbers>
#include <functional> // function
#include <algorithm> // fill, copy
#include <complex>
#include <stdexcept> // invalid_argument
#include <random> // mt19937, uniform_real_distribution

#include "ewald.hpp" // from math/helper.hpp: fasterfc, fastexp, mod
#include "../math/dft.hpp"

namespace physics
{
	static constexpr std::size_t pppm_min_order = 1, pppm_max_order = 7;

	template <typename T>
	constexpr T charge_assignment_function_coeff[pppm_max_order][pppm_max_order][pppm_max_order] =
		{
			// P = 1
			{{1}},
			// P = 2
			{
				{1.L/2, -1},
				{1.L/2, +1},
			},
			// P = 3
			{
				{1.L/8, -1.L/2, 1.L/2},
				{3.L/4, 0, -1},
				{1.L/8, +1.L/2, 1.L/2},
			},
			// P = 4
			{
				{1.L/48, -1.L/8, 1.L/4, -1.L/6},
				{23.L/48, -5.L/8, -1.L/4, +1.L/2},
				{23.L/48, +5.L/8, -1.L/4, -1.L/2},
				{1.L/48, +1.L/8, 1.L/4, +1.L/6},
			},
			// P = 5
			{
				{1.L/384, -1.L/48, 1.L/16, -1.L/12, 1.L/24},
				{19.L/96, -11.L/24, 1.L/4, +1.L/6, -1.L/6},
				{115.L/192, 0, -5.L/8, 0, 1.L/4},
				{19.L/96, +11.L/24, 1.L/4, -1.L/6, -1.L/6},
				{1.L/384, +1.L/48, 1.L/16, +1.L/12, 1.L/24},
			},
			// P = 6
			{
				{1.L/3840, -1.L/384, 1.L/96, -1.L/48, 1.L/48, -1.L/120},
				{237.L/3840, -25.L/128, 21.L/96, -1.L/16, -1.L/16, +1.L/24},
				{841.L/1920, -77.L/192, -11.L/48, +7.L/24, 1.L/24, -1.L/12},
				{841.L/1920, +77.L/192, -11.L/48, -7.L/24, 1.L/24, +1.L/12},
				{237.L/3840, +25.L/128, 21.L/96, +1.L/16, -1.L/16, -1.L/24},
				{1.L/3840, +1.L/384, 1.L/96, +1.L/48, 1.L/48, +1.L/120},
			},
			// P = 7
			{
				{1.L/46080, -1.L/3840, 1.L/768, -1.L/288, 1.L/192, -1.L/240, 1.L/720},
				{361.L/23040, -59.L/960, 37.L/384, -5.L/72, 1.L/96, +1.L/60, -1.L/120},
				{10543.L/46080, -289.L/768, 79.L/768, +43.L/288, -17.L/192, -1.L/48, 1.L/48},
				{5887.L/11520, 0, -77.L/192, 0, 21.L/144, 0, -1.L/36},
				{10543.L/46080, +289.L/768, 79.L/768, -43.L/288, -17.L/192, +1.L/48, 1.L/48},
				{361.L/23040, +59.L/960, 37.L/384, +5.L/72, 1.L/96, -1.L/60, -1.L/120},
				{1.L/46080, +1.L/3840, 1.L/768, +1.L/288, 1.L/192, +1.L/240, 1.L/720},
			},
		};

	template <typename T>
	constexpr T charge_assignment_function(T x, std::size_t pos, std::size_t order) noexcept
	{
		T res = 0;
		for (std::ptrdiff_t i = order-1; i >= 0; --i)
			res = res * x + charge_assignment_function_coeff<T>[order-1][pos][i];
		return res;
	}
	template <typename T>
	constexpr T charge_assignment_function_der(T x, std::size_t pos, std::size_t order) noexcept
	{
		T res = 0;
		for (std::ptrdiff_t i = order-1; i >= 1; --i)
			res = res * x + i*charge_assignment_function_coeff<T>[order-1][pos][i];
		return res;
	}

	template <typename T>
	const std::function<T(T)> influence_function_den[pppm_max_order] =
		{
			[](T) { return T(1); }, // P = 1
			[](T s2) { return 1 - T(2.L/3)*s2; }, // P = 2
			[](T s2) { return 1 - s2 + T(2.L/15)*s2*s2; }, // P = 3
			[](T s2) { T s4(s2*s2); return 1 - T(4.L/3)*s2 + T(2.L/5)*s4 - T(4.L/315)*s4*s2; }, // P = 4
			[](T s2) { T s4(s2*s2); return 1 - T(5.L/3)*s2 + T(7.L/9)*s4 - T(17.L/189)*s4*s2 + T(2.L/2835)*s4*s4; }, // P = 5
			[](T s2) { T s4(s2*s2), s8(s4*s4); return 1 - 2*s2 + T(19.L/15)*s4 - T(256.L/945)*s4*s2 + T(62.L/4725)*s8 - T(4.L/155'925)*s8*s2; }, // P = 6
			[](T s2)
			{
				T s4(s2*s2), s8(s4*s4);
				return 1 - T(7.L/3)*s2 + T(28.L/15)*s4 - T(16.L/27)*s4*s2 + T(26.L/405)*s8 - T(2.L/1485)*s8*s2 + T(4.L/6'081'075)*s8*s4;
			}, // P = 7
		};

	template <typename T>
	constexpr T error_expansion_coeff[pppm_max_order][pppm_max_order] =
		{
			{2.L/3}, // P = 1
			{1.L/50, 5.L/294}, // P = 2
			{1.L/588, 7.L/1440, 21.L/3872}, // P = 3
			{1.L/4320, 3.L/1936, 7601.L/2'271'360, 143.L/28800}, // P = 4
			{1.L/23232, 7601.L/13'628'160, 143.L/69120, 517'231.L/106'536'960, 106'640'677.L/11'737'571'328}, // P = 5
			{691.L/68'140'800, 13.L/57600, 47021.L/35'512'320, 9'694'607.L/2'095'994'880, 733'191'589.L/59'609'088'000, 326'190'917.L/11'700'633'600}, // P = 6
			{
				1.L/345'600, 3617.L/35'512'320, 745'739.L/838'397'952, 56'399'353.L/12'773'376'000, 25'091'609.L/1'560'084'480,
				1'755'948'832'039.L/36'229'939'200'000, 4'887'769'399.L/37'838'389'248
			}, // P = 7
		};

	template <typename T, typename State>
	struct pppm
	// Particle-particle, particle-mesh method (also known as PPPM or P^3M or P3M method)
	// based on Ewald summation, it has O(N log N) asymptotic complexity
	{
		void charge_assignment_order(std::size_t order)
		{
			if (this -> order != order)
			{
				if (order < pppm_min_order)
					throw std::invalid_argument("Error: minimum PPPM order is 1");
				if (order > pppm_max_order)
					throw std::invalid_argument("Error: maximum PPPM order is 7");
				this -> order = order;
				update = true;
			}
		}

		void cutoff_radius(T cutoff)
		{
			if (this -> cutoff != cutoff)
			{
				if (cutoff <= 0)
					throw std::invalid_argument("Error: cutoff radius must be positive");
				this -> cutoff = cutoff;
				update = true;
			}
		}

		void choose_diff_scheme(std::string scheme)
		{
			if (scheme == std::string("ik"))
				use_ik = true;
			else if (scheme == std::string("ad"))
				use_ik = false;
			else
				throw std::invalid_argument(std::string("Error: Unknown differentiation scheme: ") + scheme);
		}

		T ewald_par() const noexcept
		{
			return kappa;
		}

		template <coulomb_and_LJ_periodic<T, State> S>
		void operator()(S& s, utils::thread_pool& tp)
		{
			using std::size_t;
			using std::ptrdiff_t;
			using std::ceil;
			using std::sin;
			using std::remainder;
			using std::fabs;
			if (s.n < 2)
				return;

			/*size_t num_cells = 2;
			T cell_size = s.side / num_cells;
			while (cell_size > 1)
			{
				num_cells *= 2;
				cell_size = s.side / num_cells;
			}*/
			size_t num_cells = 1ull << size_t(std::round(std::log2(s.n)/3));
			size_t m = num_cells*num_cells*num_cells;
			size_t num_threads = tp.size();
			T cell_size = s.side / num_cells, cell_size3 = cell_size*cell_size*cell_size;
			T volume = s.side * s.side * s.side;
			if (cutoff > s.side/2)
			{
				cutoff = s.side/2;
				std::clog << "Set cutoff to " << s.side/2 << '\n';
			}

			partU.resize(num_threads);
			partvir.resize(num_threads);
			part2mesh_inds.resize(s.n);
			sorted_inds.resize(s.n);
			first_inds.resize(m);
			n_part_cell.resize(m);
			zMG.resize(m/2);
			// reinterpreting complex<T> as T[2] is well-defined per C++ standard
			zMG_real = reinterpret_cast<T*>(zMG.data());
			G.resize(m/2);
			W.resize(num_threads*order);
			dW.resize(num_threads*order);
			for (size_t i = 0; i < 3; ++i)
				electric_field[i].resize(m/2);
			for (size_t i = 0; i < 3; ++i)
				electric_field_real[i] = reinterpret_cast<T*>(electric_field[i].data());

			auto nearestmesh = [num_cells, &s](size_t i)
				{
					//return trunc((remainder(s.x[i] + T(0.5)/num_cells, s.side) / s.side + T(0.5))*num_cells/(1+eps));
					return mod(round(s.x[i]/s.side*num_cells), num_cells);
				};

			auto index2vec = [num_cells](size_t index)
				{
					point3<ptrdiff_t> vec;
					vec[2] = index % num_cells;
					index /= num_cells;
					vec[1] = index % num_cells;
					index /= num_cells;
					vec[0] = index;
					return vec;
				};

			auto index2vec_complex = [num_cells](size_t index)
				{
					point3<ptrdiff_t> vec;
					vec[2] = index % (num_cells/2);
					index /= num_cells/2;
					vec[1] = index % num_cells;
					index /= num_cells;
					vec[0] = index;
					return vec;
				};

			auto vec2index = [num_cells](point3<ptrdiff_t> vec)
				{
					return (math::mod(vec[0], ptrdiff_t(num_cells))*ptrdiff_t(num_cells) +
							math::mod(vec[1], ptrdiff_t(num_cells)))*ptrdiff_t(num_cells) +
						math::mod(vec[2], ptrdiff_t(num_cells));
				};

			for (size_t i = 0; i < s.n; ++i)
			{
				point3<size_t> a(nearestmesh(i));
				part2mesh_inds[i] = (a[0]*num_cells + a[1])*num_cells + a[2];
			}
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
			for (size_t k = 0; k <= ki(0); ++k)
				first_inds[k] = 0;
			for (size_t j = 1; j < s.n; ++j)
				for (size_t k = ki(j-1)+1; k <= ki(j); ++k)
					first_inds[k] = j;
			for (size_t k = ki(s.n-1)+1; k < m; ++k)
				first_inds[k] = s.n;
			for (size_t k = 0; k < m-1; ++k)
				n_part_cell[k] = first_inds[k+1] - first_inds[k];
			n_part_cell[m-1] = s.n - first_inds[m-1];

			if (update)
			{
				// calculate (approximate) optimal Ewald parameter
				T cutoff2 = cutoff*cutoff;
				T prev_kappa, errF;
				unsigned counter = 0, max_counter = 1000;
				do
				{
					if ((counter % 10 == 0 && counter < max_counter-25) || kappa <= 0)
					{
						do
							kappa = u_dist(mersenne_twister);
						while (kappa <= 0);
					}
					prev_kappa = kappa;
					T k2 = kappa*kappa;
					T errFr = 2/sqrt(cutoff) * math::fastexp(-cutoff2*k2);
					T errFr1 = -2*kappa*cutoff2*errFr; // 1st derivative
					T errFr2 = 2*cutoff2*errFr*(2*k2*cutoff2 - 1); // 2nd derivative
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
					T hkP = 1;
					for (unsigned i = 0; i < order; ++i)
						hkP *= hk;
					T errFk = hkP * sqrt(kappa * std::numbers::sqrt2_v<T> / std::numbers::inv_sqrtpi_v<T> * expans);
					T factor = 2*order + 1 + 2 * expans1/expans;
					T errFk1 = errFk*factor/(2*kappa); // 1st derivative
					T errFk2 = (errFk1/(2*kappa) - errFk/(2*k2))*factor + 2*errFk/k2*(expans*expans2 - expans1*expans1)/(expans*expans); // 2nd derivative
					errF = errFr + errFk;
					T errF1 = errFr1 + errFk1;
					T errF2 = errFr2 + errFk2;
					kappa -= errF1/errF2;
					//std::cout << "k = " << kappa << '\n'; // debug
					++counter;
				}
				while (fabs(kappa-prev_kappa) / fabs(kappa) > 1e-3 && counter < max_counter);
				if (counter == max_counter)
					std::clog << "Warning: Ewald parameter optimization not converged!\n";
				std::clog << "ewald par (A): " << kappa << '\n';
				std::clog << "estimated force error (kcal/(mol A)): " << errF << '\n';
				std::clog << "N_M: " << num_cells << '\n';
				std::clog << "h (A): " << cell_size << '\n';
			}

			// calculate real-space contribution to force/energy
			auto eval_r = [this, num_threads, cell_size, num_cells, ki, index2vec, vec2index, &s](size_t idx)
				{
					T uC = 0, u = 0, v = 0;
					T cutoff2 = cutoff*cutoff;
					size_t block = (s.n-1)/num_threads + 1;
					for (size_t i = idx*block; (i < s.n) && (i < (idx+1)*block); ++i)
					{
						size_t si = sorted_inds[i];
						ptrdiff_t cell_dist = ceil(cutoff/cell_size);
						point3<ptrdiff_t> cell(index2vec(part2mesh_inds[si])), kkk;
						// iterating through neighbor cells
						for (kkk[0] = -cell_dist; kkk[0] <= cell_dist; ++kkk[0])
							for (kkk[1] = -cell_dist; kkk[1] <= cell_dist; ++kkk[1])
								for (kkk[2] = -cell_dist; kkk[2] <= cell_dist; ++kkk[2])
								{
									T min_r2 = (dot(kkk, kkk) - 3) * cell_size*cell_size;
									if (min_r2 > cutoff2)
										continue; // all particles inside the cell are too far away
									size_t other_k = vec2index(cell + kkk);
									size_t first_ind = first_inds[other_k], last_ind = first_ind+n_part_cell[other_k];
									// iterating through all particles in the cell
									for (size_t j = first_ind; j < last_ind; ++j)
										if (i != j)
										{
											size_t sj = sorted_inds[j];
											auto r = remainder(s.x[si] - s.x[sj], s.side);
											T r2 = dot(r, r);
											if (r2 <= cutoff2)
												eval_ewald_r(s, u, uC, v, kappa, si, sj, r, r2);
										}
								}
					}
					partU[idx] = (u + uC)/2;
					partvir[idx] = (v + uC)/2;
				};
			for (size_t i = 0; i < num_threads; ++i)
				tp.enqueue(eval_r, i);
			// calculate mesh charges
			auto charge_assignment_center = [this, cell_size, num_cells, index2vec, &s](size_t i, size_t p)
				{
					// calculate nearest cell
					point3<ptrdiff_t> nearest_cell(index2vec(p));
					if ((order & 1) == 0) // order is even
					{
						point3<ptrdiff_t> nearest_vertex = nearest_cell;
						// calculate nearest vertex
						for (size_t j = 0; j < 3; ++j)
							if (remainder(nearest_cell[j] * cell_size - s.x[i][j], s.side) > 0)
								--nearest_vertex[j];
						return point3<T>(nearest_vertex) + T(0.5);
					}
					else
						return point3<T>(nearest_cell);
				};
			std::fill(zMG.begin(), zMG.end(), 0);
			for (size_t i = 0; i < s.n; ++i)
			{
				ptrdiff_t p = part2mesh_inds[i];
				point3<T> center(charge_assignment_center(i, p));
				point3<T> x_prime(0);
				if (order > 1)
					x_prime = remainder(s.x[i]/cell_size - center, T(num_cells));
				for (size_t k = 0; k < order; ++k)
					for (size_t j = 0; j < 3; ++j)
						W[k][j] = charge_assignment_function(x_prime[j], k, order);
				point3<ptrdiff_t> kkk, cell(mod(center, T(num_cells)));
				ptrdiff_t half_order = (order-1)/2;
				cell -= half_order;
				for (kkk[0] = 0; kkk[0] < ptrdiff_t(order); ++kkk[0])
					for (kkk[1] = 0; kkk[1] < ptrdiff_t(order); ++kkk[1])
						for (kkk[2] = 0; kkk[2] < ptrdiff_t(order); ++kkk[2])
						{
							size_t index = vec2index(cell + kkk);
							zMG_real[index] += s.z[i] * W[kkk[0]][0] * W[kkk[1]][1] * W[kkk[2]][2];
						}
			}
			tp.wait();
			for (size_t i = 0; i < num_threads; ++i)
			{
				s.U += partU[i];
				s.virial += partvir[i];
			}
			// calculate influence function
			T m1_4kappa2 = -1/(4*kappa*kappa);
			T pi_h = std::numbers::pi_v<T>/cell_size, twopi_h = 2*pi_h, twopi_L = twopi_h/num_cells;
			ptrdiff_t maxn_ = 2;
			auto influence_function_ik = [this, cell_size, m1_4kappa2, twopi_h, maxn_](point3<T> kvec, T den, size_t i)
				{
					point3<ptrdiff_t> nnn_;
					T k2 = dot(kvec, kvec), num = 0;
					for (nnn_[0] = -maxn_; nnn_[0] <= maxn_; ++nnn_[0])
					{
						point3<T> k_, s_;
						k_[0] = kvec[0] + twopi_h * nnn_[0];
						s_[0] = math::sinc(k_[0]*cell_size/2);
						for (nnn_[1] = -maxn_; nnn_[1] <= maxn_; ++nnn_[1])
						{
							k_[1] = kvec[1] + twopi_h * nnn_[1];
							s_[1] = s_[0]*math::sinc(k_[1]*cell_size/2);
							for (nnn_[2] = -maxn_; nnn_[2] <= maxn_; ++nnn_[2])
							{
								k_[2] = kvec[2] + twopi_h * nnn_[2];
								s_[2] = s_[1]*math::sinc(k_[2]*cell_size/2);
								T k_2 = dot(k_, k_);
								T U2 = 1;
								for (unsigned k = 0; k < order; ++k)
									U2 *= s_[2];
								U2 *= U2;
								num += dot(kvec, k_) * U2 / k_2 * math::fastexp(k_2*m1_4kappa2);
							}
						}
					}
					G[i] = num * 2 * math::two_pi<T>() / (k2*den*den);
				};
			auto influence_function_ad = [this, cell_size, m1_4kappa2, twopi_h, maxn_](point3<T> kvec, T den, size_t i)
				{
					point3<ptrdiff_t> nnn_;
					T num = 0, den_ = 0;
					for (nnn_[0] = -maxn_; nnn_[0] <= maxn_; ++nnn_[0])
					{
						point3<T> k_, s_;
						k_[0] = kvec[0] + twopi_h * nnn_[0];
						s_[0] = math::sinc(k_[0]*cell_size/2);
						for (nnn_[1] = -maxn_; nnn_[1] <= maxn_; ++nnn_[1])
						{
							k_[1] = kvec[1] + twopi_h * nnn_[1];
							s_[1] = s_[0]*math::sinc(k_[1]*cell_size/2);
							for (nnn_[2] = -maxn_; nnn_[2] <= maxn_; ++nnn_[2])
							{
								k_[2] = kvec[2] + twopi_h * nnn_[2];
								s_[2] = s_[1]*math::sinc(k_[2]*cell_size/2);
								T k_2 = dot(k_, k_);
								T U2 = 1;
								for (unsigned k = 0; k < order; ++k)
									U2 *= s_[2];
								U2 *= U2;
								num += U2 * math::fastexp(k_2*m1_4kappa2);
								den_ += U2 * k_2;
							}
						}
					}
					G[i] = num * 2 * math::two_pi<T>() / (den*den_);
				};
			
			auto eval_G = [this, num_threads, num_cells, cell_size, m, index2vec_complex, twopi_L, influence_function_ik, influence_function_ad](size_t idx)
				{
					size_t block = (m/2-1)/num_threads + 1;
					for (size_t i = idx*block; (i < m/2) && (i < (idx+1)*block); ++i)
					{
						if (i == 0)
						{
							G[i] = 0;
							continue;
						}
						auto nvec = index2vec_complex(i);
						nvec -= ptrdiff_t(num_cells)*(2*nvec/ptrdiff_t(num_cells));
						point3<T> kvec = twopi_L * point3<T>(nvec), s = kvec*cell_size/2;
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
				};
			if (update)
			{
				for (size_t i = 0; i < num_threads; ++i)
					tp.enqueue(eval_G, i);
				tp.wait();
				update = false;
			}
			// calculate reciprocal-space contribution to energy
			std::array<size_t, 3> nlist{num_cells, num_cells, num_cells/2};
			dft.template rfftn<3>(zMG, nlist, tp);
			{
				T u{};
				for (size_t i = 0; i < m/2; ++i)
					u += norm(zMG[i]) * G[i];
				u /= volume;
				s.U += u;
				s.virial += u;
			}
			// calculate reciprocal-space contribution to forces
			for (size_t i = 0; i < m/2; ++i)
				zMG[i] *= G[i] / cell_size3;
			if (use_ik)
			{
				point3<size_t> nnn;
				point3<T> kvec;
				for (nnn[0] = 0; nnn[0] < num_cells; ++nnn[0])
					for (nnn[1] = 0; nnn[1] < num_cells; ++nnn[1])
						for (nnn[2] = 0; nnn[2] < num_cells/2; ++nnn[2])
						{
							size_t index = (nnn[0]*num_cells + nnn[1])*(num_cells/2) + nnn[2];
							point3<ptrdiff_t> nvec(nnn);
							for (size_t j = 0; j < 2; ++j)
								if (nvec[j] == ptrdiff_t(num_cells/2))
									nvec[j] = 0;
							nvec -= ptrdiff_t(num_cells)*(2*nvec/ptrdiff_t(num_cells));
							kvec = twopi_L * point3<T>(nvec);
							std::complex<T> zGi = zMG[index] * std::complex<T>(0, 1);
							for (size_t j = 0; j < 3; ++j)
								electric_field[j][index] = kvec[j] * zGi;
						}
				for (size_t j = 0; j < 3; ++j)
					dft.template irfftn<3>(electric_field[j], nlist, tp);
			}
			else
				dft.template irfftn<3>(zMG, nlist, tp);
			auto eval_f_k = [this, num_threads, cell_size, num_cells, vec2index, charge_assignment_center, &s] (size_t idx)
				{
					size_t block = (s.n-1)/num_threads + 1;
					for (size_t i = idx*block; (i < s.n) && (i < (idx+1)*block); ++i)
					{
						ptrdiff_t p = part2mesh_inds[i];
						point3<T> center(charge_assignment_center(i, p));
						point3<T> x_prime(0);
						if (order > 1)
							x_prime = remainder(s.x[i]/cell_size - center, T(num_cells));
						for (size_t k = 0; k < order; ++k)
							for (size_t j = 0; j < 3; ++j)
								W[order*idx + k][j] = charge_assignment_function(x_prime[j], k, order);
						if (!use_ik)
							for (size_t k = 0; k < order; ++k)
								for (size_t j = 0; j < 3; ++j)
									dW[order*idx + k][j] = charge_assignment_function_der(x_prime[j], k, order);
						point3<T> ei(0);
						point3<ptrdiff_t> kkk, cell(mod(center, T(num_cells)));
						ptrdiff_t half_order = (order-1)/2;
						cell -= half_order;
						for (kkk[0] = 0; kkk[0] < ptrdiff_t(order); ++kkk[0])
							for (kkk[1] = 0; kkk[1] < ptrdiff_t(order); ++kkk[1])
								for (kkk[2] = 0; kkk[2] < ptrdiff_t(order); ++kkk[2])
								{
									size_t index = vec2index(cell + kkk);
									if (use_ik)
									{
										T Wxyz = W[order*idx + kkk[0]][0] * W[order*idx + kkk[1]][1] * W[order*idx + kkk[2]][2];
										for (size_t j = 0; j < 3; ++j)
											ei[j] -= electric_field_real[j][index] * Wxyz;
									}
									else
									{
										point3<T> Wxyz(
											dW[order*idx + kkk[0]][0] * W[order*idx + kkk[1]][1] * W[order*idx + kkk[2]][2],
											W[order*idx + kkk[0]][0] * dW[order*idx + kkk[1]][1] * W[order*idx + kkk[2]][2],
											W[order*idx + kkk[0]][0] * W[order*idx + kkk[1]][1] * dW[order*idx + kkk[2]][2]
										);
										ei -= zMG_real[index] * Wxyz;
									}
								}
						s.f[i] += s.z[i] * ei;
					}
				};
			for (size_t i = 0; i < num_threads; ++i)
				tp.enqueue(eval_f_k, i);
			tp.wait();
			// corrections
			eval_ewald_scd(s, kappa, volume);
		}

		private:

			std::vector<point3<T>> W, dW;
			std::vector<std::complex<T>> zMG, electric_field[3];
			std::vector<T> G, partU, partvir;
			std::vector<unsigned> part2mesh_inds, sorted_inds, first_inds, n_part_cell;
			math::dft<T> dft;
			std::mt19937_64 mersenne_twister = std::mt19937_64(0);
			std::uniform_real_distribution<T> u_dist = std::uniform_real_distribution<T>(0, 1);
			T *zMG_real, *electric_field_real[3];
			T kappa = 0.4, cutoff = 6;
			std::size_t order = 5;
			bool update = true, use_ik = false;
	};

} // namespace physics

#endif // PHYSICS_PPPM_H
































