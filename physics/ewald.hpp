//  Ewald summation methods
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

#ifndef PHYSICS_EWALD_H
#define PHYSICS_EWALD_H

#include <cmath> // round, ceil, log2, sqrt, sin, cos, fabs, remainder
#include <numeric> // iota
#include <numbers>
#include <functional> // function
#include <algorithm> // fill, copy
#include <complex>
#include <stdexcept> // runtime_error

#include "point.hpp"
#include "../math/helper.hpp" // fasterfc, fastexp, mod
#include "../math/dft.hpp"
#include "../utils/thread_pool.hpp"
#include "../utils/parallel_sort.hpp"

namespace physics
{
	template <typename S, typename T, typename State>
	concept coulomb_and_LJ = physical_system<S, T, State>
		&& requires(S& s, T t, std::size_t i)
		{
			i < s.n;
			t += s.z[i];
			t += s.LJ_sqrteps[i];
			t += s.LJ_halfR[i];
			s.f[i] += s.x[i] - s.x[i];
			s.U += t;
			s.virial += dot(s.x[i], s.x[i]);
		};

	template <typename S, typename T, typename State>
	concept coulomb_and_LJ_periodic = coulomb_and_LJ<S, T, State>
		&& requires(S& s, std::size_t i)
		{
			s.f[i] += remainder(s.x[i] - s.x[i], s.side);
		};

	template <typename T, typename State>
	struct direct
	// Direct summation
	{
		template <coulomb_and_LJ<T, State> S>
		void operator()(S& s)
		{
			using std::sqrt;
			using std::size_t;
			for (size_t i = 0; i < s.n; ++i)
				for (size_t j = i+1; j < s.n; ++j)
				{
					T Rij = s.LJ_halfR[i] + s.LJ_halfR[j];

					auto r = s.x[i] - s.x[j];
					T r2_ = 1/dot(r, r);
					T d_ = sqrt(r2_);
					T t = Rij * Rij * r2_;
					t = t * t * t;
					T epsijt = s.LJ_sqrteps[i] * s.LJ_sqrteps[j] * t;
					T coulomb = s.z[i] * s.z[j] * d_;
					T tot_virial = 12 * epsijt * (t - 1) + coulomb;
					auto fij = (tot_virial * r2_) * r;
					s.f[i] += fij;
					s.f[j] -= fij;
					s.U += epsijt * (t - 2) + coulomb;
					s.virial += tot_virial;
				}
		}

		template <coulomb_and_LJ<T, State> S>
		void operator()(S& s, utils::thread_pool& tp)
		{
			using std::sqrt;
			using std::size_t;
			auto num_threads = tp.size();
			if (num_threads == 1)
				return operator()(s);
			partU.resize(num_threads);
			partvir.resize(num_threads);
			auto eval_lambda = [this, num_threads, &s](size_t idx)
				{
					T u = 0, v = 0;
					size_t block = (s.n-1)/num_threads + 1;
					for (size_t i = idx*block; (i < s.n) && (i < (idx+1)*block); ++i)
						for (size_t j = 0; j < s.n; ++j)
							if (i != j)
							{
								T Rij = s.LJ_halfR[i] + s.LJ_halfR[j];

								auto r = s.x[i] - s.x[j];
								T r2_ = 1/dot(r, r);
								T d_ = sqrt(r2_);
								T t = Rij * Rij * r2_;
								t = t * t * t;
								T epsijt = s.LJ_sqrteps[i] * s.LJ_sqrteps[j] * t;
								T coulomb = s.z[i] * s.z[j] * d_;
								T tot_virial = 12 * epsijt * (t - 1) + coulomb;
								auto fij = (tot_virial * r2_) * r;
								s.f[i] += fij;
								u += epsijt * (t - 2) + coulomb;
								v += tot_virial;
							}
					partU[idx] = u/2;
					partvir[idx] = v/2;
				};
			for (size_t i = 0; i < num_threads; ++i)
				tp.enqueue(eval_lambda, i);
			tp.wait();
			for (size_t i = 0; i < num_threads; ++i)
			{
				s.U += partU[i];
				s.virial += partvir[i];
			}
		}

		private:

			std::vector<T> partU, partvir;
	};

	template <typename T, typename S>
	void eval_ewald_r(S& s, T& u, T& uC, T& v, T kappa, std::size_t i, std::size_t j,
		decltype(S::x[i])& r, T r2)
	// real space summation
	{
		T Rij = s.LJ_halfR[i] + s.LJ_halfR[j];

		T r2_ = 1/r2;
		T d = sqrt(r2);
		T d_ = 1/d;
		T kd = kappa * d;
		T t = Rij * Rij * r2_;
		t = t * t * t;
		T epsijt = s.LJ_sqrteps[i] * s.LJ_sqrteps[j] * t;
		T zij_d = s.z[i] * s.z[j] * d_;
		T coulomb = zij_d * math::fasterfc(kd);
		T lennard_jones = 12 * epsijt * (t - 1);
		auto fij = ((lennard_jones + coulomb + zij_d * (2 * std::numbers::inv_sqrtpi_v<T>) * kd * math::fastexp(-kd*kd)) * r2_) * r;
		s.f[i] += fij;
		u += epsijt * (t - 2);
		uC += coulomb;
		v += lennard_jones;
	}

	/*template <typename T, typename S>
	void eval_ewald_r_lj(S& s, T& u, T& uC, T& v, T kappa, std::size_t i, std::size_t j,
		decltype(S::x[i])& r, T r2)
	// real space summation for Ewald summation for both Coulomb and LJ interactions
	{
		T Rij = s.LJ_halfR[i] + s.LJ_halfR[j];
		T epsij = s.LJ_sqrteps[i] * s.LJ_sqrteps[j];

		T r2_ = 1/r2, r4_ = r2_*r2_, r6_ = r4_*r2_, r12_ = r6_*r6_;
		T d = sqrt(r2);
		T d_ = 1/d;
		T kd = kappa * d;
		T emkd2 = math::fastexp(-kd*kd);
		T C6ij = Rij * Rij;
		C6ij = C6ij * C6ij * C6ij;
		T C12ij = C6ij * C6ij;
		C6ij *= 2*epsij;
		C12ij *= epsij;
		T zij_d = s.z[i] * s.z[j] * d_;
		T coulomb = zij_d * math::fasterfc(kd);
		T lennard_jones = 12 * C12ij * r12_ - C6ij * ();
		auto fij = ((lennard_jones + coulomb + zij_d * (2 * std::numbers::inv_sqrtpi_v<T>) * kd * emkd2) * r2_) * r;
		s.f[i] += fij;
		u += epsijt * (t - 2);
		uC += coulomb;
		v += lennard_jones;
	}*/

	template <typename S, typename T>
	void eval_ewald_r(S& s, T& u, T& uC, T& v, T kappa, std::size_t i, std::size_t j)
	// real space summation
	{
		auto r = remainder(s.x[i] - s.x[j], s.side);
		T r2 = dot(r, r);
		eval_ewald_r(s, u, uC, v, kappa, i, j, r, r2);
	}

	template <typename S, typename T>
	void eval_ewald_scd(S& s, T kappa, T volume, T dielectric = 1)
	// self-energy, charged system and dipole corrections
	{
		using std::size_t;
		point3<T> xsum = 0;
		T z2sum = 0, zsum = 0;
		for (size_t i = 0; i < s.n; ++i)
			xsum += s.z[i] * s.x[i];
		for (size_t i = 0; i < s.n; ++i)
			zsum += s.z[i];
		for (size_t i = 0; i < s.n; ++i)
			z2sum += s.z[i]*s.z[i];
		T fact = math::two_pi<T>()/((1+2*dielectric)*volume);
		T u = fact * dot(xsum, xsum) - kappa*std::numbers::inv_sqrtpi_v<T> * z2sum - std::numbers::pi_v<T>/(2*kappa*kappa*volume)*zsum*zsum;
		s.U += u;
		s.virial += u;
		xsum *= fact*2;
		for (size_t i = 0; i < s.n; ++i)
			s.f[i] -= xsum * s.z[i];
	}

	template <typename T, typename State>
	struct ewald
	// Ewald summation
	{
		template <coulomb_and_LJ_periodic<T, State> S>
		void operator()(S& s, utils::thread_pool& tp, std::size_t maxn = 6)
		{
			using std::sqrt;
			using std::sin;
			using std::cos;
			using std::fabs;
			using std::size_t;
			size_t maxn1 = maxn+1, maxn3 = maxn1*maxn1*maxn1;
			T volume = s.side * s.side * s.side;
			T kappa = maxn/s.side;
			T cutoff = s.side/2, cutoff2 = cutoff*cutoff;

			auto num_threads = tp.size();
			partU.resize(num_threads);
			partvir.resize(num_threads);

			k.resize(maxn3);
			factor.resize(maxn3);
			csum.resize(maxn3);
			ssum.resize(maxn3);

			auto eval_1 = [this, num_threads, maxn1, maxn3, volume, kappa, cutoff2, &s](size_t idx)
				{
					T uC = 0, u = 0, v = 0;
					size_t block = (s.n-1)/num_threads + 1;
					for (size_t i = idx*block; (i < s.n) && (i < (idx+1)*block); ++i)
						for (size_t j = 0; j < s.n; ++j)
							if (i != j)
							{
								auto r = remainder(s.x[i] - s.x[j], s.side);
								T r2 = dot(r, r);
								if (r2 <= cutoff2)
									eval_ewald_r(s, u, uC, v, kappa, i, j, r, r2);
							}
					u /= 2;
					uC /= 2;
					v /= 2;
					block = (maxn3-1)/num_threads + 1;
					for (size_t nijk = idx*block+1; (nijk < maxn3) && (nijk < (idx+1)*block+1); ++nijk)
					{
						size_t n3 = nijk % maxn1;
						size_t nij = nijk / maxn1;
						size_t n2 = nij % maxn1;
						size_t n1 = nij / maxn1;
						k[nijk-1] = math::two_pi<T>()*point3<T>{n1, n2, n3}/s.side;
						T k2 = dot(k[nijk-1], k[nijk-1]);
						csum[nijk-1] = 0, ssum[nijk-1] = 0;
						for (size_t i = 0; i < s.n; ++i)
						{
							T dotkx = dot(k[nijk-1], s.x[i]);
							csum[nijk-1] += s.z[i] * cos(dotkx);
							ssum[nijk-1] += s.z[i] * sin(dotkx);
						}
						factor[nijk-1] = 4*math::two_pi<T>()*math::fastexp(-k2/(4*kappa*kappa))/(k2*volume);
						T kspace_contrib = (csum[nijk-1]*csum[nijk-1] + ssum[nijk-1]*ssum[nijk-1])*factor[nijk-1]/2;
						uC += kspace_contrib;
					}
					partU[idx] = u + uC;
					partvir[idx] = v + uC;
				};
			auto eval_2 = [this, num_threads, maxn3, &s](size_t idx)
				{
					size_t block = (s.n-1)/num_threads + 1;
					for (size_t i = idx*block; (i < s.n) && (i < (idx+1)*block); ++i)
					{
						for (size_t nijk = 1; nijk < maxn3; ++nijk)
						{
							T dotkx = dot(k[nijk-1], s.x[i]);
							s.f[i] += (factor[nijk-1] * s.z[i] * (csum[nijk-1] * sin(dotkx) - ssum[nijk-1] * cos(dotkx))) * k[nijk-1];
						}
					}
				};
			for (size_t i = 0; i < num_threads; ++i)
				tp.enqueue(eval_1, i);
			tp.wait();
			for (size_t i = 0; i < num_threads; ++i)
				tp.enqueue(eval_2, i);
			tp.wait();
			for (size_t i = 0; i < num_threads; ++i)
			{
				s.U += partU[i];
				s.virial += partvir[i];
			}
			eval_ewald_scd(s, kappa, volume);
		}

		private:

			std::vector<point3<T>> k;
			std::vector<T> partU, partvir;
			std::vector<T> factor, csum, ssum;
	};

	static constexpr std::size_t pppm_max_order = 7;

	template <typename T>
	const std::function<T(const std::vector<T>&)> charge_assignment_function[pppm_max_order][pppm_max_order] =
		{
			// P = 1
			{[](const std::vector<T>&) { return T(1); }},
			// P = 2
			{
				[](const std::vector<T>& x) { return T(0.5) - x[1]; },
				[](const std::vector<T>& x) { return T(0.5) + x[1]; },
			},
			// P = 3
			{
				[](const std::vector<T>& x) { T s(T(0.5) - x[1]); return s*s/2; },
				[](const std::vector<T>& x) { return T(0.75) - x[2]; },
				[](const std::vector<T>& x) { T s(T(0.5) + x[1]); return s*s/2; },
			},
			// P = 4
			{
				[](const std::vector<T>& x) { T s(T(0.5) - x[1]); return s*s*s/6; },
				[](const std::vector<T>& x) { return (23 - 30*x[1] - 12*x[2] + 24*x[3])/48; },
				[](const std::vector<T>& x) { return (23 + 30*x[1] - 12*x[2] - 24*x[3])/48; },
				[](const std::vector<T>& x) { T s(T(0.5) + x[1]); return s*s*s/6; },
			},
			// P = 5
			{
				[](const std::vector<T>& x) { T s(T(0.5) - x[1]), s2(s*s); return s2*s2/24; },
				[](const std::vector<T>& x) { return (19 - 44*x[1] + 24*x[2] + 16*x[3] - 16*x[4])/96; },
				[](const std::vector<T>& x) { return (115 - 120*x[2] + 48*x[4])/192; },
				[](const std::vector<T>& x) { return (19 + 44*x[1] + 24*x[2] - 16*x[3] - 16*x[4])/96; },
				[](const std::vector<T>& x) { T s(T(0.5) + x[1]), s2(s*s); return s2*s2/24; },
			},
			// P = 6
			{
				[](const std::vector<T>& x) { T s(T(0.5) - x[1]), s2(s*s); return s2*s2*s/120; },
				[](const std::vector<T>& x) { return (237 - 750*x[1] + 840*x[2] - 240*x[3] - 240*x[4] + 160*x[5])/3840; },
				[](const std::vector<T>& x) { return (841 - 770*x[1] - 440*x[2] + 560*x[3] + 80*x[4] - 160*x[5])/1920; },
				[](const std::vector<T>& x) { return (841 + 770*x[1] - 440*x[2] - 560*x[3] + 80*x[4] + 160*x[5])/1920; },
				[](const std::vector<T>& x) { return (237 + 750*x[1] + 840*x[2] + 240*x[3] - 240*x[4] - 160*x[5])/3840; },
				[](const std::vector<T>& x) { T s(T(0.5) + x[1]), s2(s*s); return s2*s2*s/120; },
			},
			// P = 7
			{
				[](const std::vector<T>& x) { T s(T(0.5) - x[1]), s3(s*s*s); return s3*s3/720; },
				[](const std::vector<T>& x) { return (361 - 1416*x[1] + 2220*x[2] - 1600*x[3] + 240*x[4] + 384*x[5] - 192*x[6])/23040; },
				[](const std::vector<T>& x) { return (10543 - 17340*x[1] + 4740*x[2] + 6880*x[3] - 4080*x[4] - 960*x[5] + 960*x[6])/46080; },
				[](const std::vector<T>& x) { return (5887 - 4620*x[2] + 1680*x[4] - 320*x[6])/11520; },
				[](const std::vector<T>& x) { return (10543 + 17340*x[1] + 4740*x[2] - 6880*x[3] - 4080*x[4] + 960*x[5] + 960*x[6])/46080; },
				[](const std::vector<T>& x) { return (361 + 1416*x[1] + 2220*x[2] + 1600*x[3] + 240*x[4] - 384*x[5] - 192*x[6])/23040; },
				[](const std::vector<T>& x) { T s(T(0.5) + x[1]), s3(s*s*s); return s3*s3/720; },
			},
		};

	template <typename T>
	const std::function<T(T)> influence_function_den[pppm_max_order] =
		{
			[](T) { return T(1); }, // P = 1
			[](T s2) { return 1 - T(2.L/3)*s2; }, // P = 2
			[](T s2) { return 1 - s2 + T(2.L/15)*s2*s2; }, // P = 3
			[](T s2) { T s4(s2*s2); return 1 - T(4.L/3)*s2 + T(2.L/5)*s4 - T(4.L/315)*s4*s2; }, // P = 4
			[](T s2) { T s4(s2*s2); return 1 - T(5.L/3)*s2 + T(7.L/9)*s4 - T(17.L/189)*s4*s2 + T(2.L/2835)*s4*s4; }, // P = 5
			[](T s2) { T s4(s2*s2), s8(s4*s4); return 1 - 2*s2 + T(19.L/15)*s4 - T(256.L/945)*s4*s2 + T(62.L/4725)*s8 - T(4.L/155925)*s8*s2; }, // P = 6
			[](T s2)
			{
				T s4(s2*s2), s8(s4*s4);
				return 1 - T(7.L/3)*s2 + T(28.L/15)*s4 - T(16.L/27)*s4*s2 + T(26.L/405)*s8 - T(2.L/1485)*s8*s2 + T(4.L/6081075)*s8*s4;
			}, // P = 7
		};

	template <typename T, typename State>
	struct pppm
	// Particle-particle, particle-mesh method (also known as PPPM or P^3M or P3M method)
	// based on Ewald summation, it has O(N log N) asymptotic complexity
	{
		template <coulomb_and_LJ_periodic<T, State> S>
		void operator()(S& s, utils::thread_pool& tp, T cutoff = 7, T kappa = .4, std::size_t order = 7)
		// cutoff is the cutoff radius in real space, while kappa is the Ewald parameter
		// order is the order of the assignment scheme
		{
			using std::size_t;
			using std::ptrdiff_t;
			using std::ceil;
			using std::sin;
			using std::remainder;
			if (s.n < 2)
				return;
			if (order < 1)
				throw std::runtime_error("Error: minimum PPPM order is 1");
			if (order > pppm_max_order)
				throw std::runtime_error("Error: maximum PPPM order is 7");

			size_t num_cells = 1ull << size_t(std::round(std::log2(s.n)/3));
			size_t m = num_cells*num_cells*num_cells;
			size_t num_threads = tp.size();
			T cell_size = s.side / num_cells, cell_size3 = cell_size*cell_size*cell_size;
			T volume = s.side * s.side * s.side;
			cutoff = std::min(cutoff, s.side/2);

			partU.resize(num_threads);
			partvir.resize(num_threads);
			part2mesh_inds.resize(s.n);
			sorted_inds.resize(s.n);
			first_inds.resize(m);
			n_part_cell.resize(m);
			zM.resize(m);
			zG.resize(m);
			G.resize(m);
			W.resize(num_threads*order);
			x_prime.resize(3*num_threads);
			for (size_t i = 0; i < 3*num_threads; ++i)
				x_prime[i].resize(order);
			for (size_t i = 0; i < 3; ++i)
				electric_field[i].resize(m);

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

			auto vec2index = [num_cells](point3<ptrdiff_t> vec)
				{
					return (math::mod(vec[0], ptrdiff_t(num_cells))*num_cells +
							math::mod(vec[1], ptrdiff_t(num_cells)))*num_cells +
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
			// calculate real-space contribution to force/energy
			auto eval_r = [this, num_threads, kappa, cutoff, cell_size, num_cells, ki, index2vec, vec2index, &s](size_t idx)
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
									size_t first_ind = first_inds[other_k], npc = n_part_cell[other_k];
									// iterating through all particles in the cell
									for (size_t j = first_ind; j < first_ind+npc; ++j)
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
			// calculate mesh charges (zM)
			auto charge_assignment_center = [order, cell_size, num_cells, index2vec, &s](size_t i, size_t p)
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
			std::fill(zM.begin(), zM.end(), 0);
			for (size_t i = 0; i < s.n; ++i)
			{
				ptrdiff_t p = part2mesh_inds[i];
				point3<T> center(charge_assignment_center(i, p));
				if (order > 1)
				{
					point3<T> xp = remainder(s.x[i]/cell_size - center, T(num_cells));
					for (size_t j = 0; j < 3; ++j)
						x_prime[j][1] = xp[j];
					for (size_t j = 0; j < 3; ++j)
						for (size_t k = 2; k < order; ++k)
							x_prime[j][k] = x_prime[j][k-1] * xp[j];
				}
				for (size_t k = 0; k < order; ++k)
					for (size_t j = 0; j < 3; ++j)
						W[k][j] = charge_assignment_function<T>[order-1][k](x_prime[j]);
				point3<ptrdiff_t> kkk, cell(mod(center, T(num_cells)));
				ptrdiff_t half_order = (order-1)/2;
				cell -= half_order;
				for (kkk[0] = 0; kkk[0] < ptrdiff_t(order); ++kkk[0])
					for (kkk[1] = 0; kkk[1] < ptrdiff_t(order); ++kkk[1])
						for (kkk[2] = 0; kkk[2] < ptrdiff_t(order); ++kkk[2])
						{
							size_t index = vec2index(cell + kkk);
							zM[index] += s.z[i] * W[kkk[0]][0] * W[kkk[1]][1] * W[kkk[2]][2];
						}
			}
			tp.wait();
			for (size_t i = 0; i < num_threads; ++i)
			{
				s.U += partU[i];
				s.virial += partvir[i];
			}
			// calculate influence function
			auto eval_G = [this, num_threads, kappa, order, num_cells, cell_size, m, index2vec, &s](size_t idx)
				{
					T pi_h = std::numbers::pi_v<T>/cell_size, twopi_h = 2*pi_h, twopi_L = twopi_h/num_cells;
					T m1_4kappa2 = -1/(4*kappa*kappa);
					point3<ptrdiff_t> nnn_;
					ptrdiff_t maxn_ = 2;
					size_t block = (m-1)/num_threads + 1;
					for (size_t i = idx*block; (i < m) && (i < (idx+1)*block); ++i)
					{
						if (i == 0)
						{
							G[i] = 0;
							continue;
						}
						auto nvec = index2vec(i);
						nvec -= ptrdiff_t(num_cells)*(2*nvec/ptrdiff_t(num_cells));
						point3<T> kvec = twopi_L * point3<T>(nvec), s = kvec*cell_size/2;
						T k2 = dot(kvec, kvec);
						for (size_t j = 0; j < 3; ++j)
							s[j] = sin(s[j]);
						T den = 1, num = 0;
						for (size_t j = 0; j < 3; ++j)
							den *= influence_function_den<T>[order-1](s[j]*s[j]);
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
									for (size_t k = 0; k < order; ++k)
										U2 *= s_[2];
									U2 *= U2;
									num += dot(kvec, k_) * U2 / k_2 * math::fastexp(k_2*m1_4kappa2);
								}
							}
						}
						G[i] = num * 2 * math::two_pi<T>() / (k2*den*den);
					}
				};
			if (update)
			{
				for (size_t i = 0; i < num_threads; ++i)
					tp.enqueue(eval_G, i);
				tp.wait();
				update = false;
			}
			// calculate reciprocal-space energy
			std::array<size_t, 3> nlist{num_cells, num_cells, num_cells};
			for (size_t i = 0; i < m; ++i)
				zG[i] = zM[i];
			dft.template fftn<3>(zG, nlist, tp);
			{
				T u{};
				for (size_t i = 0; i < m; ++i)
					u += norm(zG[i]) * G[i];
				u *= 1/(2*volume);
				s.U += u;
				s.virial += u;
			}
			// calculate reciprocal-space forces
			for (size_t i = 0; i < m; ++i)
				zG[i] *= G[i] / cell_size3;
			T twopi_L = math::two_pi<T>()/s.side;
			point3<size_t> nnn;
			point3<T> kvec;
			for (nnn[0] = 0; nnn[0] < num_cells; ++nnn[0])
				for (nnn[1] = 0; nnn[1] < num_cells; ++nnn[1])
					for (nnn[2] = 0; nnn[2] < num_cells; ++nnn[2])
					{
						size_t index = (nnn[0]*num_cells + nnn[1])*num_cells + nnn[2];
						point3<ptrdiff_t> nvec(nnn);
						for (size_t j = 0; j < 3; ++j)
							if (nvec[j] == ptrdiff_t(num_cells)/2)
								nvec[j] = 0;
						nvec -= ptrdiff_t(num_cells)*(2*nvec/ptrdiff_t(num_cells));
						kvec = twopi_L * point3<T>(nvec);
						for (size_t j = 0; j < 3; ++j)
							electric_field[j][index] = kvec[j] * zG[index];
					}
			for (size_t j = 0; j < 3; ++j)
				dft.template ifftn<3>(electric_field[j], nlist, tp);
			auto eval_f_k = [this, num_threads, order, cell_size, num_cells, vec2index, &s, &charge_assignment_center] (size_t idx)
				{
					size_t block = (s.n-1)/num_threads + 1;
					for (size_t i = idx*block; (i < s.n) && (i < (idx+1)*block); ++i)
					{
						ptrdiff_t p = part2mesh_inds[i];
						point3<T> center(charge_assignment_center(i, p));
						if (order > 1)
						{
							point3<T> xp = remainder(s.x[i]/cell_size - center, T(num_cells));
							for (size_t j = 0; j < 3; ++j)
								x_prime[3*idx + j][1] = xp[j];
							for (size_t j = 0; j < 3; ++j)
								for (size_t k = 2; k < order; ++k)
									x_prime[3*idx + j][k] = x_prime[3*idx + j][k-1] * xp[j];
						}
						for (size_t k = 0; k < order; ++k)
							for (size_t j = 0; j < 3; ++j)
								W[order*idx + k][j] = charge_assignment_function<T>[order-1][k](x_prime[3*idx + j]);
						point3<T> ei(0);
						point3<ptrdiff_t> kkk, cell(mod(center, T(num_cells)));
						ptrdiff_t half_order = (order-1)/2;
						cell -= half_order;
						for (kkk[0] = 0; kkk[0] < ptrdiff_t(order); ++kkk[0])
							for (kkk[1] = 0; kkk[1] < ptrdiff_t(order); ++kkk[1])
								for (kkk[2] = 0; kkk[2] < ptrdiff_t(order); ++kkk[2])
								{
									size_t index = vec2index(cell + kkk);
									T Wxyz = W[order*idx + kkk[0]][0] * W[order*idx + kkk[1]][1] * W[order*idx + kkk[2]][2];
									for (size_t j = 0; j < 3; ++j)
										ei[j] += electric_field[j][index].imag() * Wxyz;
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

			std::vector<point3<T>> W;
			std::vector<std::vector<T>> x_prime;
			std::vector<std::complex<T>> zG, electric_field[3];
			std::vector<T> zM, G, partU, partvir;
			std::vector<unsigned> part2mesh_inds, sorted_inds, first_inds, n_part_cell;
			math::dft<T> dft;
			bool update = true;
	};

} // namespace physics

#endif // PHYSICS_EWALD_H
































