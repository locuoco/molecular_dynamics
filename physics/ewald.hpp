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
#include <numbers>

#include "point.hpp"
#include "integrator.hpp" // physical_system
#include "../math/helper.hpp" // fasterfc, fastexp
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

	template <typename T, typename S>
	void eval_ewald_r_6(S& s, T& u, T& uC, T& v, T kappa, std::size_t i, std::size_t j,
		decltype(S::x[i])& r, T r2)
	// real space summation for Ewald summation for both Coulomb and dispersion interactions
	{
		T Rij = s.LJ_halfR[i] + s.LJ_halfR[j];
		T epsij = s.LJ_sqrteps[i] * s.LJ_sqrteps[j];

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
		T r2_ = 1/r2, r4_ = r2_*r2_, r6_ = r4_*r2_;
		T k2 = kappa*kappa, k4 = k2*k2;
		T lj12 = C12ij * r6_*r6_;
		T lennard_jones = 12 * lj12 - C6ij * (6*r6_ + 6*k2*r4_ + 3*k4*r2_ + k4*k2) * emkd2;
		auto fij = ((lennard_jones + coulomb + zij_d * (2 * std::numbers::inv_sqrtpi_v<T>) * kd * emkd2) * r2_) * r;
		s.f[i] += fij;
		u += lj12 - C6ij * (r6_ + k2*r4_ + k4*r2_/2) * emkd2;
		uC += coulomb;
		v += lennard_jones;
	}

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

} // namespace physics

#endif // PHYSICS_EWALD_H
































