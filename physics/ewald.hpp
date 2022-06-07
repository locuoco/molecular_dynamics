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

#include <cmath>
#include <thread>
#include <execution>
#include <numeric> // iota
#include <algorithm> // sort

#include "point.hpp"
#include "../math/helper.hpp" // fasterfc, fastexp
#include "../thread_pool.hpp"
#include "../math/dft.hpp"

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
		void operator()(S& s, std::size_t num_threads = 1)
		{
			using std::sqrt;
			using std::size_t;
			if (num_threads == 1)
			{
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
			else
			{
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
		}

		private:

			std::vector<T> partU, partvir;
			thread_pool<std::size_t> tp;
	};

	template <typename T, typename State>
	struct ewald
	// Ewald summation
	{
		template <coulomb_and_LJ_periodic<T, State> S>
		void operator()(S& s, std::size_t num_threads = 1, std::size_t maxn = 6)
		{
			using std::sqrt;
			using std::sin;
			using std::cos;
			using std::fabs;
			using std::size_t;
			size_t maxn1 = maxn+1, maxn3 = maxn1*maxn1*maxn1;
			T volume = s.side * s.side * s.side;
			T kappa = maxn/s.side;

			partU.resize(num_threads);
			partvir.resize(num_threads);

			k.resize(maxn3);
			factor.resize(maxn3);
			csum.resize(maxn3);
			ssum.resize(maxn3);

			auto eval_1 = [this, num_threads, maxn, maxn1, maxn3, volume, kappa, &s](size_t idx)
				{
					T uC = 0, u = 0, v = 0;
					size_t block = (s.n-1)/num_threads + 1;
					for (size_t i = idx*block; (i < s.n) && (i < (idx+1)*block); ++i)
						for (size_t j = 0; j < s.n; ++j)
							if (i != j)
							{
								T Rij = s.LJ_halfR[i] + s.LJ_halfR[j];

								auto r = remainder(s.x[i] - s.x[j], s.side);
								T r2 = dot(r, r);
								T r2_ = 1/dot(r, r);
								T d = sqrt(r2);
								T d_ = 1/d;
								T kd = kappa * d;
								T t = Rij * Rij * r2_;
								t = t * t * t;
								T epsijt = s.LJ_sqrteps[i] * s.LJ_sqrteps[j] * t;
								T zij_d = s.z[i] * s.z[j] * d_;
								T coulomb = zij_d * math::fasterfc(kd);
								T lennard_jones = 12 * epsijt * (t - 1);
								auto fij = ((lennard_jones + coulomb + zij_d * (2 / math::sqrtpi<T>()) * kd * math::fastexp(-kd*kd)) * r2_) * r;
								s.f[i] += fij;
								u += epsijt * (t - 2);
								uC += coulomb;
								v += lennard_jones;
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
			auto eval_2 = [this, num_threads, maxn1, maxn3, &s](size_t idx)
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
			point3<T> xsum = 0;
			T z2sum = 0, zsum = 0;
			for (size_t i = 0; i < s.n; ++i)
				xsum += s.z[i] * s.x[i];
			for (size_t i = 0; i < s.n; ++i)
				zsum += s.z[i];
			for (size_t i = 0; i < s.n; ++i)
				z2sum += s.z[i]*s.z[i];
			T factor = 2*math::two_pi<T>()/(3*volume);
			T u = factor * dot(xsum, xsum)/2 - kappa/math::sqrtpi<T>() * z2sum - math::pi<T>()/(2*kappa*kappa*volume)*fabs(zsum);
			s.U += u;
			s.virial += u;
			xsum *= factor;
			for (size_t i = 0; i < s.n; ++i)
				s.f[i] -= xsum * s.z[i];
		}

		private:

			std::vector<point3<T>> k;
			std::vector<T> partU, partvir;
			std::vector<T> factor, csum, ssum;
			thread_pool<std::size_t> tp;
	};

	/*template <typename T, typename State>
	struct pppm
	// Particle-particle, particle-mesh method (also known as PPPM or P^3M or P3M method)
	// based on Ewald summation, it has O(N log N) asymptotic complexity
	{
		template <coulomb_and_LJ_periodic<T, State> S>
		void operator()(S& s, std::size_t num_threads = 1)
		{
			using std::size_t;
			using std::floor;

			T eps = 1e-6;
			size_t mesh_size = 1 << 5;

			k.resize(s.n);
			indices.resize(s.n);

			for (size_t i = 0; i < s.n; ++i)
			{
				auto a = point3<size_t>((remainder(s.x[i]) / s.side + T(0.5) + eps)*mesh_size)
				k[i] = a[0]*mesh_size + a[1])*mesh_size + a[2];
			}
			auto comp = [this](size_t i, size_t j)
				{
					return k[i] < k[j];
				};
			std::iota(indices.begin(), indices.end(), 0);
			std::sort(std::execution::par_unseq, indices.begin(), indices.end(), comp);
			
		}

		private:

			std::vector<T> zM;
			std::vector<unsigned> k, indices;
	};*/

} // namespace physics

#endif // PHYSICS_EWALD_H
































