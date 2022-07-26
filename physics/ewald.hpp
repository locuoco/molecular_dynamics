//  Ewald summation method
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

#include <type_traits> // add_const_t, remove_reference_t
#include <cmath> // sqrt, sin, cos, remainder
#include <numbers> // numbers::inv_sqrtpi_v, numbers::pi_v

#include "tensor.hpp" // remainder
#include "physical_system.hpp"
#include "direct.hpp" // coulomb_and_lj
#include "../math/helper.hpp" // fasterfc, fastexp
#include "../utils/thread_pool.hpp"

namespace physics
{
	template <typename System, typename T, typename State>
	concept coulomb_and_lj_periodic = coulomb_and_lj<System, T, State>
		&& requires(System& s, std::size_t i)
	// A `coulomb_and_lj_periodic` system `s` is a `coulomb_and_lj` system so that the following
	// instructions are well-formed (compilable).
	// `side` is the side of the simulation box.
		{
			s.f[i] += remainder(s.x[i] - s.x[i], s.side);
		};

	template <typename T, typename System>
	void eval_ewald_r(System& s, T& u, T& uC, T& v, T kappa, std::size_t i, std::size_t j,
		std::add_const_t<std::remove_reference_t<decltype(System::x[i])>>& r, T r2)
	// real space summation
	{
		using std::sqrt;

		T Rij = s.lj_halfR[i] + s.lj_halfR[j];

		T r2_ = 1/r2;
		T d = sqrt(r2);
		T d_ = 1/d;
		T kd = kappa * d;
		T t = Rij * Rij * r2_;
		t = t * t * t;
		T epsijt = s.lj_sqrteps[i] * s.lj_sqrteps[j] * t;
		T zij_d = s.z[i] * s.z[j] * d_;
		T coulomb = zij_d * math::fasterfc(kd);
		T lennard_jones = 12 * epsijt * (t - 1);
		auto fij = ((lennard_jones + coulomb + zij_d * (2 * std::numbers::inv_sqrtpi_v<T>) * kd * math::fastexp(-kd*kd)) * r2_) * r;
		s.f[i] += fij;
		u += epsijt * (t - 2);
		uC += coulomb;
		v += lennard_jones;
	}

	template <typename System, typename T>
	void eval_ewald_r(System& s, T& u, T& uC, T& v, T kappa, std::size_t i, std::size_t j)
	// real space summation
	{
		auto r = remainder(s.x[i] - s.x[j], s.side);
		T r2 = dot(r, r);
		eval_ewald_r(s, u, uC, v, kappa, i, j, r, r2);
	}

	template <typename System, typename T>
	requires requires(System& s, T x)
		{
			x += s.Z;
			x += s.Z2;
		}
	void eval_ewald_scd(System& s, T kappa, T volume, T dielectric = 1)
	// self-energy, charged system and dipole corrections
	{
		using std::size_t;
		vec3<T> xsum = 0;
		for (size_t i = 0; i < s.n; ++i)
			xsum += s.z[i] * s.x[i];
		T fact = math::two_pi<T>/((1+2*dielectric)*volume);
		T u = fact * dot(xsum, xsum) - kappa*std::numbers::inv_sqrtpi_v<T> * s.Z2 - std::numbers::pi_v<T>/(2*kappa*kappa*volume)*s.Z*s.Z;
		s.potential += u;
		s.virial += u;
		xsum *= fact*2;
		for (size_t i = 0; i < s.n; ++i)
			s.f[i] -= xsum * s.z[i];
	}

	template <typename T, typename State>
	struct ewald
	// Ewald summation
	{
		T cutoff_radius() const noexcept
		{
			return cutoff;
		}

		void max_n(std::size_t maxn) noexcept
		{
			this -> maxn = maxn;
		}

		std::size_t max_n() const noexcept
		{
			return maxn;
		}

		template <coulomb_and_lj_periodic<T, State> System>
		void operator()(System& s, utils::thread_pool& tp)
		{
			using std::sin;
			using std::cos;
			using std::size_t;
			size_t maxn1 = maxn+1, maxn3 = maxn1*maxn1*maxn1;
			T volume = s.side * s.side * s.side;
			cutoff = s.side/2;
			T cutoff2 = cutoff*cutoff;
			T kappa = 3/cutoff;

			auto num_threads = tp.size();
			partpot.resize(num_threads);
			partvir.resize(num_threads);

			k.resize(maxn3);
			factor.resize(maxn3);
			csum.resize(maxn3);
			ssum.resize(maxn3);

			auto eval_1 = [this, num_threads, maxn1, maxn3, volume, kappa, cutoff2, &s](size_t idx)
				{
					T uC = 0, uD = 0, u = 0, v = 0;
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
					uD /= 2;
					v /= 2;
					block = (maxn3-1)/num_threads + 1;
					for (size_t nijk = idx*block+1; (nijk < maxn3) && (nijk < (idx+1)*block+1); ++nijk)
					{
						size_t n3 = nijk % maxn1;
						size_t nij = nijk / maxn1;
						size_t n2 = nij % maxn1;
						size_t n1 = nij / maxn1;
						k[nijk-1] = math::two_pi<T>*vec3<T>{n1, n2, n3}/s.side;
						T k2 = dot(k[nijk-1], k[nijk-1]);
						csum[nijk-1] = 0, ssum[nijk-1] = 0;
						for (size_t i = 0; i < s.n; ++i)
						{
							T dotkx = dot(k[nijk-1], s.x[i]);
							csum[nijk-1] += s.z[i] * cos(dotkx);
							ssum[nijk-1] += s.z[i] * sin(dotkx);
						}
						factor[nijk-1] = 4*math::two_pi<T>*math::fastexp(-k2/(4*kappa*kappa))/(k2*volume);
						T kspace_contrib = (csum[nijk-1]*csum[nijk-1] + ssum[nijk-1]*ssum[nijk-1])*factor[nijk-1]/2;
						uC += kspace_contrib;
					}
					partpot[idx] = u + uC + uD;
					partvir[idx] = v + uC + 6*uD;
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
			f = s.f;
			for (size_t i = 0; i < num_threads; ++i)
				tp.enqueue(eval_1, i);
			tp.wait();
			for (size_t i = 0; i < num_threads; ++i)
				tp.enqueue(eval_2, i);
			tp.wait();
			for (size_t i = 0; i < num_threads; ++i)
			{
				s.potential += partpot[i];
				s.virial += partvir[i];
			}
			eval_ewald_scd(s, kappa, volume);
			f = s.f - f;
		}

		state<T, 3> f;

		private:

			std::vector<vec3<T>> k;
			std::vector<T> partpot, partvir;
			std::vector<T> factor, csum, ssum;
			T cutoff = 0;
			std::size_t maxn = 6;
	};

} // namespace physics

#endif // PHYSICS_EWALD_H
































