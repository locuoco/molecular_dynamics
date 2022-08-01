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
#include <complex>

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
	void eval_ewald_r(System& s, T& uC, T& uLJ, T& vLJ, T kappa, std::size_t i, std::size_t j,
		std::add_const_t<std::remove_reference_t<decltype(System::x[i])>>& r, T r2)
	// Calculate the real space part of Ewald summation between two particles `i` and `j`.
	// `s` is the physical system (`coulomb_and_lj` constraints required).
	// `uC` is a variable which accumulates the electrostatic potential energy.
	// `uLJ` is a variable which accumulates the Lennard-Jones potential energy.
	// `vLJ` is a variable which accumulates the Lennard-Jones virial.
	// `kappa` is the Ewald parameter.
	// `r` is the displacement vector between the two particles, given by:
	//	r = remainder(s.x[i] - s.x[j])
	// `r2` is the squared distance between the two particles.
	// Both Coulomb and Lennard-Jones potentials and forces are calculated at
	// the same time, but Ewald summation is performed only for Coulomb potential.
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
		uC += coulomb;
		uLJ += epsijt * (t - 2);
		vLJ += lennard_jones;
	}

	template <typename T, typename System>
	void eval_ewald_r(System& s, T& uC, T& uLJ, T& vLJ, T kappa, std::size_t i, std::size_t j)
	// Calculate the real space part of Ewald summation between two particles `i` and `j`.
	// `s` is the physical system (`coulomb_and_lj` constraints required).
	// `uC` is a variable which accumulates the electrostatic potential energy.
	// `uLJ` is a variable which accumulates the Lennard-Jones potential energy.
	// `vLJ` is a variable which accumulates the Lennard-Jones virial.
	// `kappa` is the Ewald parameter.
	// This overload calculates the displacement `r` and its square `r2` between the two
	// particles and then calls eval_ewald_r(s, uC, uLJ, v, kappa, i, j, r, r2).
	{
		auto r = remainder(s.x[i] - s.x[j], s.side);
		T r2 = dot(r, r);
		eval_ewald_r(s, uC, uLJ, vLJ, kappa, i, j, r, r2);
	}

	template <typename T, typename System>
	requires requires(System& s, T x)
		{
			x += s.Z;
			x += s.Z2;
		}
	void eval_ewald_scd(System& s, T& uC, T kappa, T volume, T dielectric = 1)
	// Calculate self-energy, charged system and dipole corrections.
	// `s` is the physical system (`coulomb_and_lj` constraints required). Additional
	// requirements are s.Z which is the total charge of the system and s.Z2 which is the sum
	// square of the charges (after rescaling the charges by the square root of the Coulomb
	// constant).
	// `uC` is a variable which accumulates the electrostatic potential energy.
	// `kappa` is the Ewald parameter.
	// `volume` is the volume of the system.
	// `dielectric` is the relative dielectric constant (1 for vacuum).
	{
		using std::size_t;
		vec3<T> xsum = 0;
		for (size_t i = 0; i < s.n; ++i)
			xsum += s.z[i] * s.x[i];
		T fact = math::two_pi<T>/((1+2*dielectric)*volume);
		T u = fact*dot(xsum, xsum); // dipole correction
		u += -kappa*std::numbers::inv_sqrtpi_v<T> * s.Z2; // self-energy correction
		u += -std::numbers::pi_v<T>/(2*kappa*kappa*volume)*s.Z*s.Z; // charged system correction
		uC += u;
		s.virial += u;
		xsum *= fact*2;
		for (size_t i = 0; i < s.n; ++i)
			s.f[i] -= xsum * s.z[i];
	}

	template <typename T, typename State>
	struct ewald
	// Ewald summation for electrostatic potential. Lennard-Jones potential instead
	// is truncated at the cutoff radius (which is set as the same as the Ewald cutoff
	// parameter for convenience).
	{
		T cutoff_radius() const noexcept
		// Return the cutoff radius (Ewald cutoff parameter and LJ truncation radius
		// are the same).
		{
			return cutoff;
		}

		void max_n(std::size_t maxn) noexcept
		// Set the maximum reciprocal number, related to the maximum reciprocal vector:
		// 	kmax = 2 pi max_n / side
		{
			this->maxn = maxn;
		}

		std::size_t max_n() const noexcept
		// Return the maximum reciprocal number, related to the maximum reciprocal vector:
		// 	kmax = 2 pi max_n / side
		{
			return maxn;
		}

		template <coulomb_and_lj_periodic<T, State> System>
		void operator()(System& s, utils::thread_pool& tp)
		// Perform Ewald summation with multi-threading.
		// `s` is the physical system.
		// `tp` is a thread pool.
		{
			using std::size_t;

			maxn1 = 2*maxn+1;
			maxn3 = maxn1*maxn1*maxn1;
			T volume = s.side * s.side * s.side;
			cutoff = s.side/2;
			kappa = 7/s.side;

			num_threads = tp.size();
			partpot_coulomb.resize(num_threads);
			partpot_lj.resize(num_threads);
			partvir.resize(num_threads);

			k.resize(maxn3);
			factor.resize(maxn3);
			structure.resize(maxn3);

			f = s.f;
			energy_coulomb = energy_lj = 0;
			for (size_t i = 0; i < num_threads; ++i)
				tp.enqueue([this, i, &s] { eval_1(s, i); });
			tp.wait();
			for (size_t i = 0; i < num_threads; ++i)
				tp.enqueue([this, i, &s] { eval_2(s, i); });
			tp.wait();
			for (size_t i = 0; i < num_threads; ++i)
			{
				energy_coulomb += partpot_coulomb[i];
				energy_lj += partpot_lj[i];
				s.virial += partvir[i];
			}
			eval_ewald_scd(s, energy_coulomb, kappa, volume);
			f = s.f - f;
			s.potential += energy_coulomb + energy_lj;
		}

		state<T, 3> f; // last force evaluation computed with this class
		// (the force is also accumulated inside s.f)
		T energy_coulomb, energy_lj; // energy for electrostatic and LJ potentials

		private:

			std::vector<vec3<T>> k; // wavevectors
			std::vector<std::complex<T>> structure; // structure factor
			std::vector<T> partpot_coulomb, partpot_lj, partvir;
			std::vector<T> factor;
			T cutoff = 0, kappa;
			std::size_t maxn = 6, num_threads, maxn1, maxn3;

			template <typename System>
			void eval_1(System& s, std::size_t idx)
			// 1st phase: compute force and energy in real space (eval_ewald_r)
			// and then compute structure factor and energy in reciprocal space.
			// `s` is the physical system (`coulomb_and_lj_periodic` constraints required).
			// `idx` is the index of the thread.
			// the potential and the virial are accumulated inside `partpot_*[idx]`
			// and `partvir[idx]` respectively.
			{
				using std::size_t;
				using std::ptrdiff_t;
				using std::sin;
				using std::cos;
				T volume = s.side*s.side*s.side;
				T cutoff2 = cutoff*cutoff;
				T uC = 0, uLJ = 0, vLJ = 0;
				// calculate real-space forces, energy and virial
				size_t block = (s.n-1)/num_threads + 1;
				for (size_t i = idx*block; (i < s.n) && (i < (idx+1)*block); ++i)
					for (size_t j = 0; j < s.n; ++j)
						if (i != j)
						{
							auto r = remainder(s.x[i] - s.x[j], s.side);
							T r2 = dot(r, r);
							if (r2 <= cutoff2)
								eval_ewald_r(s, uC, uLJ, vLJ, kappa, i, j, r, r2);
						}
				uC /= 2;
				uLJ /= 2;
				vLJ /= 2;
				// calculate structure factor, and reciprocal-space energy
				block = (maxn3-1)/num_threads + 1;
				for (size_t nijk = idx*block; (nijk < maxn3) && (nijk < (idx+1)*block); ++nijk)
				{
					if (nijk == 0)
					{
						k[nijk] = 0;
						structure[nijk] = 0;
						factor[nijk] = 0;
						continue;
					}
					ptrdiff_t n3 = nijk % maxn1;
					ptrdiff_t nij = nijk / maxn1;
					ptrdiff_t n2 = nij % maxn1;
					ptrdiff_t n1 = nij / maxn1;
					if (n1 > ptrdiff_t(maxn)) n1 -= ptrdiff_t(maxn1);
					if (n2 > ptrdiff_t(maxn)) n2 -= ptrdiff_t(maxn1);
					if (n3 > ptrdiff_t(maxn)) n3 -= ptrdiff_t(maxn1);
					k[nijk] = math::two_pi<T>*vec3<T>{n1, n2, n3}/s.side;
					T k2 = dot(k[nijk], k[nijk]);
					structure[nijk] = 0;
					for (size_t i = 0; i < s.n; ++i)
					{
						T kx = -dot(k[nijk], s.x[i]);
						structure[nijk] += s.z[i] * std::complex(cos(kx), sin(kx));
					}
					factor[nijk] = 2*math::two_pi<T>*math::fastexp(-k2/(4*kappa*kappa))/(k2*volume);
					T kspace_contrib = norm(structure[nijk])*factor[nijk]/2;
					uC += kspace_contrib;
				}
				partpot_coulomb[idx] = uC;
				partpot_lj[idx] = uLJ;
				partvir[idx] = vLJ + uC;
			}

			template <typename System>
			void eval_2(System& s, std::size_t idx)
			// 2nd phase: compute forces in reciprocal space.
			// `s` is the physical system.
			// `idx` is the index of the thread.
			{
				using std::size_t;
				using std::sin;
				using std::cos;
				size_t block = (s.n-1)/num_threads + 1;
				for (size_t i = idx*block; (i < s.n) && (i < (idx+1)*block); ++i)
					for (size_t nijk = 0; nijk < maxn3; ++nijk)
					{
						T kx = -dot(k[nijk], s.x[i]);
						s.f[i] += (factor[nijk] * s.z[i] * (structure[nijk] * std::complex(-sin(kx), cos(kx))).real()) * k[nijk];
					}
			}
	};

} // namespace physics

#endif // PHYSICS_EWALD_H
































