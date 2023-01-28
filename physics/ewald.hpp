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

#include <iostream> // clog
#include <type_traits> // add_const_t, remove_reference_t
#include <cmath> // sqrt, sin, cos, remainder
#include <numbers> // numbers::inv_sqrtpi_v, numbers::pi_v
#include <complex>
#include <limits> // infinity

#include "tensor.hpp" // remainder
#include "physical_system.hpp"
#include "direct.hpp" // coulomb_and_lj
#include "../math/helper.hpp" // fasterfc, fastexp
#include "../utils/thread_pool.hpp"

namespace physics
{
	template <bool Fast = true, coulomb_and_lj_periodic System>
	void eval_ewald_r(
		System&                                                            s,
		scalar_type_of<System>&                                            uC,
		scalar_type_of<System>&                                            uLJ,
		scalar_type_of<System>&                                            vLJ,
		scalar_type_of<System>                                             kappa,
		std::size_t                                                        i,
		std::size_t                                                        j,
		std::add_const_t<std::remove_reference_t<decltype(System::x[i])>>& r,
		scalar_type_of<System>                                             r2
	)
	// Calculate the real space part of Ewald summation between two particles `i` and `j`.
	// `s` is the physical system (`coulomb_and_lj_periodic` constraints required).
	// `uC` is a variable which accumulates the electrostatic potential energy.
	// `uLJ` is a variable which accumulates the Lennard-Jones potential energy.
	// `vLJ` is a variable which accumulates the Lennard-Jones virial.
	// `kappa` is the Ewald parameter.
	// `r` is the displacement vector between the two particles, given by:
	//	r = remainder(s.x[i] - s.x[j]).
	// `r2` is the squared distance between the two particles.
	// If the template argument `Fast` is set to false, the standard functions
	// `std::exp` and `std::erfc` are used instead of `math::fastexp` and
	// `math::fasterfc` (defined in math/helper.hpp).
	// Both Coulomb and Lennard-Jones potentials and forces are calculated at
	// the same time, but Ewald summation is performed only for Coulomb potential.
	{
		using std::sqrt;
		using scalar_type = scalar_type_of<System>;

		scalar_type Rij = s.lj_halfR[i] + s.lj_halfR[j];

		scalar_type r2_ = 1/r2, d = sqrt(r2), d_ = 1/d;
		scalar_type kd = kappa * d;
		scalar_type t = Rij * Rij * r2_;
		t = t * t * t;
		scalar_type epsijt = s.lj_sqrteps[i] * s.lj_sqrteps[j] * t;
		scalar_type zij_d = s.z[i] * s.z[j] * d_;
		scalar_type coulomb;
		if constexpr (Fast)
			coulomb = zij_d * math::fasterfc(kd);
		else
			coulomb = zij_d * std::erfc(kd);
		scalar_type lennard_jones = 12 * epsijt * (t - 1);
		decltype(-r) fij;
		if constexpr (Fast)
			fij = ((lennard_jones + coulomb + zij_d * (2 * std::numbers::inv_sqrtpi_v<scalar_type>) * kd * math::fastexp(-kd*kd)) * r2_) * r;
		else
			fij = ((lennard_jones + coulomb + zij_d * (2 * std::numbers::inv_sqrtpi_v<scalar_type>) * kd * std::exp(-kd*kd)) * r2_) * r;
		s.f[i] += fij;
		uC += coulomb;
		uLJ += epsijt * (t - 2);
		vLJ += lennard_jones;
	}

	template <bool Fast = true, coulomb_and_lj_periodic System>
	void eval_ewald_r(
		System&                 s,
		scalar_type_of<System>& uC,
		scalar_type_of<System>& uLJ,
		scalar_type_of<System>& vLJ,
		scalar_type_of<System>  kappa,
		std::size_t             i,
		std::size_t             j
	)
	// Calculate the real space part of Ewald summation between two particles `i` and `j`.
	// `s` is the physical system (`coulomb_and_lj` constraints required).
	// `uC` is a variable which accumulates the electrostatic potential energy.
	// `uLJ` is a variable which accumulates the Lennard-Jones potential energy.
	// `vLJ` is a variable which accumulates the Lennard-Jones virial.
	// `kappa` is the Ewald parameter.
	// If the template argument `Fast` is set to false, the standard functions
	// `std::exp` and `std::erfc` are used instead of `math::fastexp` and
	// `math::fasterfc` (defined in math/helper.hpp).
	// This overload calculates the displacement `r` and its square `r2` between the two
	// particles and then calls eval_ewald_r(s, uC, uLJ, v, kappa, i, j, r, r2).
	{
		auto r = remainder(s.x[i] - s.x[j], s.side);
		eval_ewald_r<Fast>(s, uC, uLJ, vLJ, kappa, i, j, r, dot(r, r));
	}

	template <coulomb_and_lj_periodic System>
	void eval_ewald_scd(
		System&                 s,
		scalar_type_of<System>& uC,
		scalar_type_of<System>  kappa,
		scalar_type_of<System>  volume,
		scalar_type_of<System>  dielectric
	)
	// Calculate self-energy, charged system and dipole corrections.
	// `s` is the physical system (`coulomb_and_lj` constraints required).
	// `uC` is a variable which accumulates the electrostatic potential energy.
	// `kappa` is the Ewald parameter.
	// `volume` is the volume of the system.
	// `dielectric` is the relative dielectric constant of the space surrounding the system
	// (1 for vacuum).
	{
		using std::size_t;
		using scalar_type = scalar_type_of<System>;

		vec3<scalar_type> xsum = 0;
		for (size_t i = 0; i < s.n; ++i)
			xsum += s.z[i] * s.x[i];
		scalar_type fact = math::two_pi<scalar_type>/((1+2*dielectric)*volume);
		scalar_type u = fact*dot(xsum, xsum); // dipole correction
		u += -kappa*std::numbers::inv_sqrtpi_v<scalar_type> * s.Z2; // self-energy correction
		u += -std::numbers::pi_v<scalar_type>/(2*kappa*kappa*volume)*s.Z*s.Z; // charged system correction
		uC += u;
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

		T dielectric() const noexcept
		// Return the external relative dielectric constant (infinity by default)
		{
			return dielec;
		}
		void dielectric(T dielec) noexcept
		// Set the external relative dielectric constant (infinity by default)
		{
			this->dielec = dielec;
		}

		void max_n(std::size_t maxn) noexcept
		// Set the maximum reciprocal number, related to the maximum reciprocal vector:
		// 	kmax = 2 pi max_n / side
		{
			this->maxn = maxn;
			update = true;
		}
		std::size_t max_n() const noexcept
		// Return the maximum reciprocal number, related to the maximum reciprocal vector:
		// 	kmax = 2 pi max_n / side
		{
			return maxn;
		}

		void ewald_par(T kappa)
		// Set the Ewald parameter manually. If this method is not called, the Ewald
		// parameter is optimized automatically.
		// Throw a `std::invalid_argument` if `kappa` is not strictly positive.
		{
			if (kappa <= 0)
				throw std::invalid_argument("Error: Ewald parameter must be positive");
			this->kappa = kappa;
			manual = true;
		}
		T ewald_par() const noexcept
		// Return the Ewald parameter.
		// If `operator()` is not called yet, return 0. The actual Ewald parameter will
		// be calculated only after `operator()` is called.
		{
			return kappa;
		}

		void precise(bool flag) noexcept
		// If `flag` is set to true, the calculations are carried on with the highest
		// precision. By default it is set to false, enabling faster calculations.
		{
			fast = !flag;
		}

		void update_ewald() noexcept
		// set update flag to true (calculates optimal Ewald parameter again)
		{
			update = true;
			manual = false;
		}

		template <coulomb_and_lj_periodic System>
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
			// update Ewald paremeter every once in a while for constant-P simulations
			if (update_counter % update_max_count == (update_max_count-1))
			{
				update = true;
				++update_counter;
			}
			if (update)
				optimize_ewald_par(s);
			else if (prev_side != s.side)
			{
				// rescale Ewald parameter
				kappa *= prev_side / s.side;
				++update_counter;
			}
			num_threads = tp.size();
			partpot_r.resize(num_threads);
			partpot_k.resize(num_threads);
			partpot_lj.resize(num_threads);
			partvir_lj.resize(num_threads);

			k.resize(maxn3);
			factor.resize(maxn3);
			structure.resize(maxn3);

			for (size_t i = 0; i < num_threads; ++i)
				tp.enqueue([this, i, &s] { eval_1(s, i); });
			tp.wait();
			for (size_t i = 0; i < num_threads; ++i)
				tp.enqueue([this, i, &s] { eval_2(s, i); });
			tp.wait();
			energy_r = energy_k = energy_scd = energy_lj = 0;
			for (size_t i = 0; i < num_threads; ++i)
			{
				energy_r += partpot_r[i];
				energy_k += partpot_k[i];
				energy_lj += partpot_lj[i];
				s.virial += partvir_lj[i];
			}
			eval_ewald_scd(s, energy_scd, kappa, volume, dielec);
			energy_coulomb = energy_r + energy_k + energy_scd;
			s.potential += energy_coulomb + energy_lj;
			s.virial += energy_coulomb;
			update = false;
			prev_side = s.side;
		}

		T estimated_error, estimated_error_coulomb, estimated_error_lj; // estimated force RMS error
		T energy_coulomb, energy_lj; // energy for electrostatic and LJ potentials
		T energy_r, energy_k, energy_scd; // real, reciprocal and corrections parts of electrostatic potential

		private:

			std::vector<vec3<T>> k; // wavevectors
			std::vector<std::complex<T>> structure; // structure factor
			std::vector<T> partpot_r, partpot_k, partpot_lj, partvir_lj;
			std::vector<T> factor;
			T cutoff = 0, kappa = 0, dielec = std::numeric_limits<T>::infinity(), prev_side;
			std::size_t maxn = 6, num_threads, maxn1, maxn3;
			unsigned update_counter, update_max_count = 1000;
			bool update = true, manual = false, fast = true, verbose = false;

			template <typename System>
			void optimize_ewald_par(System& s)
			// Calculate (approximate) optimal Ewald parameter with Newton's method.
			// See Kolafa et al. (1992)
			// `s` is the physical system.
			{
				using std::abs;
				using std::sqrt;
				if (!manual)
					kappa = 7/s.side;
				T cutoff2 = cutoff*cutoff;
				T next_kappa = kappa, errF;
				T corr = T(1)/1000; // empirical correcting for overestimation
				unsigned counter = 0, max_counter = 20;
				do
				{
					kappa = next_kappa;
					T k2 = kappa*kappa;
					// real-space force error
					T errFr = math::fastexp(-2*cutoff2*k2)/cutoff;
					T errFr1 = -4*kappa*cutoff2*errFr; // 1st derivative
					T errFr2 = 4*cutoff2*errFr*(4*k2*cutoff2 - 1); // 2nd derivative
					// reciprocal-space force error
					// note that this estimation assumes spherical cutoff, which is likely to
					// overestimate the error of a cubic cutoff, so we multiply by `corr`
					T factor = (k2*s.side)/(2*std::numbers::pi_v<T>*std::numbers::pi_v<T>*maxn);
					T expo = maxn/(s.side*factor);
					T errFk = math::fastexp(-expo) * factor * corr;
					T factor2 = 2*(1 + expo) / kappa;
					T errFk1 = errFk * factor2; // 1st derivative
					T errFk2 = errFk * (factor2*factor2 - 2*(1+3*expo)/k2); // 2nd derivative
					// total force error
					errF = errFr + errFk;
					T errF1 = errFr1 + errFk1; // 1st derivative
					T errF2 = errFr2 + errFk2; // 2nd derivative

					next_kappa -= errF1/errF2;
					++counter;
				}
				while (abs((next_kappa-kappa) / kappa) > 1e-3 && counter < max_counter && !manual);

				T volume = s.side*s.side*s.side;
				estimated_error_coulomb = 2 * s.Z2 * sqrt(errF / (s.n * volume));
				estimated_error_lj = 12 / (cutoff2 * cutoff2 * cutoff) * sqrt(std::numbers::pi_v<T> * s.sumdisp62 / (11 * cutoff * s.n * volume));
				estimated_error = sqrt(estimated_error_coulomb*estimated_error_coulomb + estimated_error_lj*estimated_error_lj);

				if (verbose)
				{
					std::clog << "===== EWALD SUMMATION LOG =====\n";
					std::clog << "Ewald parameter (A^-1): " << kappa << '\n';
					std::clog << "estimated electrostatic force RMS error (kcal/(mol A)): " << estimated_error_coulomb << '\n';
					std::clog << "estimated total force RMS error (kcal/(mol A)): " << estimated_error << '\n';
					std::clog << "r_max (A): " << cutoff << "\n\n";
				}
				if (counter == max_counter)
					std::clog << "Warning: Ewald parameter optimization may have not converged!\n";
				if (estimated_error < 1e-5 && fast)
					std::clog << "Warning: estimated accuracy may not be reliable if `precise` is not set to true.\n";
			}

			template <typename System>
			void eval_1(System& s, std::size_t idx)
			// 1st phase: compute force and energy in real space (eval_ewald_r)
			// and then compute structure factor and energy in reciprocal space.
			// `s` is the physical system (`coulomb_and_lj_periodic` constraints required).
			// `idx` is the index of the thread.
			// the potential and the virial are accumulated inside `partpot_*[idx]`
			// and `partvir_*[idx]` respectively.
			{
				using std::size_t;
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
							{
								if (fast)
									eval_ewald_r<true>(s, uC, uLJ, vLJ, kappa, i, j, r, r2);
								else
									eval_ewald_r<false>(s, uC, uLJ, vLJ, kappa, i, j, r, r2);
							}
						}
				partpot_r[idx] = uC/2;
				partpot_lj[idx] = uLJ/2;
				partvir_lj[idx] = vLJ/2;
				uC = 0;
				// calculate structure factor, and reciprocal-space energy
				block = (maxn3-1)/num_threads + 1;
				for (size_t nijk = idx*block; (nijk < maxn3) && (nijk < (idx+1)*block); ++nijk)
				{
					structure[nijk] = 0;
					if (nijk == 0)
					{
						k[nijk] = 0;
						factor[nijk] = 0;
						continue;
					}
					vec3<int> nn;
					nn[2] = nijk % maxn1;
					int nij = nijk / maxn1;
					nn[1] = nij % maxn1;
					nn[0] = nij / maxn1;
					nn -= int(maxn1)*(nn/int(maxn+1));

					k[nijk] = math::two_pi<T>*vec3<T>(nn)/s.side;
					T k2 = dot(k[nijk], k[nijk]);
					for (size_t i = 0; i < s.n; ++i)
					{
						T kx = -dot(k[nijk], s.x[i]);
						structure[nijk] += s.z[i] * std::complex(cos(kx), sin(kx));
					}
					if (fast)
						factor[nijk] = 2*math::two_pi<T>*math::fastexp(-k2/(4*kappa*kappa))/(k2*volume);
					else
						factor[nijk] = 2*math::two_pi<T>*std::exp(-k2/(4*kappa*kappa))/(k2*volume);
					uC += norm(structure[nijk])*factor[nijk];
				}
				partpot_k[idx] = uC/2;
			}

			template <typename System>
			void eval_2(System& s, std::size_t idx)
			// 2nd phase: compute reciprocal-space contribution to forces.
			// `s` is the physical system.
			// `idx` is the index of the thread.
			{
				using std::size_t;
				using std::sin;
				using std::cos;
				size_t block = (s.n-1)/num_threads + 1;
				for (size_t i = idx*block; (i < s.n) && (i < (idx+1)*block); ++i)
				{
					decltype(-s.f[i]) ei;
					ei = 0;
					for (size_t nijk = 0; nijk < maxn3; ++nijk)
					{
						T kx = dot(k[nijk], s.x[i]);
						ei += (factor[nijk] * (structure[nijk] * std::complex(sin(kx), -cos(kx))).real()) * k[nijk];
					}
					s.f[i] += ei * s.z[i];
				}
			}
	};

} // namespace physics

#endif // PHYSICS_EWALD_H































