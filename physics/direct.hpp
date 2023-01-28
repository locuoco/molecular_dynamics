//  Direct summation for non-bonded forces
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

#ifndef PHYSICS_DIRECT_H
#define PHYSICS_DIRECT_H

#include <cmath> // sqrt, remainder

#include "tensor.hpp" // remainder
#include "physical_system.hpp"
#include "../utils/thread_pool.hpp"

namespace physics
{
	template <typename System>
	concept coulomb_and_lj = physical_system<System>
		&& requires(System& s, scalar_type_of<System> t, std::size_t i)
	// A `coulomb_and_lj` system `s` is a `physical_system` so that the following
	// instructions are well-formed (compilable).
	// `n` is the number of atoms.
	// `z[i]` is the charge for particle `i`, rescaled by the square root of the Coulomb constant.
	// `lj_sqrteps[i]` is the square root of the Lennard-Jones energy parameter for particle `i`.
	// `lj_halfR[i]` is half the minimum LJ potential distance parameter for particle `i`.
	// `f[i]` is the force for particle `i`.
	// `potential` is the potential energy.
	// `virial` is the virial (used to calculate the pressure).
	// The Lennard-Jones parameters for any interaction pair can be calculated on the fly
	// using the Lorentz-Berthelot mixing rules:
	//	lj_epsij = s.lj_sqrteps[i] * s.lj_sqrteps[j];
	//	lj_Rij = s.lj_halfR[i] + s.lj_halfR[j];
		{
			i < s.n;
			t += s.z[i];
			t += s.lj_sqrteps[i];
			t += s.lj_halfR[i];
			s.f[i] += s.x[i] - s.x[i];
			s.potential += t;
			s.virial += dot(s.x[i], s.x[i]);
		};

	template <typename System>
	concept coulomb_and_lj_periodic = coulomb_and_lj<System>
		&& requires(System& s, scalar_type_of<System> t, std::size_t i)
	// A `coulomb_and_lj_periodic` system `s` is a `coulomb_and_lj` system so that the following
	// instructions are well-formed (compilable).
	// `side` must be the side of the simulation box.
	// `sumdisp62` must be the sum square of all dispersion coefficients for each atom.
	// `Z` must be the sum of all charges (rescaled by the square root of kC).
	// `Z2` must be the sum square of all charges.
		{
			s.f[i] += remainder(s.x[i] - s.x[i], s.side);
			t *= s.sumdisp62;
			t += s.Z;
			t += s.Z2;
		};

	template <bool Bidirectional, coulomb_and_lj System>
	void eval_direct(
		System&                                                            s,
		scalar_type_of<System>&                                            uC,
		scalar_type_of<System>&                                            uLJ,
		scalar_type_of<System>&                                            vLJ,
		std::size_t                                                        i,
		std::size_t                                                        j,
		std::add_const_t<std::remove_reference_t<decltype(System::x[i])>>& r,
		scalar_type_of<System>                                             r2
	)
	// compute force between two particles `i` and `j`. It represents a single
	// iteration of a direct summation.
	// `s` is the physical system (`coulomb_and_lj` constraints required).
	// `uC` is a variable which accumulates the electrostatic potential energy.
	// `uLJ` is a variable which accumulates the Lennard-Jones potential energy.
	// `vLJ` is a variable which accumulates the Lennard-Jones virial.
	// `r` is the displacement vector between the two particles, given by:
	//	r = remainder(s.x[i] - s.x[j]).
	// `r2` is the squared distance between the two particles.
	// The `Bidirectional` template non-type parameter must be specified. If
	// it is `true`, the force for particle `j` will be calculated along with
	// the one for `i` (exploiting Newton's third law). Otherwise, the force
	// for particle `j` will not be calculated (this might be desiderable in
	// a parallel implementation).
	// Both Coulomb and Lennard-Jones potentials and forces are calculated
	// at the same time.
	{
		using std::sqrt;

		scalar_type_of<System> Rij = s.lj_halfR[i] + s.lj_halfR[j];

		scalar_type_of<System> r2_ = 1/r2;
		scalar_type_of<System> d_ = sqrt(r2_);
		scalar_type_of<System> t = Rij * Rij * r2_;
		t = t * t * t;
		scalar_type_of<System> epsijt = s.lj_sqrteps[i] * s.lj_sqrteps[j] * t;
		scalar_type_of<System> coulomb = s.z[i] * s.z[j] * d_;
		scalar_type_of<System> lennard_jones = 12 * epsijt * (t - 1);
		auto fij = ((lennard_jones + coulomb) * r2_) * r;
		s.f[i] += fij;
		if constexpr (Bidirectional)
			s.f[j] -= fij;
		uC += coulomb;
		uLJ += epsijt * (t - 2);
		vLJ += lennard_jones;
	}

	template <bool Bidirectional, coulomb_and_lj System>
	void eval_direct(
		System&                 s,
		scalar_type_of<System>& uC,
		scalar_type_of<System>& uLJ,
		scalar_type_of<System>& vLJ,
		std::size_t             i,
		std::size_t             j
	)
	// compute force between two particles `i` and `j`. It represents a single
	// iteration of a direct summation.
	// `s` is the physical system (`coulomb_and_lj` constraints required).
	// `uC` is a variable which accumulates the electrostatic potential energy.
	// `uLJ` is a variable which accumulates the Lennard-Jones potential energy.
	// `vLJ` is a variable which accumulates the Lennard-Jones virial.
	// The `Bidirectional` template non-type parameter must be specified. If
	// it is `true`, the force for particle `j` will be calculated along with
	// the one for `i` (exploiting Newton's third law). Otherwise, the force
	// for particle `j` will not be calculated (this might be desiderable in
	// a parallel implementation).
	// This overload calculates the displacement `r` and its square `r2` between the two
	// particles and then calls eval_direct<Bidirectional>(s, uC, uLJ, v, i, j, r, r2).
	{
		auto r = s.x[i] - s.x[j];
		scalar_type_of<System> r2 = dot(r, r);
		eval_direct<Bidirectional>(s, uC, uLJ, vLJ, i, j, r, r2);
	}

	template <typename T, typename State>
	struct direct
	// Direct summation of non-bonded forces (electrostatic and Lennard-Jones)
	// assuming non-periodic boundaries.
	{
		static constexpr T cutoff_radius() noexcept
		// There is no cutoff radius as all contributes are being calculated.
		// Return 0.
		{
			return 0;
		}

		template <coulomb_and_lj System>
		void operator()(System& s)
		// single-threaded implementation
		// `s` is the physical system.
		{
			using std::size_t;
			energy_coulomb = energy_lj = 0;
			for (size_t i = 0; i < s.n; ++i)
				for (size_t j = i+1; j < s.n; ++j)
					eval_direct<true>(s, energy_coulomb, energy_lj, s.virial, i, j);
			s.potential += energy_coulomb + energy_lj;
			s.virial += energy_coulomb;
		}

		template <coulomb_and_lj System>
		void operator()(System& s, utils::thread_pool& tp)
		// multi-threaded implementation
		// `s` is the physical system (`coulomb_and_lj` constraints required).
		// `tp` is a thread pool.
		{
			using std::size_t;
			auto num_threads = tp.size();
			// calling single-threaded implementation if the number of threads is less than 2
			if (num_threads <= 2)
				return (*this)(s);
			partpot_coulomb.resize(num_threads);
			partpot_lj.resize(num_threads);
			partvir_lj.resize(num_threads);
			auto eval_lambda = [this, num_threads, &s](size_t idx)
				// lambda function calculating forces, potentials and virials
				// for thread `idx`.
				{
					T uC = 0, uLJ = 0, vLJ = 0;
					size_t block = (s.n-1)/num_threads + 1;
					for (size_t i = idx*block; (i < s.n) && (i < (idx+1)*block); ++i)
						for (size_t j = 0; j < s.n; ++j)
							if (i != j)
								eval_direct<false>(s, uC, uLJ, vLJ, i, j);
					partpot_coulomb[idx] = uC/2;
					partpot_lj[idx] = uLJ/2;
					partvir_lj[idx] = vLJ/2;
				};
			for (size_t i = 0; i < num_threads; ++i)
				tp.enqueue(eval_lambda, i);
			tp.wait();
			energy_coulomb = energy_lj = 0;
			for (size_t i = 0; i < num_threads; ++i)
			{
				energy_coulomb += partpot_coulomb[i];
				energy_lj += partpot_lj[i];
				s.virial += partvir_lj[i];
			}
			s.potential += energy_coulomb + energy_lj;
			s.virial += energy_coulomb;
		}

		T energy_coulomb, energy_lj; // energy for electrostatic and LJ potentials

		private:

			std::vector<T> partpot_coulomb, partpot_lj, partvir_lj;
	};

	template <typename T, typename State>
	struct direct_periodic
	// Direct summation of non-bonded forces (electrostatic and Lennard-Jones)
	// assuming periodic boundaries. Only a finite amount of periodic images will be considered.
	// The summation will be performed along cubic shells. Note that energy may be very
	// slow to converge or not converge at all to the correct value.
	{
		static constexpr T cutoff_radius() noexcept
		// Cutoff radius undefined, since summation is performed expanding from cubes.
		// Return 0.
		{
			return 0;
		}

		std::size_t max_n_shell() noexcept
		// return the maximum distance in terms of cubic shells considered
		{
			return maxn;
		}
		void max_n_shell(std::size_t maxn) noexcept
		// set the maximum distance in terms of cubic shells considered
		{
			this->maxn = maxn;
		}

		template <coulomb_and_lj_periodic System>
		void operator()(System& s, utils::thread_pool& tp)
		// multi-threaded implementation
		// `s` is the physical system (`coulomb_and_lj_periodic` constraints required).
		// `tp` is a thread pool.
		{
			using std::size_t;
			num_threads = tp.size();
			partpot_coulomb.resize(num_threads);
			partpot_lj.resize(num_threads);
			partvir_lj.resize(num_threads);
			for (size_t i = 0; i < num_threads; ++i)
				tp.enqueue([this, i, &s] { eval(s, i); } );
			tp.wait();
			energy_coulomb = energy_lj = 0;
			for (size_t i = 0; i < num_threads; ++i)
			{
				energy_coulomb += partpot_coulomb[i];
				energy_lj += partpot_lj[i];
				s.virial += partvir_lj[i];
			}
			s.potential += energy_coulomb + energy_lj;
			s.virial += energy_coulomb;
		}

		T energy_coulomb, energy_lj; // energy for electrostatic and LJ potentials

		private:

			std::vector<T> partpot_coulomb, partpot_lj, partvir_lj;
			std::size_t maxn = 10, num_threads;

			template <typename System>
			void eval(System& s, std::size_t idx)
			// Compute forces.
			// `s` is the physical system.
			// `idx` is the index of the thread.
			{
				using std::size_t;
				T uC = 0, uLJ = 0, vLJ = 0;
				size_t block = (s.n-1)/num_threads + 1;
				for (size_t i = idx*block; (i < s.n) && (i < (idx+1)*block); ++i)
					for (size_t j = 0; j < s.n; ++j)
						if (i != j)
						{
							auto r = remainder(s.x[i] - s.x[j], s.side);
							T r2 = dot(r, r);
							eval_direct<false>(s, uC, uLJ, vLJ, i, j, r, r2);
						}
				for (size_t i = idx*block; (i < s.n) && (i < (idx+1)*block); ++i)
					for (size_t j = 0; j < s.n; ++j)
					{
						auto r0 = remainder(s.x[i] - s.x[j], s.side);
						for (int k1 = -int(maxn); k1 <= int(maxn); ++k1)
							for (int k2 = -int(maxn); k2 <= int(maxn); ++k2)
								for (int k3 = -int(maxn); k3 <= int(maxn); ++k3)
								{
									if (k1 == 0 && k2 == 0 && k3 == 0)
										continue;
									decltype(-s.x[i]) k(k1, k2, k3);
									auto r = r0 + k * s.side;
									T r2 = dot(r, r);
									eval_direct<false>(s, uC, uLJ, vLJ, i, j, r, r2);
								}
					}
				partpot_coulomb[idx] = uC/2;
				partpot_lj[idx] = uLJ/2;
				partvir_lj[idx] = vLJ/2;
			}
	};

} // namespace physics

#endif // PHYSICS_DIRECT_H
































