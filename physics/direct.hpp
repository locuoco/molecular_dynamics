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
	template <typename System, typename T, typename State>
	concept coulomb_and_lj = physical_system<System, T, State>
		&& requires(System& s, T t, std::size_t i)
	// A `coulomb_and_lj` system `s` is a `physical_system` so that the following
	// instructions are well-formed (compilable).
	// `n` is the number of atoms.
	// `z[i]` is the rescaled charge for particle `i`.
	// `lj_sqrteps[i]` is the square root of the Lennard-Jones energy parameter for particle `i`.
	// `lj_halfR[i]` is half the minimum LJ potential distance parameter for particle `i`.
	// `f[i]` is the force for particle `i`.
	// `potential` is the potential energy.
	// `virial` is the virial.
		{
			i < s.n;
			t += s.z[i];
			t += s.lj_sqrteps[i];
			t += s.lj_halfR[i];
			s.f[i] += s.x[i] - s.x[i];
			s.potential += t;
			s.virial += dot(s.x[i], s.x[i]);
		};

	template <bool Bidirectional, typename T, typename System>
	void eval_direct(System& s, T& u, T& v, std::size_t i, std::size_t j)
	// compute force between two particles `i` and `j`. It represents a single
	// iteration of a direct summation.
	// `s` is the physical system.
	// `u` is a variable which accumulates the potential energy.
	// `v` is a variable which accumulates the virial.
	// The `Bidirectional` template non-type parameter must be specified. If
	// it is `true`, the force for particle `j` will be calculated along with
	// the one for `i` (exploiting Newton's third law). Otherwise, the force
	// for particle `j` will not be calculated (this might be desiderable in
	// a parallel implementation).
	// Both Coulomb and Lennard-Jones potentials and forces are calculated
	// at the same time.
	{
		using std::sqrt;

		T Rij = s.lj_halfR[i] + s.lj_halfR[j];

		auto r = s.x[i] - s.x[j];
		T r2_ = 1/dot(r, r);
		T d_ = sqrt(r2_);
		T t = Rij * Rij * r2_;
		t = t * t * t;
		T epsijt = s.lj_sqrteps[i] * s.lj_sqrteps[j] * t;
		T coulomb = s.z[i] * s.z[j] * d_;
		T tot_virial = 12 * epsijt * (t - 1) + coulomb;
		auto fij = (tot_virial * r2_) * r;
		s.f[i] += fij;
		if constexpr (Bidirectional)
			s.f[j] -= fij;
		u += epsijt * (t - 2) + coulomb;
		v += tot_virial;
	}

	template <typename T, typename State>
	struct direct
	// Direct summation of non-bonded forces (electrostatic and Lennard-Jones)
	// assuming non-periodic boundaries.
	{
		static T cutoff_radius() noexcept
		// there is no cutoff radius as all contributes are being calculated
		{
			return 0;
		}

		template <coulomb_and_lj<T, State> System>
		void operator()(System& s)
		// single-threaded implementation
		// `s` is the physical system.
		{
			using std::size_t;
			for (size_t i = 0; i < s.n; ++i)
				for (size_t j = i+1; j < s.n; ++j)
					eval_direct<true>(s, s.potential, s.virial, i, j);
		}

		template <coulomb_and_lj<T, State> System>
		void operator()(System& s, utils::thread_pool& tp)
		// multi-threaded implementation
		// `s` is the physical system.
		// `tp` is a thread pool.
		{
			using std::size_t;
			auto num_threads = tp.size();
			// calling single-threaded implementation if the number of threads is less than 2
			if (num_threads <= 2)
				return operator()(s);
			partpot.resize(num_threads);
			partvir.resize(num_threads);
			auto eval_lambda = [this, num_threads, &s](size_t idx)
				{
					T u = 0, v = 0;
					size_t block = (s.n-1)/num_threads + 1;
					for (size_t i = idx*block; (i < s.n) && (i < (idx+1)*block); ++i)
						for (size_t j = 0; j < s.n; ++j)
							if (i != j)
								eval_direct<false>(s, u, v, i, j);
					partpot[idx] = u/2;
					partvir[idx] = v/2;
				};
			for (size_t i = 0; i < num_threads; ++i)
				tp.enqueue(eval_lambda, i);
			tp.wait();
			for (size_t i = 0; i < num_threads; ++i)
			{
				s.potential += partpot[i];
				s.virial += partvir[i];
			}
		}

		private:

			std::vector<T> partpot, partvir;
	};

} // namespace physics

#endif // PHYSICS_DIRECT_H
































