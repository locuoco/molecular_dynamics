//  Energy minimizers
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

#ifndef PHYSICS_ENERGY_MINIMIZER_H
#define PHYSICS_ENERGY_MINIMIZER_H

#include "physical_system.hpp"
#include "integrator/integrator.hpp" // leapfrog
#include "tensor.hpp" // rms

namespace physics
{
	template <typename System>
	struct energy_minimizer_base
	{
		virtual ~energy_minimizer_base() = default;

		virtual void step(System& s) = 0;

		void minimize(System& s, std::size_t n_steps)
		// Minimize for `n_steps` steps.
		// `s` is the system for which the energy is minimized.
		// Same as:
		//	for (std::size_t i = 0; i < n_steps; ++i) step(s);
		{
			for (std::size_t i = 0; i < n_steps; ++i)
				step(s);
		}

		template <typename S = System>
		requires requires(S& s, scalar_type_of<S> t)
			{
				t = rms(s.f);
			}
		bool minimize(System& s, scalar_type_of<System> f_rms, std::size_t max_steps)
		// Minimize until convergence (RMS of forces smaller than `f_rms`).
		// `s` is the system for which the energy is minimized.
		// `max_steps` is the maximum number of steps.
		// return whether the minimization converged.
		{
			bool converged = true;
			for (std::size_t i = 0; (converged = rms(s.f) >= f_rms) && i < max_steps; ++i)
				step(s);
			return converged;
		}
	};

	// All energy minimizers derive from `energy_minimizer_base` and will have a method
	// called `step` with the following signature:

	//		void step(System& s);

	// where `System` is the template argument of the integrator class, corresponding to
	// the class of the system to integrate.

	// concepts associated to these template classes (can be used as constraints in template
	// type arguments)
	template <typename Minimizer, typename System>
	concept energy_minimizer = std::is_base_of_v<energy_minimizer_base<System>, Minimizer>;

	template <having_coordinates System, template <typename> typename IntegT = leapfrog>
	requires requires(System& s, scalar_type_of<System> t)
		{
			t = dot(s.f, s.p);
			t = norm(s.f);
			t = norm(s.p);
		}
	struct fire : energy_minimizer_base<System>
	// FIRE (Fast Inertial Relaxation Engine) minimizer
	// see E. Bitzek, P. Koskinen, F. Gahler, M. Moseler, P. Gumbsch, "Structural Relaxation Made Simple"
	{
		fire(scalar_type init_dt, const IntegT<System>& integ = IntegT<System>())
		// constructor:
		// `init_dt` is the initial time-step.
		// `integ` is the integrator to be used.
			: dt_start(init_dt), integ(integ)
		{
			reset();
		}

		void step(System& s) override
		{
			integ.step(s, dt);

			scalar_type P = dot(s.f, s.p);
			s.p = (1 - alpha) * s.p + alpha*norm(s.p)/norm(s.f) * s.f;
			if (P > 0)
			{
				++n_positive;
				if (n_positive > n_min)
				{
					dt = std::min(f_inc*dt, max_dt);
					alpha *= f_alpha;
				}
			}
			else
			{
				n_positive = 0;
				dt *= f_dec;
				s.p = 0;
				alpha = alpha_start;
			}
		}

		void reset()
		// reset minimizer (for a new minimization)
		{
			dt = dt_start;
			max_dt = 10*dt_start;
			alpha = alpha_start;
			n_positive = 0;
		}

		private:

			using scalar_type = scalar_type_of<System>;

			IntegT<System> integ;
			std::size_t n_positive;
			scalar_type dt_start, dt, max_dt, alpha;

			static constexpr std::size_t n_min = 5;
			static constexpr scalar_type f_inc = scalar_type(1.1), f_dec = scalar_type(0.5);
			static constexpr scalar_type alpha_start = scalar_type(0.5), f_alpha = scalar_type(0.99);
	};

} // namespace physics

#endif // PHYSICS_ENERGY_MINIMIZER_H






























