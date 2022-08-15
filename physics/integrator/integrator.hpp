//  Symplectic and stochastic integrators
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

#ifndef PHYSICS_INTEGRATOR_INTEGRATOR_H
#define PHYSICS_INTEGRATOR_INTEGRATOR_H

#include <cmath> // exp, sqrt

#include "integrator_base.hpp"

namespace physics
{
	template <typename T, typename State>
	struct symplectic_euler : symplectic_integrator_base<T, State>
	// Symplectic Euler method (1st order, 1 stage)
	// It is equivalent to leapfrog after correcting the initial conditions
	{
		template <having_coordinates<T, State> System>
		void step(System& s, T dt, bool = false) const
		{
			s.p += s.force() * dt;
			s.x += s.vel() * dt;

			s.t += dt;
		}
	};

	template <typename T, typename State>
	struct leapfrog : symplectic_integrator_base<T, State>
	// Leapfrog method (2nd order, 1 stage)
	{
		template <having_coordinates<T, State> System>
		void step(System& s, T dt, bool first_step = false) const
		{
			s.p += s.force(first_step) * (dt/2);
			s.x += s.vel() * dt;
			s.p += s.force() * (dt/2);

			s.t += dt;
		}
	};

	template <typename T, typename State>
	struct multi_timestep_leapfrog : symplectic_integrator_base<T, State>
	// Multi-timestep leapfrog method (2nd order, 1 long stage, n_short short stages)
	// If the force is divided into a high frequency and a low frequency part,
	// and the low frequency is much more expensive to compute, then this method
	// may be useful to reduce compute time or to increase accuracy.
	// In the microcanonical (NVE) ensemble this method can better conserve the
	// total energy without reducing the timestep
	{
		multi_timestep_leapfrog(std::size_t n_short = 10) : n_short(n_short)
		{}

		template <having_coordinates<T, State> System>
		requires requires(System& s, T dt)
			{
				s.p += s.force_long() * dt;
				s.p += s.force_short() * dt;
				s.p += s.force_long(false) * dt;
			}
		void step(System& s, T dt, bool first_step = false) const
		{
			T deltat = dt/n_short;
			s.p += s.force_long(first_step) * (dt/2);
			s.p += s.force_short() * (deltat/2);
			for (std::size_t i = 0; i < n_short-1; ++i)
			{
				s.x += s.vel() * deltat;
				s.p += s.force_short() * deltat;
			}
			s.x += s.vel() * deltat;
			s.p += s.force_short() * (deltat/2);
			s.p += s.force_long() * (dt/2);

			s.t += dt;
		}

		private:

			std::size_t n_short;
	};

	template <typename System, typename T, typename State>
	concept having_coordinates_damped = having_coordinates<System, T, State>
		&& requires(System& s, T dt, std::size_t i)
	// A system `having_coordinates_damped` requires to be a `having_coordinates`,
	// and has `gamma` (the damping coefficient, which must be dereferenceable) defined.
	// The following operations are required to be defined.
		{
			{s.gamma[i]} -> std::convertible_to<T>;
			s.p[i] = s.p[i] * dt;
			i < s.n;
		};

	template <typename T, typename State>
	struct damped_leapfrog : integrator_base<T, State>
	// damped integrator (2nd order?, 1 stage)
	// dp = f dt - gamma p dt
	{
		template <having_coordinates_damped<T, State> System>
		void step(System& s, T dt, bool first_step = false) const
		{
			using std::exp;

			s.p += s.force(first_step) * (dt/2);
			s.x += s.vel() * (dt/2);
			for (std::size_t i = 0; i < s.n; ++i)
				s.p[i] *= exp(-s.gamma[i] * dt);
			s.x += s.vel() * (dt/2);
			s.p += s.force() * (dt/2);

			s.t += dt;
		}
	};

	template <typename System, typename T, typename State>
	concept having_coordinates_stochastic = having_coordinates_damped<System, T, State>
		&& requires(System& s, T dt, std::size_t i)
	// A system `having_coordinates_stochastic` requires to be a `having_coordinates_damped`,
	// and has `D` (the diffusion coefficient, which is a scalar), `noise` and `rand` defined. `rand`
	// is a method that must generate random numbers for the `noise` vector (which must be
	// deferenceable and its elements should be convertible to s.p[i]).
	// The following operations are required to be defined.
		{
			{s.D} -> std::convertible_to<T>;
			s.p[i] = s.p[i] * dt + s.noise[i] * dt;
			s.rand();
		};

	template <typename T, typename State>
	struct stochastic_leapfrog : stochastic_integrator_base<T, State>
	// stochastic leapfrog (2nd order?, 1 stage)
	// dp = f dt - gamma p dt + sigma dw
	// with sigma dw = sqrt(2 gamma) * noise
	// For Langevin's equation the following equations must be satisfied:
	// 		gamma = k T / m D [related to diffusion coefficient],
	// 		sigma = sqrt(2 gamma m k T) = k T sqrt(2 / D) [fluctuation-dissipation theorem],
	//		dw : Wiener process differential
	// H = H0(x, p) + H1(p) xi(t) <- stochastic hamiltonian (xi : gaussian process)
/*
	ALLEN, TILDESLEY
	COMPUTER SIMULATION OF LIQUIDS
	2017
	P. 383-389
*/
	{
		template <having_coordinates_stochastic<T, State> System>
		void step(System& s, T dt, bool first_step = false) const
		{
			using std::exp;
			using std::sqrt;

			s.p += s.force(first_step) * (dt/2);
			s.x += s.vel() * (dt/2);
			s.rand();
			if (s.D == 0) // gamma -> inf
				s.p = s.noise; // s.noise = sqrt(m k T) * gaussian noise
			else
				for (std::size_t i = 0; i < s.n; ++i)
					s.p[i] = exp(-s.gamma[i] * dt) * s.p[i] + sqrt(1 - exp(-2 * s.gamma[i] * dt)) * s.noise[i];
				// p = exp(-gamma dt) p + sqrt((1 - exp(-2 gamma dt)) m k T) g
				// expanding the exp (gamma small): p += -gamma p dt + sqrt(2 gamma m k T dt) g
				// g: gaussian noise with mean = 0, std.dev = 1
			s.x += s.vel() * (dt/2);
			s.p += s.force() * (dt/2);

			s.t += dt;
		}
	};

	template <typename T, typename State>
	struct isokinetic_leapfrog : integrator_base<T, State>
	// (gaussian) isokinetic integrator (2nd order, 1 stage)
	// it assumes quadratic kinetic energy
	// not symplectic, but still time-reversible
	// it conserves the kinetic energy rather than the total hamiltonian
	// at equilibrium, its configurations sample a canonical ensemble rather than a microcanonical one
	// unfortunately, momenta do not follow the MB distribution
	// dp = f dt - xi p dt
	// xi = sum(p * f / m) / sum(p^2 / m)
	{
		template <having_coordinates<T, State> System>
		void step(System& s, T dt, bool first_step = false) const
		{
			kick(s, dt/2, first_step);
			s.x += s.vel() * dt;
			kick(s, dt/2);

			s.t += dt;
		}

		private:

			template <having_coordinates<T, State> System>
			requires requires(System& s, T t, std::size_t i)
				{
					t += dot(s.p[i], s.p[i]) / s.m[i];
					t += dot(s.p[i], s.f[i]) / s.m[i];
					t += dot(s.f[i], s.f[i]) / s.m[i];
				}
			void kick(System& s, T tau, bool eval = true) const
			{
				using std::size_t;
				using std::sqrt;
				using std::sinh;
				using std::cosh;

				T den = 0, xi0 = 0, omega02 = 0, omega0, a, b;
				for (size_t i = 0; i < s.n; ++i)
					den += dot(s.p[i], s.p[i]) / s.m[i];
				s.force(eval);
				if (den)
				{
					for (size_t i = 0; i < s.n; ++i)
						xi0 += dot(s.p[i], s.f[i]) / s.m[i];
					xi0 /= den;
					for (size_t i = 0; i < s.n; ++i)
						omega02 += dot(s.f[i], s.f[i]) / s.m[i];
					omega02 /= den;
					omega0 = sqrt(omega02);

					if (omega0)
					{
						a = cosh(omega0*tau) + xi0*sinh(omega0*tau)/omega0;
						b = sinh(omega0*tau)/omega0 + xi0*(cosh(omega0*tau)-1)/omega02;
					}
					else
					{
						a = 1 + xi0*tau;
						b = tau*(1 + xi0*tau/2);
					}

					s.p = (s.p + s.force(false) * b) / a;
				}
				else
					s.p += s.force(false) * tau;
			}
	};

	template <typename T, typename State>
	struct nose_hoover : integrator_base<T, State>
	// Nosé-Hoover thermostats chain integrator (2nd order, 1 stage)
	// It approximates a canonical (NVT) ensemble
	// See Allen, Computer simulation of liquids, 2017, pp. 134-139
	{
		nose_hoover(std::size_t n_th = 10, T tau = .2) : tau_relax(tau), n_th(n_th), p_th(n_th), m_th(n_th)
		{}

		template <having_coordinates<T, State> System>
		requires requires(System& s, T t)
			{
				t += 2*s.kinetic_energy() - s.dof*s.kT_ref();
			}
		void step(System& s, T dt, bool first_step = false)
		{
			using std::size_t;
			using std::exp;
			if (n_th == 0)
				throw std::runtime_error("Error: number of thermostats must be greater than 0.");

			p_th.resize(n_th);
			m_th.resize(n_th);

			T tau2 = s.kT_ref() * tau_relax*tau_relax;
			m_th[0] = s.dof * tau2;
			for (size_t i = 1; i < n_th; ++i)
				m_th[i] = tau2;

			// U4
			p_th[n_th-1] += G(s,n_th-1) * (dt/2);
			for (size_t i = n_th-1; i --> 0; )
				if (p_th[i+1])
				{
					T a = exp(-xi(i+1)*(dt/2));
					p_th[i] = a * p_th[i] + (1 - a) * G(s,i) / xi(i+1);
				}
				else
					p_th[i] += G(s,i) * (dt/2);
			// U3
			s.p *= exp(-xi(0) * (dt/2));
			// U2
			s.p += s.force(first_step) * (dt/2);
			// U1
			s.x += s.vel() * dt;
			// U2
			s.p += s.force() * (dt/2);
			// U3
			s.p *= exp(-xi(0) * (dt/2));
			// U4
			for (size_t i = 0; i < n_th-1; ++i)
				if (p_th[i+1])
				{
					T a = exp(-xi(i+1)*(dt/2));
					p_th[i] = a * p_th[i] + (1 - a) * G(s,i) / xi(i+1);
				}
				else
					p_th[i] += G(s,i) * (dt/2);
			p_th[n_th-1] += G(s,n_th-1) * (dt/2);

			s.t += dt;
		}

		T tau_relax;
		std::size_t n_th;

		private:

			std::vector<T> p_th, m_th;

			T xi(std::size_t i) const
			{
				return p_th[i]/m_th[i];
			}

			template <having_coordinates<T, State> System>
			T G(System& s, std::size_t i) const
			{
				if (i == 0)
					return 2*s.kinetic_energy() - s.dof*s.kT_ref();
				else
					return p_th[i-1]*xi(i-1) - s.kT_ref();
			}
	};

	template <typename T, typename State>
	struct martyna_tobias_klein : integrator_base<T, State>
	// Martyna-Tobias-Klein equations integrator (2nd order, 1 stage)
	// It approximates an isothermal-isobaric (NPT) ensemble
	// See Allen et al., "Computer simulation of liquids", 2017, pp. 140-144
	// and Tuckerman et al., "A Liouville-operator derived measure-preserving integrator
	// for molecular dynamics simulations in the isothermal-isobaric ensemble", 2006
	// Here, the velocity of the particles is assumed to be the same as p/m instead of dx/dt
	// (they are not the same due to position rescaling).
	{
		martyna_tobias_klein(std::size_t n_th = 10, T tau_th = T(.2), T tau_ba = 5)
			: tau_th(tau_th), tau_ba(tau_ba), n_th(n_th), p_th(n_th), m_th(n_th), p_ba(n_th), m_ba(n_th)
		{}

		template <having_coordinates<T, State> System>
		requires requires(System& s, T t)
			{
				t += 2*s.kinetic_energy() - s.dof*s.kT_ref();
				t += ((s.pressure() - s.pressure_ref)*s.volume() + 2 * s.kinetic_energy() / s.n)*t;
				s.side *= t;
			}
		void step(System& s, T dt, bool first_step = false)
		{
			using std::size_t;
			using std::exp;
			using std::log;
			if (n_th == 0)
				throw std::runtime_error("Error: number of thermostats must be greater than 0.");

			p_th.resize(n_th);
			m_th.resize(n_th);
			p_ba.resize(n_th);
			m_ba.resize(n_th);

			T tau2 = s.kT_ref() * tau_th*tau_th;
			m_th[0] = s.dof * tau2;
			for (size_t i = 1; i < n_th; ++i)
				m_th[i] = tau2;

			tau2 = s.kT_ref() * tau_ba*tau_ba;
			m_strain = s.dof * tau2;
			for (size_t i = 0; i < n_th; ++i)
				m_ba[i] = tau2;

			// U4
			U4_forward(s, dt/2);
			// U3
			s.p *= exp(-xi_th(0) * (dt/2));
			p_strain *= exp(-xi_ba(0) * (dt/2));
			// U2
			p_strain += G_strain(s) * (dt/2);
			if (p_strain)
			{
				T ax = (1 + T(1)/s.n) * xi_strain();
				T a = exp(-ax*(dt/2));
				s.p = a * s.p + (1 - a) / ax * s.force(first_step);
			}
			else
				s.p += s.force(first_step) * (dt/2);
			// U1
			if (p_strain)
			{
				T a = exp(xi_strain()*dt);
				s.x = a * s.x + (a - 1) / xi_strain() * s.vel();
				s.side *= a;
			}
			else
				s.x += s.vel() * dt;
			// U2
			if (p_strain)
			{
				T ax = (1 + T(1)/s.n) * xi_strain();
				T a = exp(-ax*(dt/2));
				s.p = a * s.p + (1 - a) / ax * s.force();
			}
			else
				s.p += s.force() * (dt/2);
			p_strain += G_strain(s) * (dt/2);
			// U3
			p_strain *= exp(-xi_ba(0) * (dt/2));
			s.p *= exp(-xi_th(0) * (dt/2));
			// U4
			U4_backward(s, dt/2);

			s.t += dt;
		}

		T tau_th, tau_ba;
		std::size_t n_th;

		private:

			std::vector<T> p_th, m_th, p_ba, m_ba;
			T p_strain, m_strain;

			T xi_th(std::size_t i) const
			{
				return p_th[i]/m_th[i];
			}

			T xi_ba(std::size_t i) const
			{
				return p_ba[i]/m_ba[i];
			}

			T xi_strain() const
			{
				return p_strain/m_strain;
			}

			template <having_coordinates<T, State> System>
			T G_th(System& s, std::size_t i) const
			{
				if (i == 0)
					return 2*s.kinetic_energy() - s.dof*s.kT_ref();
				else
					return p_th[i-1]*xi_th(i-1) - s.kT_ref();
			}

			template <having_coordinates<T, State> System>
			T G_ba(const System& s, std::size_t i) const
			{
				if (i == 0)
					return p_strain*xi_strain() - s.kT_ref();
				else
					return p_ba[i-1]*xi_ba(i-1) - s.kT_ref();
			}

			template <having_coordinates<T, State> System>
			T G_strain(System& s) const
			{
				return (s.dof*(s.pressure() - s.pressure_ref)*s.volume() + 2 * s.kinetic_energy()) / s.n;
			}

			template <having_coordinates<T, State> System>
			void U4_forward(System& s, T deltat)
			{
				// U4 thermostats coupled to particles momenta
				p_th[n_th-1] += G_th(s,n_th-1) * deltat;
				for (size_t i = n_th-1; i --> 0; )
					if (p_th[i+1])
					{
						T a = exp(-xi_th(i+1)*deltat);
						p_th[i] = a * p_th[i] + (1 - a) * G_th(s,i) / xi_th(i+1);
					}
					else
						p_th[i] += G_th(s,i) * deltat;
				// U4 thermostats coupled to barostat momentum
				p_ba[n_th-1] += G_ba(s,n_th-1) * deltat;
				for (size_t i = n_th-1; i --> 0; )
					if (p_ba[i+1])
					{
						T a = exp(-xi_ba(i+1)*deltat);
						p_ba[i] = a * p_ba[i] + (1 - a) * G_ba(s,i) / xi_ba(i+1);
					}
					else
						p_ba[i] += G_ba(s,i) * deltat;
			}

			template <having_coordinates<T, State> System>
			void U4_backward(System& s, T deltat)
			{
				// U4 thermostats coupled to barostat momentum
				for (size_t i = 0; i < n_th-1; ++i)
					if (p_ba[i+1])
					{
						T a = exp(-xi_ba(i+1)*deltat);
						p_ba[i] = a * p_ba[i] + (1 - a) * G_ba(s,i) / xi_ba(i+1);
					}
					else
						p_ba[i] += G_ba(s,i) * deltat;
				p_ba[n_th-1] += G_ba(s,n_th-1) * deltat;
				// U4 thermostats coupled to particles momenta
				for (size_t i = 0; i < n_th-1; ++i)
					if (p_th[i+1])
					{
						T a = exp(-xi_th(i+1)*deltat);
						p_th[i] = a * p_th[i] + (1 - a) * G_th(s,i) / xi_th(i+1);
					}
					else
						p_th[i] += G_th(s,i) * deltat;
				p_th[n_th-1] += G_th(s,n_th-1) * deltat;
			}
	};

	template <typename T, typename State>
	struct pefrl : symplectic_integrator_base<T, State>
	// Position-extended Forest-Ruth-like (4th order, 4 stages)
	// OMELYAN, MRYGLOD, FOLK
	// OPTIMIZED FOREST-RUTH- AND SUZUKI-LIKE ALGORITHMS FOR INTEGRATION
	// OF MOTION IN MANY-BODY SYSTEMS
	// 2008
	// important notice: the final force (and potential energy)
	// will not be synchronized with the final state of the system,
	// so energy may need to be computed separately from the force.
	// The VEFRL variant does not have this annoyance.
	{
		template <having_coordinates<T, State> System>
		void step(System& s, T dt, bool first_step = false) const
		{
			s.x += s.vel(first_step) * (xi * dt);
			s.p += s.force() * ((1 - 2*lambda) * dt/2);
			s.x += s.vel() * (chi * dt);
			s.p += s.force() * (lambda * dt);
			s.x += s.vel() * ((1 - 2*(chi + xi)) * dt);
			s.p += s.force() * (lambda * dt);
			s.x += s.vel() * (chi * dt);
			s.p += s.force() * ((1 - 2*lambda) * dt/2);
			s.x += s.vel() * (xi * dt);

			s.t += dt;
		}

		private:

			static constexpr T xi = 0.1786178958448091L;
			static constexpr T lambda = -0.2123418310626054L;
			static constexpr T chi = -0.6626458266981849e-1L;
	};

	template <typename T, typename State>
	struct vefrl : symplectic_integrator_base<T, State>
	// Velocity-extended Forest-Ruth-like (4th order, 4 stages)
	// OMELYAN, MRYGLOD, FOLK
	// OPTIMIZED FOREST-RUTH- AND SUZUKI-LIKE ALGORITHMS FOR INTEGRATION
	// OF MOTION IN MANY-BODY SYSTEMS
	// 2008
	{
		template <having_coordinates<T, State> System>
		void step(System& s, T dt, bool first_step = false) const
		{
			s.p += s.force(first_step) * (xi * dt);
			s.x += s.vel() * ((1 - 2*lambda) * dt/2);
			s.p += s.force() * (chi * dt);
			s.x += s.vel() * (lambda * dt);
			s.p += s.force() * ((1 - 2*(chi + xi)) * dt);
			s.x += s.vel() * (lambda * dt);
			s.p += s.force() * (chi * dt);
			s.x += s.vel() * ((1 - 2*lambda) * dt/2);
			s.p += s.force() * (xi * dt);

			s.t += dt;
		}

		private:

			static constexpr T xi = 0.1644986515575760L;
			static constexpr T lambda = -0.2094333910398989e-1L;
			static constexpr T chi = 0.1235692651138917e+1L;
	};

} // namespace physics

#endif // PHYSICS_INTEGRATOR_INTEGRATOR_H






























