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

#include <vector>
#include <cmath> // exp, sqrt

#include "integrator_base.hpp"

namespace physics
{
	template <having_coordinates System>
	struct symplectic_euler : symplectic_integrator_base<System>
	// Symplectic Euler method (1st order, 1 stage)
	// It is equivalent to leapfrog after correcting the first and the last steps
	{
		using symplectic_integrator_base<System>::first_step;

		System& s;

		symplectic_euler(System& s) : s(s) {}
		// constructor:
		// 	`s` is the system to integrate.

		void step(scalar_type_of<System> dt) override
		{
			s.p += s.force(first_step) * dt;
			s.x += s.vel() * dt;

			s.force();

			s.t += dt;
			first_step = false;
		}
	};

	template <having_coordinates System>
	struct leapfrog : symplectic_integrator_base<System>
	// Leapfrog method (2nd order, 1 stage)
	{
		using symplectic_integrator_base<System>::first_step;

		System& s;

		leapfrog(System& s) : s(s) {}
		// constructor:
		// 	`s` is the system to integrate.

		void step(scalar_type_of<System> dt) override
		{
			s.p += s.force(first_step) * (dt/2);
			s.x += s.vel() * dt;
			s.p += s.force() * (dt/2);

			s.t += dt;
			first_step = false;
		}
	};

	template <having_coordinates System>
	requires requires(System& s, scalar_type_of<System> dt)
		{
			s.p += s.force_long() * dt;
			s.p += s.force_short() * dt;
			s.p += s.force_long(false) * dt;
		}
	struct multi_timestep_leapfrog : symplectic_integrator_base<System>
	// Multi-timestep leapfrog method (2nd order, 1 long stage, n_short short stages)
	// If the force is divided into a high frequency and a low frequency part,
	// and the low frequency is much more expensive to compute, then this method
	// may be useful to reduce compute time or to increase accuracy.
	// In the microcanonical (NVE) ensemble this method can better conserve the
	// total energy without reducing the timestep
	{
		using symplectic_integrator_base<System>::first_step;

		System& s;

		multi_timestep_leapfrog(System& s, std::size_t n_short = 10) : s(s), n_short(n_short)
		// constructor:
		// 	`s` is the system to integrate.
		//	`n_short` is the number of "short" stages during one "long" stage
		{}

		void step(scalar_type_of<System> dt) override
		{
			scalar_type_of<System> dt_short = dt/n_short;
			s.p += s.force_long(first_step) * (dt/2);
			s.p += s.force_short() * (dt_short/2);
			for (std::size_t i = 0; i < n_short-1; ++i)
			{
				s.x += s.vel() * dt_short;
				s.p += s.force_short() * dt_short;
			}
			s.x += s.vel() * dt_short;
			s.p += s.force_short() * (dt_short/2);
			s.p += s.force_long() * (dt/2);

			s.t += dt;
			first_step = false;
		}

		private:

			std::size_t n_short;
	};

	template <having_coordinates System>
	struct pefrl : symplectic_integrator_base<System>
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
		using symplectic_integrator_base<System>::first_step;

		System &s;

		pefrl(System& s) : s(s) {}
		// constructor:
		// 	`s` is the system to integrate.

		void step(scalar_type_of<System> dt) override
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
			first_step = false;
		}

		private:

			static constexpr scalar_type_of<System> xi = 0.1786178958448091L;
			static constexpr scalar_type_of<System> lambda = -0.2123418310626054L;
			static constexpr scalar_type_of<System> chi = -0.6626458266981849e-1L;
	};

	template <having_coordinates System>
	struct vefrl : symplectic_integrator_base<System>
	// Velocity-extended Forest-Ruth-like (4th order, 4 stages)
	// OMELYAN, MRYGLOD, FOLK
	// OPTIMIZED FOREST-RUTH- AND SUZUKI-LIKE ALGORITHMS FOR INTEGRATION
	// OF MOTION IN MANY-BODY SYSTEMS
	// 2008
	{
		using symplectic_integrator_base<System>::first_step;

		System& s;

		vefrl(System& s) : s(s) {}
		// constructor:
		// 	`s` is the system to integrate.

		void step(scalar_type_of<System> dt) override
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
			first_step = false;
		}

		private:

			static constexpr scalar_type_of<System> xi = 0.1644986515575760L;
			static constexpr scalar_type_of<System> lambda = -0.2094333910398989e-1L;
			static constexpr scalar_type_of<System> chi = 0.1235692651138917e+1L;
	};

	template <typename System>
	concept having_coordinates_damped = having_coordinates<System>
		&& requires(System& s, scalar_type_of<System> dt, std::size_t i)
	// A system `having_coordinates_damped` requires to be a `having_coordinates`,
	// and has `gamma` (the damping coefficient, which must be dereferenceable) defined.
	// The following operations are required to be defined.
		{
			{s.gamma[i]} -> std::convertible_to<scalar_type_of<System>>;
			s.p[i] = s.p[i] * dt;
			i < s.n;
		};

	template <having_coordinates_damped System>
	struct damped_leapfrog : integrator_base<System>
	// damped integrator (2nd order?, 1 stage)
	// dp = f dt - gamma p dt
	{
		using integrator_base<System>::first_step;

		System& s;

		damped_leapfrog(System& s) : s(s) {}
		// constructor:
		// 	`s` is the system to integrate.

		void step(scalar_type_of<System> dt) override
		{
			using std::exp;

			s.p += s.force(first_step) * (dt/2);
			s.x += s.vel() * (dt/2);
			for (std::size_t i = 0; i < s.n; ++i)
				s.p[i] *= exp(-s.gamma[i] * dt);
			s.x += s.vel() * (dt/2);
			s.p += s.force() * (dt/2);

			s.t += dt;
			first_step = false;
		}
	};

	template <typename System>
	concept having_coordinates_stochastic = having_coordinates_damped<System>
		&& requires(System& s, scalar_type_of<System> dt, std::size_t i)
	// A system `having_coordinates_stochastic` requires to be a `having_coordinates_damped`,
	// and has `D` (the diffusion coefficient, which is a scalar), `noise` and `rand` defined. `rand`
	// is a method that must generate random numbers for the `noise` vector (which must be
	// deferenceable and its elements should be convertible to s.p[i]). In particular, for
	// Langevin's equations, the following should be defined: s.noise = sqrt(m k T) * gaussian(0, 1).
	// The following operations are required to be defined.
		{
			{s.D} -> std::convertible_to<scalar_type_of<System>>;
			s.p[i] = s.p[i] * dt + s.noise[i] * dt;
			s.rand();
		};

	template <having_coordinates_stochastic System>
	struct stochastic_leapfrog : stochastic_integrator_base<System>
	// stochastic leapfrog (1 stage)
	// using the following equations:
	// dp = f dt - gamma p dt + sigma dw
	// where sigma dw = sqrt(2 gamma) * noise
	// For Langevin's equations, the following equations must be satisfied:
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
		using stochastic_integrator_base<System>::first_step;

		System& s;

		stochastic_leapfrog(System& s) : s(s) {}
		// constructor:
		// 	`s` is the system to integrate.

		void step(scalar_type_of<System> dt) override
		{
			using std::exp;
			using std::sqrt;

			s.p += s.force(first_step) * (dt/2);
			s.x += s.vel() * (dt/2);
			s.rand();
			if (s.D == 0) // i.e., gamma -> inf
				s.p = s.noise; // s.noise = sqrt(m k T) * gaussian(0, 1)
			else
				for (std::size_t i = 0; i < s.n; ++i)
					s.p[i] = exp(-s.gamma[i] * dt) * s.p[i] + sqrt(1 - exp(-2 * s.gamma[i] * dt)) * s.noise[i];
				// p = exp(-gamma dt) p + sqrt((1 - exp(-2 gamma dt)) m k T) g
				// expanding the exp (gamma small): p += -gamma p dt + sqrt(2 gamma m k T dt) g
				// g: gaussian noise with mean = 0, std.dev = 1
			s.x += s.vel() * (dt/2);
			s.p += s.force() * (dt/2);

			s.t += dt;
			first_step = false;
		}
	};

	template <having_coordinates System>
	requires requires(System& s, scalar_type_of<System> t, std::size_t i)
		{
			t += dot(s.p[i], s.p[i]) / s.m[i];
			t += dot(s.p[i], s.f[i]) / s.m[i];
			t += dot(s.f[i], s.f[i]) / s.m[i];
		}
	struct isokinetic_leapfrog : integrator_base<System>
	// (gaussian) isokinetic integrator (2nd order, 1 stage)
	// it assumes quadratic kinetic energy
	// not symplectic, but still time-reversible
	// it conserves the kinetic energy rather than the total hamiltonian
	// at equilibrium, its configurations sample a canonical ensemble
	// unfortunately, momenta do not follow the MB distribution
	// dp = f dt - xi p dt
	// xi = sum(p * f / m) / sum(p^2 / m)
	{
		using integrator_base<System>::first_step;

		System& s;

		isokinetic_leapfrog(System& s) : s(s) {}
		// constructor:
		// 	`s` is the system to integrate.

		void step(scalar_type_of<System> dt) override
		{
			kick(dt/2, first_step);
			s.x += s.vel() * dt;
			kick(dt/2);

			s.t += dt;
			first_step = false;
		}

		private:

			void kick(scalar_type_of<System> tau, bool eval = true) const
			{
				using std::size_t;
				using std::sqrt;
				using std::sinh;
				using std::cosh;

				scalar_type_of<System> den = 0, xi0 = 0, omega02 = 0, omega0, a, b;
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

	template <having_coordinates System>
	requires requires(System& s, scalar_type_of<System> t)
		{
			t = 2*s.kinetic_energy() - s.dof*s.kT_ref();
		}
	struct nose_hoover : integrator_base<System>
	// Nosé-Hoover thermostats chain integrator (2nd order, 1 stage)
	// It approximates a canonical (NVT) ensemble.
	// See Allen, Computer simulation of liquids, 2017, pp. 134-139
	{
		using integrator_base<System>::first_step;

		System& s;

		nose_hoover(System& s, std::size_t n_th = 10, scalar_type_of<System> tau = .2L)
		// constructor:
		// 	`s` is the system to integrate.
		//	`n_th` is the number of thermostats in the chain
		//	`tau` is the characteristic time of the oscillations
			: s(s), tau_relax(tau), n_th(n_th), p_th(n_th), m_th(n_th)
		{
			for (auto& p : p_th)
				p = 0;
		}

		void step(scalar_type_of<System> dt) override
		// throw a `std::runtime_error` if member variable `n_th` is equal to 0.
		// (n_th == 0 would be the same as leapfrog integration).
		{
			using std::size_t;
			using std::exp;
			if (n_th == 0)
				throw std::runtime_error("Error: number of thermostats must be greater than 0.");

			p_th.resize(n_th);
			m_th.resize(n_th);

			scalar_type_of<System> tau2 = s.kT_ref() * tau_relax*tau_relax;
			m_th[0] = s.dof * tau2;
			for (size_t i = 1; i < n_th; ++i)
				m_th[i] = tau2;

			// U4
			p_th[n_th-1] += G(n_th-1) * (dt/2);
			for (size_t i = n_th-1; i --> 0; )
				if (p_th[i+1])
				{
					scalar_type_of<System> a = exp(-xi(i+1)*(dt/2));
					p_th[i] = a * p_th[i] + (1 - a) * G(i) / xi(i+1);
				}
				else
					p_th[i] += G(i) * (dt/2);
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
					scalar_type_of<System> a = exp(-xi(i+1)*(dt/2));
					p_th[i] = a * p_th[i] + (1 - a) * G(i) / xi(i+1);
				}
				else
					p_th[i] += G(i) * (dt/2);
			p_th[n_th-1] += G(n_th-1) * (dt/2);

			s.t += dt;
			first_step = false;
		}

		scalar_type_of<System> tau_relax;
		std::size_t n_th;

		private:

			std::vector<scalar_type_of<System>> p_th, m_th;

			scalar_type_of<System> xi(std::size_t i) const
			{
				return p_th[i]/m_th[i];
			}

			scalar_type_of<System> G(std::size_t i) const
			{
				if (i == 0)
					return 2*s.kinetic_energy() - s.dof*s.kT_ref();
				else
					return p_th[i-1]*xi(i-1) - s.kT_ref();
			}
	};

	template <having_coordinates System>
	requires requires(System& s, scalar_type_of<System> t)
		{
			t = 2*s.kinetic_energy() - s.dof*s.kT_ref();
			t = (s.dof*(s.pressure() - s.pressure_ref)*s.volume() + 2 * s.kinetic_energy()) / s.n;
			s.side *= t;
		}
	struct mtk : integrator_base<System>
	// Martyna-Tobias-Klein equations integrator (2nd order, 1 stage)
	// It approximates an isothermal-isobaric (NPT) ensemble if n_th > 0.
	// See Allen et al., "Computer simulation of liquids", 2017, pp. 140-144
	// and Tuckerman et al., "A Liouville-operator derived measure-preserving integrator
	// for molecular dynamics simulations in the isothermal-isobaric ensemble", 2006
	// Here, the velocity of the particles is assumed to be the same as p/m instead of dx/dt
	// (they are not the same due to position rescaling).
	{
		using integrator_base<System>::first_step;

		System& s;

		mtk(System &s, std::size_t n_th = 10, scalar_type_of<System> tau_th = .2L, scalar_type_of<System> tau_ba = 5)
		// constructor:
		// 	`s` is the system to integrate.
		//	`n_th` is the number of thermostats in the chain.
		//	`tau_th` is the characteristic time of the oscillations of the thermostats coupled
		//		to particles momenta.
		//	`tau_ba` is the characteristic time of the oscillations of the thermostats coupled
		//		to the barostat momentum and of the barostat.
			: s(s), tau_th(tau_th), tau_ba(tau_ba), n_th(n_th), p_th(n_th), m_th(n_th), p_ba(n_th), m_ba(n_th)
		{
			for (auto& p : p_th)
				p = 0;
			for (auto& p : p_ba)
				p = 0;
		}

		void step(scalar_type_of<System> dt) override
		{
			using std::size_t;
			using std::exp;

			p_th.resize(n_th == 0 ? 1 : n_th);
			m_th.resize(n_th == 0 ? 1 : n_th);
			p_ba.resize(n_th == 0 ? 1 : n_th);
			m_ba.resize(n_th == 0 ? 1 : n_th);

			scalar_type_of<System> tau2 = s.kT_ref() * tau_th*tau_th;
			m_th[0] = s.dof * tau2;
			for (size_t i = 1; i < n_th; ++i)
				m_th[i] = tau2;

			tau2 = s.kT_ref() * tau_ba*tau_ba;
			m_strain = s.dof * tau2;
			for (size_t i = 0; i < n_th; ++i)
				m_ba[i] = tau2;

			if (n_th > 0)
			{
				// U4
				U4_forward(dt/2);
				// U3
				s.p *= exp(-xi_th(0) * (dt/2));
				p_strain *= exp(-xi_ba(0) * (dt/2));
			}
			// U2
			p_strain += G_strain() * (dt/2);
			if (p_strain)
			{
				scalar_type_of<System> ax = (1 + scalar_type_of<System>(1)/s.n) * xi_strain();
				scalar_type_of<System> a = exp(-ax*(dt/2));
				s.p = a * s.p + (1 - a) / ax * s.force(first_step);
			}
			else
				s.p += s.force(first_step) * (dt/2);
			// U1
			if (p_strain)
			{
				scalar_type_of<System> a = exp(xi_strain()*dt);
				s.x = a * s.x + (a - 1) / xi_strain() * s.vel();
				s.side *= a;
			}
			else
				s.x += s.vel() * dt;
			// U2
			if (p_strain)
			{
				scalar_type_of<System> ax = (1 + scalar_type_of<System>(1)/s.n) * xi_strain();
				scalar_type_of<System> a = exp(-ax*(dt/2));
				s.p = a * s.p + (1 - a) / ax * s.force();
			}
			else
				s.p += s.force() * (dt/2);
			p_strain += G_strain() * (dt/2);

			if (n_th > 0)
			{
				// U3
				p_strain *= exp(-xi_ba(0) * (dt/2));
				s.p *= exp(-xi_th(0) * (dt/2));
				// U4
				U4_backward(dt/2);
			}

			s.t += dt;
			first_step = false;
		}

		scalar_type_of<System> tau_th, tau_ba;
		std::size_t n_th;

		private:

			std::vector<scalar_type_of<System>> p_th, m_th, p_ba, m_ba;
			scalar_type_of<System> p_strain = 0, m_strain = 0;

			scalar_type_of<System> xi_th(std::size_t i) const
			{
				return p_th[i]/m_th[i];
			}

			scalar_type_of<System> xi_ba(std::size_t i) const
			{
				return p_ba[i]/m_ba[i];
			}

			scalar_type_of<System> xi_strain() const
			{
				return p_strain/m_strain;
			}

			scalar_type_of<System> G_th(std::size_t i) const
			{
				if (i == 0)
					return 2*s.kinetic_energy() - s.dof*s.kT_ref();
				else
					return p_th[i-1]*xi_th(i-1) - s.kT_ref();
			}

			scalar_type_of<System> G_ba(std::size_t i) const
			{
				if (i == 0)
					return p_strain*xi_strain() - s.kT_ref();
				else
					return p_ba[i-1]*xi_ba(i-1) - s.kT_ref();
			}

			scalar_type_of<System> G_strain() const
			{
				return (s.dof*(s.pressure() - s.pressure_ref)*s.volume() + 2 * s.kinetic_energy()) / s.n;
			}

			void U4_forward(scalar_type_of<System> deltat)
			{
				using std::exp;
				// U4 thermostats coupled to particles momenta
				p_th[n_th-1] += G_th(n_th-1) * deltat;
				for (size_t i = n_th-1; i --> 0; )
					if (p_th[i+1])
					{
						scalar_type_of<System> a = exp(-xi_th(i+1)*deltat);
						p_th[i] = a * p_th[i] + (1 - a) * G_th(i) / xi_th(i+1);
					}
					else
						p_th[i] += G_th(i) * deltat;
				// U4 thermostats coupled to barostat momentum
				p_ba[n_th-1] += G_ba(n_th-1) * deltat;
				for (size_t i = n_th-1; i --> 0; )
					if (p_ba[i+1])
					{
						scalar_type_of<System> a = exp(-xi_ba(i+1)*deltat);
						p_ba[i] = a * p_ba[i] + (1 - a) * G_ba(i) / xi_ba(i+1);
					}
					else
						p_ba[i] += G_ba(i) * deltat;
			}

			void U4_backward(scalar_type_of<System> deltat)
			{
				using std::exp;
				// U4 thermostats coupled to barostat momentum
				for (size_t i = 0; i < n_th-1; ++i)
					if (p_ba[i+1])
					{
						scalar_type_of<System> a = exp(-xi_ba(i+1)*deltat);
						p_ba[i] = a * p_ba[i] + (1 - a) * G_ba(i) / xi_ba(i+1);
					}
					else
						p_ba[i] += G_ba(i) * deltat;
				p_ba[n_th-1] += G_ba(n_th-1) * deltat;
				// U4 thermostats coupled to particles momenta
				for (size_t i = 0; i < n_th-1; ++i)
					if (p_th[i+1])
					{
						scalar_type_of<System> a = exp(-xi_th(i+1)*deltat);
						p_th[i] = a * p_th[i] + (1 - a) * G_th(i) / xi_th(i+1);
					}
					else
						p_th[i] += G_th(i) * deltat;
				p_th[n_th-1] += G_th(n_th-1) * deltat;
			}
	};

} // namespace physics

#endif // PHYSICS_INTEGRATOR_INTEGRATOR_H






























