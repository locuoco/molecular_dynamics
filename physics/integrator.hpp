//  Integrators
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

#ifndef PHYSICS_INTEGRATOR_H
#define PHYSICS_INTEGRATOR_H

#include <concepts> // convertible_to
#include <valarray>
#include <type_traits> // is_base_of
#include <cmath> // exp, sqrt

#include "point.hpp"

namespace physics
{
	template <typename T, typename State>
	struct integrator_base {};

	template <typename T, typename State>
	struct symplectic_integrator_base : integrator_base<T, State> {};
	// Assumption: H = T(p) + V(x), with f = -grad V(x) and v = grad T(p)
	// Furthermore, H is explicitly independent on time
	// unless stated otherwise

	template <typename T, typename State>
	struct stochastic_integrator_base : integrator_base<T, State> {};

	template <typename Integ, typename T, typename State>
	concept integrator = std::is_base_of_v<integrator_base<T, State>, Integ>;

	template <typename Integ, typename T, typename State>
	concept symplectic_integrator = std::is_base_of_v<symplectic_integrator_base<T, State>, Integ>;

	template <typename Integ, typename T, typename State>
	concept stochastic_integrator = std::is_base_of_v<stochastic_integrator_base<T, State>, Integ>;

	template <typename T, typename State>
	struct physical_system_base
	{
		virtual const State& force(bool = true) = 0;
		virtual const State& vel(bool = true) = 0;
	};

	template <typename S, typename T, typename State>
	concept physical_system = std::is_base_of_v<physical_system_base<T, State>, S>;

	template <typename S, typename T, typename State>
	concept having_coordinates = physical_system<S, T, State>
		&& requires(S& s, T dt)
		{
			s.x += s.vel() * dt;
			s.p += s.force() * dt;
			s.t += dt;
		};

	template <typename S, typename T, typename State>
	concept having_coordinates_damped = having_coordinates<S, T, State>
		&& requires(S& s, T dt, std::size_t i)
		{
			{s.gamma[i]} -> std::convertible_to<T>;
			s.p[i] = s.p[i] * dt;
			i < s.n;
		};

	template <typename S, typename T, typename State>
	concept having_coordinates_stochastic = having_coordinates_damped<S, T, State>
		&& requires(S& s, T dt, std::size_t i)
		{
			{s.D} -> std::convertible_to<T>;
			s.p[i] = s.p[i] * dt + s.noise[i] * dt;
			s.rand();
		};

	template <typename T, typename State>
	struct symplectic_euler : symplectic_integrator_base<T, State>
	// Symplectic Euler method (1st order, 1 stage)
	// It is equivalent to leapfrog after correcting the initial conditions
	{
		template <having_coordinates<T, State> S>
		void step(S& s, T dt, bool = false) const
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
		template <having_coordinates<T, State> S>
		void step(S& s, T dt, bool first_step = false) const
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
	// will be useful to reduce compute time or to increase accuracy
	// In the microcanonical (NVE) ensemble this method can better conserve the
	// total energy without reducing the timestep
	{
		multi_timestep_leapfrog(std::size_t n_short = 10) : n_short(n_short)
		{}

		template <having_coordinates<T, State> S>
		requires requires(S& s, T dt)
			{
				s.p += s.force_long() * dt;
				s.p += s.force_short() * dt;
				s.p += s.force_long(false) * dt;
			}
		void step(S& s, T dt, bool first_step = false) const
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

	template <typename T, typename State>
	struct stochastic_leapfrog : stochastic_integrator_base<T, State>
	// stochastic leapfrog (2nd order?, 1 stage)
	// dp = f dt - gamma p dt + sigma dw
	// with:
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
		template <having_coordinates_stochastic<T, State> S>
		void step(S& s, T dt, bool first_step = false) const
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
	struct damped_leapfrog : integrator_base<T, State>
	// damped integrator (2nd order?, 1 stage)
	// dp = f dt - gamma p dt
	{
		template <having_coordinates_damped<T, State> S>
		void step(S& s, T dt, bool first_step = false) const
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
		template <having_coordinates<T, State> S>
		void step(S& s, T dt, bool first_step = false) const
		{
			kick(s, dt/2, first_step);
			s.x += s.vel() * dt;
			kick(s, dt/2);

			s.t += dt;
		}

		private:

			template <having_coordinates<T, State> S>
			requires requires(S& s, T t, std::size_t i)
				{
					t += dot(s.p[i], s.p[i]) / s.m[i];
					t += dot(s.p[i], s.f[i]) / s.m[i];
					t += dot(s.f[i], s.f[i]) / s.m[i];
				}
			void kick(S& s, T tau, bool eval = true) const
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
	{
		nose_hoover(std::size_t n_th = 10, T tau = 1) : p_th(n_th), m_th(n_th), tau_relax(tau), n_th(n_th)
		{}

		template <having_coordinates<T, State> S>
		requires requires(S& s, T t)
			{
				t += s.dof*(s.temperature() - s.temperature_fixed);
			}
		void step(S& s, T dt, bool first_step = false)
		{
			using std::size_t;
			using std::exp;

			T tau2 = s.temperature_fixed * tau_relax*tau_relax;
			m_th[0] = s.dof * tau2;
			for (size_t i = 1; i < n_th; ++i)
				m_th[i] = tau2;

			p_th[n_th-1] += G(s,n_th-1) * (dt/2);
			for (size_t i = n_th-1; i --> 0; )
			{
				if (p_th[i])
				{
					T a = exp(-xi(i+1)*(dt/2));
					p_th[i] = p_th[i] * a + G(s,i) * (1 - a) / xi(i+1);
				}
				else
					p_th[i] += G(s,i) * (dt/2);
			}
			s.p *= exp(-xi(0) * (dt/2));
			s.p += s.force(first_step) * (dt/2);
			s.x += s.vel() * dt;
			s.p += s.force() * (dt/2);
			s.p *= exp(-xi(0) * (dt/2));
			for (size_t i = 0; i < n_th-1; ++i)
			{
				if (p_th[i])
				{
					T a = exp(-xi(i+1)*(dt/2));
					p_th[i] = p_th[i] * a + G(s,i) * (1 - a) / xi(i+1);
				}
				else
					p_th[i] += G(s,i) * (dt/2);
			}
			p_th[n_th-1] += G(s,n_th-1) * (dt/2);

			s.t += dt;
		}

		private:

			std::vector<T> p_th, m_th;
			T tau_relax;
			std::size_t n_th;

			T xi(std::size_t i) const
			{
				return p_th[i]/m_th[i];
			}

			template <having_coordinates<T, State> S>
			T G(const S& s, std::size_t i) const
			{
				if (i == 0)
					return s.dof*(s.temperature() - s.temperature_fixed);
				else
					return p_th[i-1]*xi(i-1) - s.temperature_fixed;
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
	// so they may need to be recomputed. The VEFRL variant
	// does not have this annoyance.
	{
		template <having_coordinates<T, State> S>
		void step(S& s, T dt, bool first_step = false) const
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
		template <having_coordinates<T, State> S>
		void step(S& s, T dt, bool first_step = false) const
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

	// COMPOSITION SCHEMES
	// can be used as long as the 2nd order integrator for a given hamiltonian is known
	// useful to build higher-order symplectic methods
	// composing the 2nd-order isokinetic integrator should give higher-order isokinetic integrators

/*
	KAHAN, LI
	COMPOSITION CONSTANTS FOR RAISING THE ORDERS OF
	UNCONVENTIONAL SCHEMES FOR ORDINARY
	DIFFERENTIAL EQUATIONS
	1997

	HAIRER, LUBICH, WANNER
	GEOMETRIC NUMERICAL INTEGRATION
	STRUCTURE-PRESERVING ALGORITHMS FOR ORDINARY
	DIFFERENTIAL EQUATIONS
	2006
	P. 156-158
*/

	template <typename T, typename State, template <typename, typename> typename IntegT, std::size_t Stages>
	struct composition_scheme_base : IntegT<T, State>
	{
		template <typename ... Ts>
		composition_scheme_base(Ts ... pars) : d{T(pars)...} {}

		template <having_coordinates<T, State> S>
		void step(S& s, T dt, bool first_step = false) const
		{
			for (size_t i = 0; i < Stages; ++i)
				IntegT<T, State>::step(s, d[std::min(i, Stages-i-1)] * dt, first_step);
		}

		private:

			const std::array<T, (Stages+1)/2> d;
	};

	template <typename Integ, typename T, typename State, template <typename, typename> typename IntegT, std::size_t Stages>
	concept composition_scheme = std::is_base_of_v<composition_scheme_base<T, State, IntegT, Stages>, Integ>;

	template <typename T, typename State, template <typename, typename> typename IntegT = leapfrog>
	struct forest_ruth : composition_scheme_base<T, State, IntegT, 3>
	// Forest-Ruth method (4th order, 3 stages)
	// independently obtained by Yoshida
	{
		forest_ruth()
			: composition_scheme_base<T, State, IntegT, 3>
			{
				1.3512071919596576340476878089715L, // 1 / (2 - cbrt(2))
				-1.7024143839193152680953756179429L // -cbrt(2) / (2 - cbrt(2))
			}
		{}
	};

	template <typename T, typename State, template <typename, typename> typename IntegT = leapfrog>
	struct suzuki4 : composition_scheme_base<T, State, IntegT, 5>
	// Suzuki method (4th order, 5 stages)
	{
		suzuki4()
			: composition_scheme_base<T, State, IntegT, 5>
			{
				0.41449077179437573714235406286076L, // 1 / (4 - cbrt(4))
				0.41449077179437573714235406286076L, // 1 / (4 - cbrt(4))
				-0.65796308717750294856941625144305L // -cbrt(4) / (4 - cbrt(4))
			}
		{}
	};

	template <typename T, typename State, template <typename, typename> typename IntegT = leapfrog>
	struct kahan_li4a : composition_scheme_base<T, State, IntegT, 5>
	// Kahan-Li method (4th order, 5 stages)
	{
		kahan_li4a()
			: composition_scheme_base<T, State, IntegT, 5>
			{
				0.78867513459481288225457439025098L, // (3 + sqrt(3)) / 6
				0.21132486540518711774542560974902L, // (3 - sqrt(3)) / 6
				-1
			}
		{}
	};

	template <typename T, typename State, template <typename, typename> typename IntegT = leapfrog>
	struct kahan_li4b : composition_scheme_base<T, State, IntegT, 5>
	// Kahan-Li method (4th order, 5 stages)
	{
		kahan_li4b()
			: composition_scheme_base<T, State, IntegT, 5>
			{
				0.21132486540518711774542560974902L, // (3 - sqrt(3)) / 6
				0.78867513459481288225457439025098L, // (3 + sqrt(3)) / 6
				-1
			}
		{}
	};

	template <typename T, typename State, template <typename, typename> typename IntegT = leapfrog>
	struct yoshida6 : composition_scheme_base<T, State, IntegT, 7>
	// Yoshida method (6th order, 7 stages)
	{
		yoshida6()
			: composition_scheme_base<T, State, IntegT, 7>
			{
				0.78451361047755726382L,
				0.23557321335935813368L,
				-1.1776799841788710069L,
				1.3151863206839112189L
			}
		{}
	};

	template <typename T, typename State, template <typename, typename> typename IntegT = leapfrog>
	struct kahan_li6a : composition_scheme_base<T, State, IntegT, 9>
	// Kahan-Li method (6th order, 9 stages)
	{
		kahan_li6a()
			: composition_scheme_base<T, State, IntegT, 9>
			{
				0.39216144400731413928L,
				0.33259913678935943860L,
				-0.70624617255763935981L,
				0.082213596293550800230L,
				0.79854399093482996340L
			}
		{}
	};

	template <typename T, typename State, template <typename, typename> typename IntegT = leapfrog>
	struct kahan_li6b : composition_scheme_base<T, State, IntegT, 9>
	// Kahan-Li method (6th order, 9 stages)
	{
		kahan_li6b()
			: composition_scheme_base<T, State, IntegT, 9>
			{
				0.39103020330868478817L,
				0.33403728961113601749L,
				-0.70622728118756134346L,
				0.081877549648059445768L,
				0.79856447723936218406L
			}
		{}
	};

	template <typename T, typename State, template <typename, typename> typename IntegT = leapfrog>
	struct suzuki8 : composition_scheme_base<T, State, IntegT, 15>
	// Suzuki method (8th order, 15 stages)
	{
		suzuki8()
			: composition_scheme_base<T, State, IntegT, 15>
			{
				0.74167036435061295344822780L,
				-0.40910082580003159399730010L,
				0.19075471029623837995387626L,
				-0.57386247111608226665638773L,
				0.29906418130365592384446354L,
				0.33462491824529818378495798L,
				0.31529309239676659663205666L,
				-0.79688793935291635401978884L
			}
		{}
	};

	template <typename T, typename State, template <typename, typename> typename IntegT = leapfrog>
	struct kahan_li8a : composition_scheme_base<T, State, IntegT, 17>
	// Kahan-Li method (8th order, 17 stages)
	{
		kahan_li8a()
			: composition_scheme_base<T, State, IntegT, 17>
			{
				0.13020248308889008088L,
				0.56116298177510838456L,
				-0.38947496264484728641L,
				0.15884190655515560090L,
				-0.39590389413323757734L,
				0.18453964097831570709L,
				0.25837438768632204729L,
				0.29501172360931029887L,
				-0.60550853383003451170L
			}
		{}
	};

	template <typename T, typename State, template <typename, typename> typename IntegT = leapfrog>
	struct kahan_li8b : composition_scheme_base<T, State, IntegT, 17>
	// Kahan-Li method (8th order, 17 stages)
	{
		kahan_li8b()
			: composition_scheme_base<T, State, IntegT, 17>
			{
				0.12713692773487857916L,
				0.56170253798880269972L,
				-0.38253471994883018888L,
				0.16007605629464743119L,
				-0.40181637432680696673L,
				0.18736671654227849724L,
				0.26070870920779240570L,
				0.29039738812516162389L,
				-0.60607448323584816258L
			}
		{}
	};

	template <typename T, typename State, template <typename, typename> typename IntegT = leapfrog>
	struct kahan_li10a : composition_scheme_base<T, State, IntegT, 31>
	// Kahan-Li method (10th order, 31 stages)
	{
		kahan_li10a()
			: composition_scheme_base<T, State, IntegT, 31>
			{
				-0.48159895600253002870L,
				0.0036303931544595926879L,
				0.50180317558723140279L,
				0.28298402624506254868L,
				0.80702967895372223806L,
				-0.026090580538592205447L,
				-0.87286590146318071547L,
				-0.52373568062510581643L,
				0.44521844299952789252L,
				0.18612289547097907887L,
				0.23137327866438360633L,
				-0.52191036590418628905L,
				0.74866113714499296793L,
				0.066736511890604057532L,
				-0.80360324375670830316L,
				0.91249037635867994571L
			}
		{}
	};

	template <typename T, typename State, template <typename, typename> typename IntegT = leapfrog>
	struct kahan_li10b : composition_scheme_base<T, State, IntegT, 31>
	// Kahan-Li method (10th order, 31 stages)
	{
		kahan_li10b()
			: composition_scheme_base<T, State, IntegT, 31>
			{
				0.27338476926228452782L,
				0.44587846502560283997L,
				0.83219642847136307126L,
				-0.83396868554957942879L,
				0.27891843057015194293L,
				0.89032738045702532006L,
				0.056681514845245709418L,
				-0.85737420814978887722L,
				-0.46789492554836586111L,
				-0.47919009182398264249L,
				0.16724074680043708909L,
				-0.87443151263376143307L,
				-0.49873481853620165786L,
				0.58930536608974918851L,
				0.83458937790882729775L,
				0.28614352562198582747L
			}
		{}
	};

	template <typename T, typename State, template <typename, typename> typename IntegT = leapfrog>
	struct kahan_li10c : composition_scheme_base<T, State, IntegT, 33>
	// Kahan-Li method (10th order, 33 stages)
	{
		kahan_li10c()
			: composition_scheme_base<T, State, IntegT, 33>
			{
				0.070428877682658066880L,
				0.87415651735353949041L,
				0.055414604963802442707L,
				-0.066800477898797011598L,
				-0.62641308958799555593L,
				0.23682621087528762872L,
				-0.42221063403170054210L,
				0.24222942201040859249L,
				0.047374515478601436594L,
				0.54386826052472423338L,
				-0.93252230928447264311L,
				0.16960179883676464855L,
				0.71608567578450563608L,
				-0.80016730247310573512L,
				0.23778185292256770747L,
				-0.32330301550863943389L,
				0.95529818470370207691L
			}
		{}
	};

	template <typename T, typename State, template <typename, typename> typename IntegT = leapfrog>
	struct kahan_li10d : composition_scheme_base<T, State, IntegT, 33>
	// Kahan-Li method (10th order, 33 stages)
	{
		kahan_li10d()
			: composition_scheme_base<T, State, IntegT, 33>
			{
				0.12282427644721572094L,
				0.77644680890696440342L,
				0.14881514553734297479L,
				-0.17239125953506067249L,
				-0.54745995781852463787L,
				0.14512932327306927479L,
				-0.31564555153114460562L,
				0.12086865089833871979L,
				0.17910277517866344258L,
				0.44263408813993245949L,
				-0.81935337479593697464L,
				0.13445474141752884045L,
				0.64444239169016646538L,
				-0.71930149370201612557L,
				0.21036902497348664610L,
				-0.26908194941570516294L,
				0.83629272067135846284L
			}
		{}
	};

	// RUNGE-KUTTA SCHEMES

	template <typename T, typename State>
	struct runge_kutta_base {};

	template <typename T, typename State, std::size_t Stages>
	struct explicit_runge_kutta_base : runge_kutta_base<T, State>
	{
		template <typename ... Ts>
		explicit_runge_kutta_base(Ts ... pars) : pars{T(pars)...} {}

		template <having_coordinates<T, State> S>
		void step(S& s, T dt, bool first_step = false) const
		{
			using std::size_t;

			prev_x = s.x;
			prev_p = s.p;
			prev_t = s.t;

			kx[0] = s.vel();
			kp[0] = s.force(first_step);
			for (size_t i = 1; i < Stages; ++i)
			{
				s.x = 0;
				s.p = 0;
				// sum smaller contributes (in dt) first to minimize rounding errors
				for (size_t j = 0; j < i; ++j)
				{
					T mult = pars[(i-1)*Stages + j+1] * dt;
					s.x += kx[j] * mult;
					s.p += kp[j] * mult;
				}
				s.x += prev_x;
				s.p += prev_p;
				s.t = prev_t + pars[(i-1)*Stages] * dt;

				kx[i] = s.vel();
				kp[i] = s.force();
			}

			s.x = 0;
			s.p = 0;
			for (size_t i = 0; i < Stages; ++i)
			{
				T mult = pars[(Stages-1)*Stages + i] * dt;
				s.x += kx[i] * mult;
				s.p += kp[i] * mult;
			}
			s.x += prev_x;
			s.p += prev_p;
			s.t = prev_t + dt;
			s.force();
		}

		private:

			const std::array<T, Stages*Stages> pars;
			State kx[Stages], kp[Stages], prev_x, prev_p;
			T prev_t;
	};

	template <typename T, typename State>
	struct euler : explicit_runge_kutta_base<T, State, 1>
	// Euler method (1st order, 1 stage)
	{
		euler()
			: explicit_runge_kutta_base<T, State, 1>{1}
		{}
	};

	template <typename T, typename State>
	struct rk2 : explicit_runge_kutta_base<T, State, 2>
	// Parametrized Runge-Kutta 2 method (2nd order, 2 stages)
	{
		rk2(long double a)
			: explicit_runge_kutta_base<T, State, 2>
			{
				a, a,
				1 - 1/(2*a), 1/(2*a)
			}
		{}
	};

	template <typename T, typename State>
	struct midpoint : rk2<T, State>
	// Midpoint method (2nd order, 2 stages)
	{
		midpoint() : rk2<T, State>(.5L)
		{}
	};

	template <typename T, typename State>
	struct heun2 : rk2<T, State>
	// Heun method (2nd order, 2 stages)
	{
		heun2() : rk2<T, State>(1)
		{}
	};

	template <typename T, typename State>
	struct ralston2 : rk2<T, State>
	// Ralston method (2nd order, 2 stages)
	{
		ralston2() : rk2<T, State>(2/3.L)
		{}
	};

	template <typename T, typename State>
	struct rk4 : explicit_runge_kutta_base<T, State, 4>
	// Classical Runge-Kutta 4 method (4th order, 4 stages)
	{
		rk4()
			: explicit_runge_kutta_base<T, State, 4>
			{
				.5L, .5L, 0, 0,
				.5L, 0, .5L, 0,
				1, 0, 0, 1,
				1/6.L, 1/3.L, 1/3.L, 1/6.L
			}
		{}
	};

	template <typename T, typename State>
	struct rk4_38 : explicit_runge_kutta_base<T, State, 4>
	// 3/8-rule Runge-Kutta 4 method (4th order, 4 stages)
	{
		rk4_38()
			: explicit_runge_kutta_base<T, State, 4>
			{
				1/3.L, 1/3.L, 0, 0,
				2/3.L, -1/3.L, 1, 0,
				1, 1, -1, 1,
				1/8.L, 3/8.L, 3/8.L, 1/8.L
			}
		{}
	};

	template <typename T, typename State>
	struct ralston4 : explicit_runge_kutta_base<T, State, 4>
	// Ralston method (4th order, 4 stages)
	{
		ralston4()
			: explicit_runge_kutta_base<T, State, 4>
			{
				.4L, .4L, 0, 0,
				.45573725L, .29697761L, .15875964L, 0,
				1, .2181004L, -3.05096516L, 3.83286476L,
				.17476028L, -.55148066L, 1.2055356L, .17118478L
			}
		{}
	};

	template <typename T, typename State>
	struct butcher6 : explicit_runge_kutta_base<T, State, 7>
	// Butcher method (6th order, 7 stages)
	{
		butcher6()
			: explicit_runge_kutta_base<T, State, 7>
			{
				1/3.L, 1/3.L, 0, 0, 0, 0, 0,
				2/3.L, 0, 2/3.L, 0, 0, 0, 0,
				1/3.L, 1/12.L, 1/3.L, -1/12.L, 0, 0, 0,
				.5L, -1/16.L, 9/8.L, -3/16.L, -3/8.L, 0, 0,
				.5L, 0, 9/8.L, -3/8.L, -3/4.L, .5L, 0,
				1, 9/44.L, -9/11.L, 63/44.L, 18/11.L, -16/11.L,
				11/120.L, 0, 27/40.L, 27/40.L, -4/15.L, -4/15.L, 11/120.L
			}
		{}
	};

	template <typename T, typename State>
	struct verner8 : explicit_runge_kutta_base<T, State, 11>
	// Verner method (8th order, 11 stages)
	{
		verner8()
			: explicit_runge_kutta_base<T, State, 11>
			{
				.5L, .5L, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				.5L, .25L, .25L, 0, 0, 0, 0, 0, 0, 0, 0,
				(7+s21)/14, 1/7.L, (-7-3*s21)/98, (21+5*s21)/49, 0, 0, 0, 0, 0, 0, 0,
				(7+s21)/14, (11+s21)/84, 0, (18+4*s21)/63, (21-s21)/252, 0, 0, 0, 0, 0, 0,
				5.L, (5+s21)/48, 0, (9+s21)/36, (-231+14*s21)/360, (63-7*s21)/80, 0, 0, 0, 0, 0,
				(7-s21)/14, (10-s21)/42, 0, (-432+92*s21)/315, (633-145*s21)/90, (-504+115*s21)/70, (63-13*s21)/35, 0, 0, 0, 0,
				(7-s21)/14, 1/14.L, 0, 0, 0, (14-3*s21)/126, (13-3*s21)/63, 1/9.L, 0, 0, 0,
				.5L, 1/32.L, 0, 0, 0, (91-21*s21)/576, 11/72.L, (-385-75*s21)/1152, (63+13*s21)/128, 0, 0,
				(7+s21)/14, 1/14.L, 0, 0, 0, 1/9.L, (-733-147*s21)/2205, (515+111*s21)/504, (-51-11*s21)/56, (132+28*s21)/245, 0,
				1, 0, 0, 0, 0, (-42+7*s21)/18, (-18+28*s21)/45, (-273-53*s21)/72, (301+53*s21)/72, (28-28*s21)/45, (49-7*s21)/18,
				1/20.L, 0, 0, 0, 0, 0, 0, 49/180.L, 16/45.L, 49/180.L, 1/20.L
			}
		{}

		private:

			static constexpr long double s21 = 4.582575694955840006588047193728L; // sqrt(21)
	};

	// RUNGE-KUTTA-NYSTROM SCHEMES

	// H = p^2/2m + V(q)
	// only quadratic kinetic energy!

	template <typename T, typename State, std::size_t Stages>
	struct runge_kutta_nystrom_base : runge_kutta_base<T, State>, symplectic_integrator_base<T, State>
	{
		template <typename ... Ts>
		runge_kutta_nystrom_base(Ts ... pars) : pars{T(pars)...} {}

		template <having_coordinates<T, State> S>
		void step(S& s, T dt, bool first_step = false) const
		{
			using std::size_t;

			prev_x = s.x;
			prev_v = s.vel();
			prev_p = s.p;

			k[0] = s.force(first_step);
			for (size_t i = 1; i < Stages; ++i)
			{
				s.x = 0;
				// sum smaller contributes (in dt^2 and in dt) first to minimize rounding errors
				for (size_t j = 0; j < i; ++j)
					s.x += k[j] * (pars[(i-1)*(Stages+1) + j+1] * dt * dt);
				s.x += prev_v * (pars[(i-1)*(Stages+1)] * dt);
				s.x += prev_x;

				k[i] = s.force();
			}

			s.x = 0;
			s.p = 0;
			for (size_t i = 0; i < Stages; ++i)
			{
				s.x += k[i] * (pars[(Stages-1)*Stages + i] * dt * dt);
				s.p += k[i] * (pars[Stages*Stages + i] * dt);
			}
			s.x += prev_v * dt;
			s.x += prev_x;
			s.p += prev_p;
			s.t += dt;
			s.force();
		}

		private:

			const std::array<T, (Stages+1)*Stages> pars;
			State k[Stages], prev_x, prev_v, prev_p;
	};

	// TODO

} // namespace physics

#endif // PHYSICS_INTEGRATOR_H






























