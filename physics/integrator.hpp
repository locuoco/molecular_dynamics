#ifndef PHYSICS_INTEGRATOR_H
#define PHYSICS_INTEGRATOR_H

#include "point.hpp"

#include <concepts> // convertible_to
#include <valarray>
#include <type_traits> // is_base_of
#include <cmath> // exp, sqrt

namespace physics
{
	template <typename T>
	struct physical_system_base
	{
		virtual void accel() = 0;
	};

	template <typename S, typename T>
	concept physical_system = std::is_base_of_v<physical_system_base<T>, S>;

	template <typename S, typename T>
	concept having_coordinates = physical_system<S, T>
		&& requires(S& s, T dt)
		{
			s.x += s.v * dt;
			s.v += s.a * dt;
			s.t += dt;
			s.accel();
		};

	template <typename S, typename T>
	concept having_coordinates_damped = having_coordinates<S, T>
		&& requires(S& s, T dt, std::size_t i)
		{
			{s.gamma[i]} -> std::convertible_to<T>;
			s.v[i] = s.v[i] * dt;
			i < s.n;
		};

	template <typename S, typename T>
	concept having_coordinates_stochastic = having_coordinates_damped<S, T>
		&& requires(S& s, T dt, std::size_t i)
		{
			{s.D} -> std::convertible_to<T>;
			s.v[i] = s.v[i] * dt + s.noise[i] * dt;
			s.rand();
		};

	struct integrator_base {};

	struct symplectic_integrator_base : integrator_base {};
	// Assumption: H = T(p) + V(x) = p^2/2m + V(x), with a = -(1/m) grad V(x) and v = p/m
	// Furthermore, H is explicitly independent on time
	// unless stated otherwise

	struct stochastic_integrator_base : integrator_base {};

	template <typename Integ>
	concept integrator = std::is_base_of_v<integrator_base, Integ>;

	template <typename Integ>
	concept symplectic_integrator = std::is_base_of_v<symplectic_integrator_base, Integ>;

	template <typename Integ>
	concept stochastic_integrator = std::is_base_of_v<stochastic_integrator_base, Integ>;

	template <typename T>
	struct symplectic_euler : symplectic_integrator_base
	// Symplectic Euler method (1st order, 1 stage)
	// It is equivalent to leapfrog for long-term simulations, despite being lower order
	{
		template <having_coordinates<T> S>
		void step(S& s, T dt) const
		{
			s.accel();
			s.v += s.a * dt;
			s.x += s.v * dt;

			s.t += dt;
		}
	};

	template <typename T>
	struct leapfrog : symplectic_integrator_base
	// Leapfrog method (2nd order, 1 stage)
	{
		template <having_coordinates<T> S>
		void step(S& s, T dt) const
		{
			s.x += s.v * (dt/2);
			s.accel();
			s.v += s.a * dt;
			s.x += s.v * (dt/2);

			s.t += dt;
		}
	};

	template <typename T>
	struct stochastic_leapfrog : symplectic_integrator_base, stochastic_integrator_base
	// stochastic leapfrog (2nd order?, 1 stage)
	// dp	= f dt - gamma p dt + sigma dw
	// dv	= a dt - gamma v dt + (sigma/m) dw
	//		= a dt + (-chi v dt + sigma dw) / m
	// with a = f/m, v = p/m,
	// 		gamma = k T / m D [related to diffusion coefficient],
	// 		sigma = sqrt(2 gamma m k T) = k T sqrt(2 / D) [fluctuation-dissipation theorem],
	//		dw : Wiener process differential
	//		chi = m * gamma = k T / D  ===>   sigma = chi * sqrt(2 D)
	// H = H0(x, p) + H1(p) xi(t) <- stochastic hamiltonian (xi : gaussian process)
/*
	ALLEN, TILDESLEY
	COMPUTER SIMULATION OF LIQUIDS
	2017
	P. 383-389
*/
	{
		template <having_coordinates_stochastic<T> S>
		void step(S& s, T dt) const
		{
			using std::exp;
			using std::sqrt;

			s.v += s.a * (dt/2);
			s.x += s.v * (dt/2);
			s.rand();
			if (s.D == 0) // gamma -> inf
				s.v = s.noise; // s.noise = sqrt(k T / m) * gaussian noise
			else
				for (std::size_t i = 0; i < s.n; ++i)
					s.v[i] = exp(-s.gamma[i] * dt) * s.v[i] + sqrt(1 - exp(-2 * s.gamma[i] * dt)) * s.noise[i];
				// v = exp(-gamma dt) v + sqrt((1 - exp(-2 gamma dt)) k T / m) g
				// expanding the exp (gamma small): v += -gamma v dt + sqrt(2 gamma k T dt / m) g
				//									v += -(chi v dt + sigma sqrt(dt) g) / m
				// g: gaussian noise with mean = 0, std.dev = 1
			s.x += s.v * (dt/2);
			s.accel();
			s.v += s.a * (dt/2);

			s.t += dt;
		}
	};

	template <typename T>
	struct damped_leapfrog : integrator_base
	// damped integrator (2nd order?, 1 stage)
	// dx = dp/m, dp = f dt - gamma p dt
	//			  dv = a dt - gamma v dt
	// with a = f/m, v = p/m,
	{
		template <having_coordinates_damped<T> S>
		void step(S& s, T dt) const
		{
			using std::exp;

			s.v += s.a * (dt/2);
			s.x += s.v * (dt/2);
			for (std::size_t i = 0; i < s.n; ++i)
				s.v *= exp(-s.gamma[i] * dt);
			s.x += s.v * (dt/2);
			s.accel();
			s.v += s.a * (dt/2);

			s.t += dt;
		}
	};

	template <typename T>
	struct pefrl : symplectic_integrator_base
	// Position-extended Forest-Ruth-like (4th order, 4 stages)
	// OMELYAN, MRYGLOD, FOLK
	// OPTIMIZED FOREST-RUTH- AND SUZUKI-LIKE ALGORITHMS FOR INTEGRATION
	// OF MOTION IN MANY-BODY SYSTEMS
	// 2008
	{
		template <having_coordinates<T> S>
		void step(S& s, T dt) const
		{
			s.x += s.v * (xi * dt);
			s.accel();
			s.v += s.a * ((1 - 2*lambda) * dt/2);
			s.x += s.v * (chi * dt);
			s.accel();
			s.v += s.a * (lambda * dt);
			s.x += s.v * ((1 - 2*(chi + xi)) * dt);
			s.accel();
			s.v += s.a * (lambda * dt);
			s.x += s.v * (chi * dt);
			s.accel();
			s.v += s.a * ((1 - 2*lambda) * dt/2);
			s.x += s.v * (xi * dt);

			s.t += dt;
		}

		private:

			static constexpr T xi = 0.1786178958448091L;
			static constexpr T lambda = -0.2123418310626054L;
			static constexpr T chi = -0.6626458266981849e-1L;
	};

	// COMPOSITION SCHEMES
	// can be used as long as the symplectic integrator of 2nd order for a given hamiltonian is known

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

	template <typename T, symplectic_integrator Integ, std::size_t Stages>
	struct composition_scheme_base : Integ
	{
		template <typename ... Ts>
		composition_scheme_base(Ts ... pars) : d{T(pars)...} {}

		template <having_coordinates<T> S>
		void step(S& s, T dt) const
		{
			for (size_t i = 0; i < Stages; ++i)
				Integ::step(s, d[std::min(i, Stages-i-1)] * dt);
		}

		private:

			const std::array<T, (Stages+1)/2> d;
	};

	template <typename Integ, typename T, typename Integ2, std::size_t Stages>
	concept composition_scheme = std::is_base_of_v<composition_scheme_base<T, Integ2, Stages>, Integ>;

	template <typename T, symplectic_integrator Integ = leapfrog<T>>
	struct forest_ruth : composition_scheme_base<T, Integ, 3>
	// Forest-Ruth method (4th order, 3 stages)
	// independently obtained by Yoshida
	{
		forest_ruth()
			: composition_scheme_base<T, Integ, 3>{1.3512071919596576340476878089715L, // 1 / (2 - cbrt(2))
												   -1.7024143839193152680953756179429L} // -cbrt(2) / (2 - cbrt(2))
		{}
	};

	template <typename T, symplectic_integrator Integ = leapfrog<T>>
	struct suzuki4 : composition_scheme_base<T, Integ, 5>
	// Suzuki method (4th order, 5 stages)
	{
		suzuki4()
			: composition_scheme_base<T, Integ, 5>{0.41449077179437573714235406286076L, // 1 / (4 - cbrt(4))
												   0.41449077179437573714235406286076L, // 1 / (4 - cbrt(4))
												   -0.65796308717750294856941625144305L} // -cbrt(4) / (4 - cbrt(4))
		{}
	};

	template <typename T, symplectic_integrator Integ = leapfrog<T>>
	struct kahan_li4a : composition_scheme_base<T, Integ, 5>
	// Kahan-Li method (4th order, 5 stages)
	{
		kahan_li4a()
			: composition_scheme_base<T, Integ, 5>{0.78867513459481288225457439025098L, // (3 + sqrt(3)) / 6
												   0.21132486540518711774542560974902L, // (3 - sqrt(3)) / 6
												   -1}
		{}
	};

	template <typename T, symplectic_integrator Integ = leapfrog<T>>
	struct kahan_li4b : composition_scheme_base<T, Integ, 5>
	// Kahan-Li method (4th order, 5 stages)
	{
		kahan_li4b()
			: composition_scheme_base<T, Integ, 5>{0.21132486540518711774542560974902L, // (3 - sqrt(3)) / 6
												   0.78867513459481288225457439025098L, // (3 + sqrt(3)) / 6
												   -1}
		{}
	};

	template <typename T, symplectic_integrator Integ = leapfrog<T>>
	struct yoshida6 : composition_scheme_base<T, Integ, 7>
	// Yoshida method (6th order, 7 stages)
	{
		yoshida6()
			: composition_scheme_base<T, Integ, 7>{0.78451361047755726382L,
												   0.23557321335935813368L,
												   -1.1776799841788710069L,
												   1.3151863206839112189L}
		{}
	};

	template <typename T, symplectic_integrator Integ = leapfrog<T>>
	struct kahan_li6a : composition_scheme_base<T, Integ, 9>
	// Kahan-Li method (6th order, 9 stages)
	{
		kahan_li6a()
			: composition_scheme_base<T, Integ, 9>{0.39216144400731413928L,
												   0.33259913678935943860L,
												   -0.70624617255763935981L,
												   0.082213596293550800230L,
												   0.79854399093482996340L}
		{}
	};

	template <typename T, symplectic_integrator Integ = leapfrog<T>>
	struct kahan_li6b : composition_scheme_base<T, Integ, 9>
	// Kahan-Li method (6th order, 9 stages)
	{
		kahan_li6b()
			: composition_scheme_base<T, Integ, 9>{0.39103020330868478817L,
												   0.33403728961113601749L,
												   -0.70622728118756134346L,
												   0.081877549648059445768L,
												   0.79856447723936218406L}
		{}
	};

	template <typename T, symplectic_integrator Integ = leapfrog<T>>
	struct suzuki8 : composition_scheme_base<T, Integ, 15>
	// Suzuki method (8th order, 15 stages)
	{
		suzuki8()
			: composition_scheme_base<T, Integ, 15>{0.74167036435061295344822780L,
													-0.40910082580003159399730010L,
													0.19075471029623837995387626L,
													-0.57386247111608226665638773L,
													0.29906418130365592384446354L,
													0.33462491824529818378495798L,
													0.31529309239676659663205666L,
													-0.79688793935291635401978884L}
		{}
	};

	template <typename T, symplectic_integrator Integ = leapfrog<T>>
	struct kahan_li8a : composition_scheme_base<T, Integ, 17>
	// Kahan-Li method (8th order, 17 stages)
	{
		kahan_li8a()
			: composition_scheme_base<T, Integ, 17>{0.13020248308889008088L,
													0.56116298177510838456L,
													-0.38947496264484728641L,
													0.15884190655515560090L,
													-0.39590389413323757734L,
													0.18453964097831570709L,
													0.25837438768632204729L,
													0.29501172360931029887L,
													-0.60550853383003451170L}
		{}
	};

	template <typename T, symplectic_integrator Integ = leapfrog<T>>
	struct kahan_li8b : composition_scheme_base<T, Integ, 17>
	// Kahan-Li method (8th order, 17 stages)
	{
		kahan_li8b()
			: composition_scheme_base<T, Integ, 17>{0.12713692773487857916L,
													0.56170253798880269972L,
													-0.38253471994883018888L,
													0.16007605629464743119L,
													-0.40181637432680696673L,
													0.18736671654227849724L,
													0.26070870920779240570L,
													0.29039738812516162389L,
													-0.60607448323584816258L}
		{}
	};

	template <typename T, symplectic_integrator Integ = leapfrog<T>>
	struct kahan_li10a : composition_scheme_base<T, Integ, 31>
	// Kahan-Li method (10th order, 31 stages)
	{
		kahan_li10a()
			: composition_scheme_base<T, Integ, 31>{-0.48159895600253002870L,
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
													0.91249037635867994571L}
		{}
	};

	template <typename T, symplectic_integrator Integ = leapfrog<T>>
	struct kahan_li10b : composition_scheme_base<T, Integ, 31>
	// Kahan-Li method (10th order, 31 stages)
	{
		kahan_li10b()
			: composition_scheme_base<T, Integ, 31>{0.27338476926228452782L,
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
													0.28614352562198582747L}
		{}
	};

	template <typename T, symplectic_integrator Integ = leapfrog<T>>
	struct kahan_li10c : composition_scheme_base<T, Integ, 33>
	// Kahan-Li method (10th order, 33 stages)
	{
		kahan_li10c()
			: composition_scheme_base<T, Integ, 33>{0.070428877682658066880L,
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
													0.95529818470370207691L}
		{}
	};

	template <typename T, symplectic_integrator Integ = leapfrog<T>>
	struct kahan_li10d : composition_scheme_base<T, Integ, 33>
	// Kahan-Li method (10th order, 33 stages)
	{
		kahan_li10d()
			: composition_scheme_base<T, Integ, 33>{0.12282427644721572094L,
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
													0.83629272067135846284L}
		{}
	};

	// RUNGE-KUTTA SCHEMES

	struct runge_kutta_base {};

	template <typename T, std::size_t Stages>
	struct explicit_runge_kutta_base : runge_kutta_base
	{
		template <typename ... Ts>
		explicit_runge_kutta_base(Ts ... pars) : pars{T(pars)...} {}

		template <having_coordinates<T> S>
		void step(S& s, T dt) const
		{
			using std::size_t;

			prev_x = s.x;
			prev_v = s.v;
			prev_t = s.t;

			kx[0] = s.v;
			s.accel();
			kv[0] = s.a;
			for (size_t i = 1; i < Stages; ++i)
			{
				s.x = 0;
				s.v = 0;
				// sum smaller contributes (in dt) first to minimize rounding errors
				for (size_t j = 0; j < i; ++j)
				{
					T mult = pars[(i-1)*Stages + j+1] * dt;
					s.x += kx[j] * mult;
					s.v += kv[j] * mult;
				}
				s.x += prev_x;
				s.v += prev_v;
				s.t = prev_t + pars[(i-1)*Stages] * dt;

				kx[i] = s.v;
				s.accel();
				kv[i] = s.a;
			}

			s.x = 0;
			s.v = 0;
			for (size_t i = 0; i < Stages; ++i)
			{
				T mult = pars[(Stages-1)*Stages + i] * dt;
				s.x += kx[i] * mult;
				s.v += kv[i] * mult;
			}
			s.x += prev_x;
			s.v += prev_v;
			s.t = prev_t + dt;
		}

		private:

			const std::array<T, Stages*Stages> pars;
			std::valarray<T> kx[Stages], kv[Stages], prev_x, prev_v;
			T prev_t;
	};

	template <typename T>
	struct euler : explicit_runge_kutta_base<T, 1>
	// Euler method (1st order, 1 stage)
	{
		euler()
			: explicit_runge_kutta_base<T, 1>{1}
		{}
	};

	template <typename T>
	struct rk2 : explicit_runge_kutta_base<T, 2>
	// Parametrized Runge-Kutta 2 method (2nd order, 2 stages)
	{
		rk2(long double a)
			: explicit_runge_kutta_base<T, 2>{a, a,
											  1 - 1/(2*a), 1/(2*a)}
		{}
	};

	template <typename T>
	struct midpoint : rk2<T>
	// Midpoint method (2nd order, 2 stages)
	{
		midpoint()
			: rk2<T>(.5L)
		{}
	};

	template <typename T>
	struct heun2 : rk2<T>
	// Heun method (2nd order, 2 stages)
	{
		heun2()
			: rk2<T>(1)
		{}
	};

	template <typename T>
	struct ralston2 : rk2<T>
	// Ralston method (2nd order, 2 stages)
	{
		ralston2()
			: rk2<T>(2/3.L)
		{}
	};

	template <typename T>
	struct rk4 : explicit_runge_kutta_base<T, 4>
	// Classical Runge-Kutta 4 method (4th order, 4 stages)
	{
		rk4()
			: explicit_runge_kutta_base<T, 4>{.5L, .5L, 0, 0,
											  .5L, 0, .5L, 0,
											  1, 0, 0, 1,
											  1/6.L, 1/3.L, 1/3.L, 1/6.L}
		{}
	};

	template <typename T>
	struct rk4_38 : explicit_runge_kutta_base<T, 4>
	// 3/8-rule Runge-Kutta 4 method (4th order, 4 stages)
	{
		rk4_38()
			: explicit_runge_kutta_base<T, 4>{1/3.L, 1/3.L, 0, 0,
											  2/3.L, -1/3.L, 1, 0,
											  1, 1, -1, 1,
											  1/8.L, 3/8.L, 3/8.L, 1/8.L}
		{}
	};

	template <typename T>
	struct ralston4 : explicit_runge_kutta_base<T, 4>
	// Ralston method (4th order, 4 stages)
	{
		ralston4()
			: explicit_runge_kutta_base<T, 4>{.4L, .4L, 0, 0,
											  .45573725L, .29697761L, .15875964L, 0,
											  1, .2181004L, -3.05096516L, 3.83286476L,
											  .17476028L, -.55148066L, 1.2055356L, .17118478L}
		{}
	};

	template <typename T>
	struct butcher6 : explicit_runge_kutta_base<T, 7>
	// Butcher method (6th order, 7 stages)
	{
		butcher6()
			: explicit_runge_kutta_base<T, 7>{1/3.L, 1/3.L, 0, 0, 0, 0, 0,
											  2/3.L, 0, 2/3.L, 0, 0, 0, 0,
											  1/3.L, 1/12.L, 1/3.L, -1/12.L, 0, 0, 0,
											  .5L, -1/16.L, 9/8.L, -3/16.L, -3/8.L, 0, 0,
											  .5L, 0, 9/8.L, -3/8.L, -3/4.L, .5L, 0,
											  1, 9/44.L, -9/11.L, 63/44.L, 18/11.L, -16/11.L,
											  11/120.L, 0, 27/40.L, 27/40.L, -4/15.L, -4/15.L, 11/120.L}
		{}
	};

	template <typename T>
	struct verner8 : explicit_runge_kutta_base<T, 11>
	// Verner method (8th order, 11 stages)
	{
		verner8()
		: explicit_runge_kutta_base<T, 11>{.5L, .5L, 0, 0, 0, 0, 0, 0, 0, 0, 0,
										   .5L, .25L, .25L, 0, 0, 0, 0, 0, 0, 0, 0,
										   (7+s21)/14, 1/7.L, (-7-3*s21)/98, (21+5*s21)/49, 0, 0, 0, 0, 0, 0, 0,
										   (7+s21)/14, (11+s21)/84, 0, (18+4*s21)/63, (21-s21)/252, 0, 0, 0, 0, 0, 0,
										   5.L, (5+s21)/48, 0, (9+s21)/36, (-231+14*s21)/360, (63-7*s21)/80, 0, 0, 0, 0, 0,
										   (7-s21)/14, (10-s21)/42, 0, (-432+92*s21)/315, (633-145*s21)/90, (-504+115*s21)/70, (63-13*s21)/35, 0, 0, 0, 0,
										   (7-s21)/14, 1/14.L, 0, 0, 0, (14-3*s21)/126, (13-3*s21)/63, 1/9.L, 0, 0, 0,
										   .5L, 1/32.L, 0, 0, 0, (91-21*s21)/576, 11/72.L, (-385-75*s21)/1152, (63+13*s21)/128, 0, 0,
										   (7+s21)/14, 1/14.L, 0, 0, 0, 1/9.L, (-733-147*s21)/2205, (515+111*s21)/504, (-51-11*s21)/56, (132+28*s21)/245, 0,
										   1, 0, 0, 0, 0, (-42+7*s21)/18, (-18+28*s21)/45, (-273-53*s21)/72, (301+53*s21)/72, (28-28*s21)/45, (49-7*s21)/18,
										   1/20.L, 0, 0, 0, 0, 0, 0, 49/180.L, 16/45.L, 49/180.L, 1/20.L}
		{}

		private:

			static constexpr long double s21 = 4.582575694955840006588047193728L; // sqrt(21)
	};

	// RUNGE-KUTTA-NYSTROM SCHEMES

	template <typename T, std::size_t Stages>
	struct runge_kutta_nystrom_base : runge_kutta_base, symplectic_integrator_base
	{
		template <typename ... Ts>
		runge_kutta_nystrom_base(Ts ... pars) : pars{T(pars)...} {}

		template <having_coordinates<T> S>
		void step(S& s, T dt) const
		{
			using std::size_t;

			prev_x = s.x;
			prev_v = s.v;

			s.accel();
			k[0] = s.a;
			for (size_t i = 1; i < Stages; ++i)
			{
				s.x = 0;
				// sum smaller contributes (in dt^2 and in dt) first to minimize rounding errors
				for (size_t j = 0; j < i; ++j)
					s.x += k[j] * (pars[(i-1)*(Stages+1) + j+1] * dt * dt);
				s.x += prev_v * (pars[(i-1)*(Stages+1)] * dt);
				s.x += prev_x;

				s.accel();
				k[i] = s.a;
			}

			s.x = 0;
			s.v = 0;
			for (size_t i = 0; i < Stages; ++i)
			{
				s.x += k[i] * (pars[(Stages-1)*Stages + i] * dt * dt);
				s.v += k[i] * (pars[Stages*Stages + i] * dt);
			}
			s.x += prev_v * dt;
			s.x += prev_x;
			s.v += prev_v;
			s.t += dt;
		}

		private:

			const std::array<T, (Stages+1)*Stages> pars;
			std::valarray<T> k[Stages], prev_x, prev_v;
	};

} // namespace physics

#endif // PHYSICS_INTEGRATOR_H






























