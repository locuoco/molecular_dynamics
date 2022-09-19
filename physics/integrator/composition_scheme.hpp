//  Composition schemes
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

#ifndef PHYSICS_INTEGRATOR_COMPOSITION_SCHEME_H
#define PHYSICS_INTEGRATOR_COMPOSITION_SCHEME_H

#include "integrator.hpp"

// COMPOSITION SCHEMES
// they can be used as long as the 2nd order integrator is known.
// Useful to build higher-order symplectic methods.
// Composing the 2nd-order isokinetic integrator should give higher-order isokinetic integrators

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

namespace physics
{
	template <having_coordinates System, template <typename> typename IntegT, std::size_t Stages>
	struct composition_scheme_base : virtual IntegT<System>
	// all composition schemes inherit from this class.
	// `IntegT` is the 2nd order integrator which is going to be composed.
	// `Stages` is the number of stages of the method (number of force evaluations).
	{
		template <typename ... Ts>
		requires (sizeof...(Ts) == (Stages+1)/2)
		composition_scheme_base(Ts ... pars) : d{scalar_type_of<System>(pars)...} {}
		// constructor: `pars...` is a variadic argument which contains
		// the parameters for the composition scheme

		virtual ~composition_scheme_base() = default;

		void step(System& s, scalar_type_of<System> dt) override
		{
			for (size_t i = 0; i < Stages; ++i)
				IntegT<System>::step(s, d[std::min(i, Stages-i-1)] * dt);
		}

		private:

			const std::array<scalar_type_of<System>, (Stages+1)/2> d;
	};

	template <typename Integ, typename System, template <typename> typename IntegT, std::size_t Stages>
	concept composition_scheme = std::is_base_of_v<composition_scheme_base<System, IntegT, Stages>, Integ>;

	template <typename System, template <typename> typename IntegT = leapfrog>
	struct forest_ruth : composition_scheme_base<System, IntegT, 3>
	// Forest-Ruth method (4th order, 3 stages)
	// independently obtained by Yoshida
	{
		forest_ruth()
			: composition_scheme_base<System, IntegT, 3>
			{
				1.3512071919596576340476878089715L, // 1 / (2 - cbrt(2))
				-1.7024143839193152680953756179429L // -cbrt(2) / (2 - cbrt(2))
			}
		{}
	};

	template <typename System, template <typename> typename IntegT = leapfrog>
	struct suzuki4 : composition_scheme_base<System, IntegT, 5>
	// Suzuki method (4th order, 5 stages)
	{
		suzuki4()
			: composition_scheme_base<System, IntegT, 5>
			{
				0.41449077179437573714235406286076L, // 1 / (4 - cbrt(4))
				0.41449077179437573714235406286076L, // 1 / (4 - cbrt(4))
				-0.65796308717750294856941625144305L // -cbrt(4) / (4 - cbrt(4))
			}
		{}
	};

	template <typename System, template <typename> typename IntegT = leapfrog>
	struct kahan_li4a : composition_scheme_base<System, IntegT, 5>
	// Kahan-Li method (4th order, 5 stages)
	{
		kahan_li4a()
			: composition_scheme_base<System, IntegT, 5>
			{
				0.78867513459481288225457439025098L, // (3 + sqrt(3)) / 6
				0.21132486540518711774542560974902L, // (3 - sqrt(3)) / 6
				-1
			}
		{}
	};

	template <typename System, template <typename> typename IntegT = leapfrog>
	struct kahan_li4b : composition_scheme_base<System, IntegT, 5>
	// Kahan-Li method (4th order, 5 stages)
	{
		kahan_li4b()
			: composition_scheme_base<System, IntegT, 5>
			{
				0.21132486540518711774542560974902L, // (3 - sqrt(3)) / 6
				0.78867513459481288225457439025098L, // (3 + sqrt(3)) / 6
				-1
			}
		{}
	};

	template <typename System, template <typename> typename IntegT = leapfrog>
	struct yoshida6 : composition_scheme_base<System, IntegT, 7>
	// Yoshida method (6th order, 7 stages)
	{
		yoshida6()
			: composition_scheme_base<System, IntegT, 7>
			{
				0.78451361047755726382L,
				0.23557321335935813368L,
				-1.1776799841788710069L,
				1.3151863206839112189L
			}
		{}
	};

	template <typename System, template <typename> typename IntegT = leapfrog>
	struct kahan_li6a : composition_scheme_base<System, IntegT, 9>
	// Kahan-Li method (6th order, 9 stages)
	{
		kahan_li6a()
			: composition_scheme_base<System, IntegT, 9>
			{
				0.39216144400731413928L,
				0.33259913678935943860L,
				-0.70624617255763935981L,
				0.082213596293550800230L,
				0.79854399093482996340L
			}
		{}
	};

	template <typename System, template <typename> typename IntegT = leapfrog>
	struct kahan_li6b : composition_scheme_base<System, IntegT, 9>
	// Kahan-Li method (6th order, 9 stages)
	{
		kahan_li6b()
			: composition_scheme_base<System, IntegT, 9>
			{
				0.39103020330868478817L,
				0.33403728961113601749L,
				-0.70622728118756134346L,
				0.081877549648059445768L,
				0.79856447723936218406L
			}
		{}
	};

	template <typename System, template <typename> typename IntegT = leapfrog>
	struct suzuki8 : composition_scheme_base<System, IntegT, 15>
	// Suzuki method (8th order, 15 stages)
	{
		suzuki8()
			: composition_scheme_base<System, IntegT, 15>
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

	template <typename System, template <typename> typename IntegT = leapfrog>
	struct kahan_li8a : composition_scheme_base<System, IntegT, 17>
	// Kahan-Li method (8th order, 17 stages)
	{
		kahan_li8a()
			: composition_scheme_base<System, IntegT, 17>
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

	template <typename System, template <typename> typename IntegT = leapfrog>
	struct kahan_li8b : composition_scheme_base<System, IntegT, 17>
	// Kahan-Li method (8th order, 17 stages)
	{
		kahan_li8b()
			: composition_scheme_base<System, IntegT, 17>
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

	template <typename System, template <typename> typename IntegT = leapfrog>
	struct kahan_li10a : composition_scheme_base<System, IntegT, 31>
	// Kahan-Li method (10th order, 31 stages)
	{
		kahan_li10a()
			: composition_scheme_base<System, IntegT, 31>
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

	template <typename System, template <typename> typename IntegT = leapfrog>
	struct kahan_li10b : composition_scheme_base<System, IntegT, 31>
	// Kahan-Li method (10th order, 31 stages)
	{
		kahan_li10b()
			: composition_scheme_base<System, IntegT, 31>
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

	template <typename System, template <typename> typename IntegT = leapfrog>
	struct kahan_li10c : composition_scheme_base<System, IntegT, 33>
	// Kahan-Li method (10th order, 33 stages)
	{
		kahan_li10c()
			: composition_scheme_base<System, IntegT, 33>
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

	template <typename System, template <typename> typename IntegT = leapfrog>
	struct kahan_li10d : composition_scheme_base<System, IntegT, 33>
	// Kahan-Li method (10th order, 33 stages)
	{
		kahan_li10d()
			: composition_scheme_base<System, IntegT, 33>
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

} // namespace physics

#endif // PHYSICS_INTEGRATOR_COMPOSITION_SCHEME_H






























