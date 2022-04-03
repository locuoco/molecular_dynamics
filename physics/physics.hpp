#ifndef PHYSICS_PHYSICS_H
#define PHYSICS_PHYSICS_H

#include "point.hpp"
#include "integrator.hpp"
#include "molecule.hpp"

#include <chrono> // steady_clock, duration_cast, microseconds
#include <random> // mt19937_64, uniform_real_distribution
#include <algorithm> // generate

inline void test()
{
	const int N = 5'000'000;
	physics::state<double, 3> x(N), y(N), z(N);
	std::valarray<double> a(3*N), b(3*N), c(3*N);
	std::mt19937_64 mersenne_twister(0);
	std::uniform_real_distribution<double> dist(-1,1);
	auto gen = [&dist, &mersenne_twister]()
	{
		return physics::point<double, 3>(dist(mersenne_twister),
										 dist(mersenne_twister),
										 dist(mersenne_twister));
	};
	auto gen2 = [&dist, &mersenne_twister]()
	{
		return dist(mersenne_twister);
	};
	std::ranges::generate(x, gen);
	std::ranges::generate(y, gen);
	std::ranges::generate(a, gen2);
	std::ranges::generate(b, gen2);

	z = x + y; // warm up

	auto begin = std::chrono::steady_clock::now();

	z = x + y;

	auto end = std::chrono::steady_clock::now();
	
	std::cout << "Time elapsed for valarray +: "
				  << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() * 1.e-6
				  << " [s]" << std::endl;

	begin = std::chrono::steady_clock::now();

	for (int i = 0; i < N; ++i)
		z[i] = x[i] + y[i];
	
	end = std::chrono::steady_clock::now();
	
	std::cout << "Time elapsed for cycle +: "
				  << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() * 1.e-6
				  << " [s]" << std::endl;

	c = a + b; // warm up

	begin = std::chrono::steady_clock::now();

	c = a + b;
	
	end = std::chrono::steady_clock::now();
	
	std::cout << "Time elapsed for valarray (scalar) +: "
				  << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() * 1.e-6
				  << " [s]" << std::endl;

	//std::cout << physics::water<float>.bonds[0][0] << '\n';
}

#endif // PHYSICS_PHYSICS_H
































