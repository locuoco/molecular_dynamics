//  Parallel sorting algorithm tests
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

#include <iostream> // cout, endl
#include <vector>
#include <algorithm> // sort, generate
#include <functional> // greater
#include <random>
#include <chrono>
#include <cassert>

/*

g++ parallel_sort.cpp -o parallel_sort -std=c++20 -Wall -Wextra -pedantic -Ofast -pthread -fmax-errors=1

*/

#include "../thread_pool.hpp"
#include "../parallel_sort.hpp"

utils::thread_pool tp;
std::mt19937 mersenne_twister;

template <typename T>
std::uniform_real_distribution<T> dist(-1, 1);

void test_parallel_sort_correct()
{
	for (std::size_t n = 16; n <= (1 << 20); n <<= 1)
	{
		std::vector<double> x(n), xref;

		std::generate(begin(x), end(x), [&]{ return dist<double>(mersenne_twister); });
		xref = x;

		utils::parallel_sort(begin(x), end(x), tp);
		std::sort(begin(xref), end(xref));

		assert(x == xref);
	}
}

void test_parallel_sort_reverse_correct()
{
	tp.resize(6);
	for (std::size_t n = 21; n <= (1 << 20); n <<= 1)
	{
		std::vector<double> x(n), xref;

		std::generate(begin(x), end(x), [&]{ return dist<double>(mersenne_twister); });
		xref = x;

		utils::parallel_sort(begin(x), end(x), tp, std::greater{});
		std::sort(begin(xref), end(xref), std::greater{});

		assert(x == xref);
	}
}

void test_parallel_sort_perf()
{
	std::cout << " ======   utils::parallel_sort   ====== " << std::endl;
	std::size_t n_loops = 10;
	for (std::size_t n = 128; n <= (1 << 20); n <<= 1)
	{
		std::vector<double> x(n);

		decltype(std::chrono::steady_clock::now()) start, finish;
		double timing = 0;

		for (size_t i = 0; i < n_loops+1; ++i)
		{
			std::generate(begin(x), end(x), [&]{ return dist<double>(mersenne_twister); });

			start = std::chrono::steady_clock::now();
			utils::parallel_sort(begin(x), end(x), tp);
			finish = std::chrono::steady_clock::now();

			if (i)
				timing += std::chrono::duration<double>(finish-start).count();
		}
		timing /= n_loops;
		std::cout << "n = " << n << " -- t = " << timing << "s -- x_0 = " << x[0] << std::endl;
	}
}

void test_sort_perf()
{
	std::cout << " ======   std::sort   ====== " << std::endl;
	std::size_t n_loops = 10;
	for (std::size_t n = 128; n <= (1 << 20); n <<= 1)
	{
		std::vector<double> x(n);

		decltype(std::chrono::steady_clock::now()) start, finish;
		double timing = 0;

		for (std::size_t i = 0; i < n_loops+1; ++i)
		{
			std::generate(begin(x), end(x), [&]{ return dist<double>(mersenne_twister); });

			start = std::chrono::steady_clock::now();
			std::sort(begin(x), end(x));
			finish = std::chrono::steady_clock::now();

			if (i)
				timing += std::chrono::duration<double>(finish-start).count();
		}
		timing /= n_loops;
		std::cout << "n = " << n << " -- t = " << timing << "s -- x_0 = " << x[0] << std::endl;
	}
}

int main()
{
	test_parallel_sort_correct();
	test_parallel_sort_reverse_correct();
	test_parallel_sort_perf();
	test_sort_perf();

	return 0;
}
























