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

Compilation:
g++ parallel_sort.cpp -o parallel_sort -std=c++20 -Wall -Wextra -pedantic -Ofast -pthread -fmax-errors=1

*/

#include "../../utils/thread_pool.hpp"
#include "../../utils/parallel_sort.hpp"

utils::thread_pool tp;
std::mt19937_64 mersenne_twister;
// mersenne_twister will behave in the same way for all compilers/runs

auto gen()
// generate a pseudo-random number in [-1, 1]
{
	double two63 = 1ull << 63; // 2 ** 63
	return mersenne_twister()/two63-1;
}

void test_gen()
// test that the values generated with `gen` are between -1 and 1
{
	mersenne_twister.seed(1234);
	std::valarray<double> x(1000);

	std::ranges::generate(x, gen);

	assert((std::abs(x) <= 1).min());
}

void test_parallel_sort_correct(std::size_t n = 1 << 20)
// oracle test between `utils::parallel_sort` and `std::sort`.
// `n` is the number of elements of the sequence used in the test.
{
	mersenne_twister.seed(1234);
	std::vector<double> x(n), xref;

	std::ranges::generate(x, gen);
	xref = x;

	utils::parallel_sort(begin(x), end(x), tp);
	std::ranges::sort(xref);

	assert(x == xref);
}

void test_parallel_sort_reverse_correct(std::size_t n = 1 << 20)
// oracle test between `utils::parallel_sort` and `std::sort` (reversed/decreasing order).
// `n` is the number of elements of the sequence used in the test.
{
	mersenne_twister.seed(1234);
	std::vector<double> x(n), xref;

	std::ranges::generate(x, gen);
	xref = x;

	utils::parallel_sort(begin(x), end(x), tp, std::greater{});
	std::ranges::sort(xref, std::greater{});

	assert(x == xref);
}

int main()
{
	test_gen();
	test_parallel_sort_correct();
	test_parallel_sort_reverse_correct();

	std::cout << "All tests passed successfully!" << std::endl;

	return 0;
}
























