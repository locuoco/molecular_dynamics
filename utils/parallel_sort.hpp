//  Parallel sorting algorithm based on std::sort and std::inplace_merge
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

#ifndef UTILS_PARALLEL_SORT_H
#define UTILS_PARALLEL_SORT_H

#include <iostream> // clog
#include <algorithm> // sort, inplace_merge, min
#include <iterator> // random_access_iterator
#include <functional> // less

#include "thread_pool.hpp"
#include "../math/helper.hpp"

namespace utils
{
	template <std::random_access_iterator It, typename Compare = std::less<>>
	void parallel_sort(It first, It last, thread_pool& tp, Compare comp = std::less{})
	// PARALLEL SORT
	// it uses a thread pool `tp` to sort many parts of the range [first, last) concurrently, and
	// then merge them using an inplace algorithm from the C++ standard library (std::inplace_merge).
	// `comp` is the comparator with respect to which the range is sorted, which by default
	// is `std::less{}` (ascending order).
	// Theoretical complexity: O(n + (n/k) log n) where k is the number of processors
	{
		using std::size_t;
		if (last - first <= 0)
			return; // nothing to sort
		size_t n_blocks = tp.size();
		// if the number of threads is not a power of 2, we increase it to the nearest
		// power of 2 (hopefully without increasing overhead).
		// Having a power-of-2 number of threads makes the algorithm easier to implement.
		if (!math::is_pow_of_2(n_blocks))
		{
			do ++n_blocks;
			while (!math::is_pow_of_2(n_blocks));
			std::clog << "Warning: Number of threads in thread pool is not a power of 2. "
				"Increasing to the nearest power of 2...: " << n_blocks << '\n';
			tp.resize(n_blocks);
		}
		size_t n = last - first;
		size_t block_size = (n-1)/n_blocks+1;
		// divide the range into `n_blocks` blocks and sort them concurrently (using the thread pool).
		// This part has O((n/k) log n) complexity where k is the number of processors
		for (size_t i = 0; i < n_blocks; ++i)
		{
			size_t start = i*block_size;
			size_t finish = std::min(n, (i+1)*block_size);
			tp.enqueue(std::sort<It, Compare>, first + start, first + finish, comp);
		}
		tp.wait();

		// merge the sorted parts hierarchically with `std::inplace_merge`, which is an optimized algorithm with
		// linear complexity which sorts two sorted subranges into one (it is the basic ingredient of merge-sort).
		// This part has O(n) complexity.
		for (size_t div = n_blocks/2; div > 1; block_size *= 2, div /= 2)
		{
			for (size_t i = 0; i < div; ++i)
			{
				size_t start = 2*i*block_size;
				size_t middle = (2*i+1)*block_size;
				size_t finish = std::min(n, 2*(i+1)*block_size);
				tp.enqueue(std::inplace_merge<It, Compare>, first + start, first + middle, first + finish, comp);
			}
			tp.wait();
		}
		// last iteration of in-place merge is not done concurrently.
		// This part has O(n) complexity
		size_t start = 0;
		size_t middle = block_size;
		size_t finish = std::min(n, 2*block_size);
		std::inplace_merge(first + start, first + middle, first + finish, comp);
	}

} // namespace utils

#endif // UTILS_PARALLEL_SORT_H

























