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

#include <algorithm> // sort, inplace_merge, min
#include <functional> // less
#include <iterator> // random_access_iterator

#include "thread_pool.hpp"

namespace utils
{
	template <std::random_access_iterator It, typename ... Others>
	void parallel_sort(It first, It last, thread_pool& tp, Others ... others)
	{
		using std::size_t;
		if (last - first <= 0)
			return;
		size_t num_threads = tp.size();
		if ((num_threads&(num_threads-1)) != 0)
		{
			do ++num_threads;
			while ((num_threads&(num_threads-1)) != 0);
			tp.resize(num_threads);
		}
		size_t n = last - first;
		size_t n_blocks = num_threads;
		size_t block = (n-1)/n_blocks+1;
		for (size_t i = 0; i < n_blocks; ++i)
		{
			size_t start = i*block;
			size_t finish = std::min(n, (i+1)*block);
			tp.enqueue(std::sort<It, Others...>, first + start, first + finish, others...);
		}
		tp.wait();
		for (size_t div = n_blocks/2; div > 1; block *= 2, div /= 2)
		{
			for (size_t i = 0; i < div; ++i)
			{
				size_t start = 2*i*block;
				size_t middle = (2*i+1)*block;
				size_t finish = std::min(n, 2*(i+1)*block);
				tp.enqueue(std::inplace_merge<It, Others...>, first + start, first + middle, first + finish, others...);
			}
			tp.wait();
		}
		size_t start = 0;
		size_t middle = block;
		size_t finish = std::min(n, 2*block);
		std::inplace_merge(first + start, first + middle, first + finish, others...);
	}

} // namespace utils

#endif // UTILS_PARALLEL_SORT_H

























