//  Thread pool implementation based on threads
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

#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <thread>
#include <condition_variable>
#include <mutex>
#include <queue>
#include <deque>
#include <tuple>
#include <stop_token>
#include <utility> // move, forward, index_sequence, make_index_sequence
#include <atomic>

template <typename ... Args>
struct task
{
	using func_type = std::function<void(Args...)>;

	func_type f;
	std::tuple<Args...> args;
};

template <typename ... Args>
class notification_queue
{
	std::size_t wait_counter;
	std::deque<task<Args...>> tasks;
	std::mutex mutex;
	std::condition_variable main;
	std::condition_variable_any cond;
	// using std::condition_variable_any for std::stop_token support in C++20

	public:

		notification_queue() : wait_counter{}
		{}

		bool pop(task<Args...>& t, std::stop_token stok)
		// extract a task to perform, and wait if there are none
		// returns false if a stop has been requested
		{
			std::unique_lock lock(mutex);
			++wait_counter;
			main.notify_all();
			cond.wait(lock, stok, [this]{ return !tasks.empty(); });
			// during waiting, mutex is unlocked, then it is locked afterwards/during condition checks
			if (stok.stop_requested())
				return false;
			--wait_counter;
			t = tasks.front();
			tasks.pop_front();
			return true;
		}

		template <typename F, typename ... Brgs>
		void push(F&& f, Brgs&& ... args)
		// put a new task in the queue
		{
			task<Args...> t{std::forward<F>(f), {std::forward<Brgs>(args)...}};
			std::unique_lock lock(mutex);
			tasks.push_back(t);
			lock.unlock();
			cond.notify_one();
		}

		bool busy(std::size_t num_threads)
		// check if there is any task to be completed
		{
			std::unique_lock lock(mutex);
			return tasks.empty() && wait_counter == num_threads;
		}

		void wait(std::size_t num_threads)
		// wait for all tasks to be completed
		// num_threads is the number of threads associated to this queue*
		{
			std::unique_lock lock(mutex);
			main.wait(lock, [this, num_threads] { return tasks.empty() && wait_counter >= num_threads; } );
			if (wait_counter > num_threads)
				throw("wait_counter is wrong!"); // this should never happen, but who knows
		}

		// * For a single-queue thread pool it will be the total number of threads
};

template <typename ... Args>
class thread_pool
// a simple multi-queue thread pool
// using a thread pool rather than direct std::thread calls lead to less overhead
// since threads are reused for many tasks rather than being constructed and destroyed
// the code is very simple and not necessarily optimal
// std::async may implement a thread pool but this is not guaranteed for all STL implementations
{
	std::atomic<std::size_t> index;
	std::deque<notification_queue<Args...>> queues;
	std::deque<std::jthread> threads;
	// jthreads are threads that, on destruction, automatically request a stop and join
	// threads is declared last so that it will be destroyed first, by standard

	template <std::size_t ... Ns>
	void call_impl(task<Args...>& t, std::index_sequence<Ns...>)
	{
		t.f(get<Ns>(t.args)...);
	}

	void call(task<Args...>& t)
	{
		call_impl(t, std::make_index_sequence<sizeof...(Args)>{});
	}

	void thread_loop(std::stop_token stok, std::size_t i)
	{
		while (true)
		{
			task<Args...> t;
			if (!queues[i].pop(t, stok))
				break;

			call(t);
		}
	}

	public:

		thread_pool(std::size_t num_threads = std::thread::hardware_concurrency())
		{
			resize(num_threads);
		}

		std::size_t size() const noexcept
		// return the number of threads (i.e. the size of the thread pool)
		{
			return threads.size();
		}

		void resize(std::size_t num_threads)
		// change the number of threads
		{
			std::size_t old_num_threads = size();
			threads.resize(num_threads);
			queues.resize(num_threads);
			for (std::size_t i = old_num_threads; i < num_threads; ++i)
				threads[i] = std::jthread([this, i](std::stop_token stok) { thread_loop(stok, i); });
		}

		template <typename F, typename ... Brgs>
		void enqueue(F&& f, Brgs&& ... args)
		// put a new task in a queue to be performed by one thread
		{
			auto i = index++;
			queues[i % size()].push(std::forward<F>(f), std::forward<Brgs>(args)...);
		}

		bool busy()
		// check if there is any task to be completed
		{
			bool ret = false;
			for (auto& queue : queues)
				ret &= queue.busy(1);
			return ret;
		}

		void wait()
		// wait for all tasks to be completed
		{
			for (auto& queue : queues)
				queue.wait(1);
		}
};

#endif // THREAD_POOL_H

























