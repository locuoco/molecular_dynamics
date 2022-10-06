//  Thread pool implementation based on std::jthread
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

#ifndef UTILS_THREAD_POOL_H
#define UTILS_THREAD_POOL_H

#include <thread>
#include <condition_variable>
#include <mutex> // mutex, unique_lock, lock_guard
#include <deque>
#include <tuple>
#include <stop_token>
#include <utility> // move, forward, index_sequence, make_index_sequence, declval
#include <atomic>
#include <memory> // unique_ptr
#include <functional> // function

namespace utils
{
	struct task_base
	{
		virtual void perform() = 0;
		// perform a task

		virtual operator bool() = 0;
		// check if the task is well-formed

		virtual ~task_base() = default;
		// virtual trivial destructor
	};

	template <typename Func>
	class task_wrapper : public task_base
	// needed for type erasure, so that the `task` class is not a template class.
	// The advantage of this is that a std::vector of tasks may be created independently
	// of the number of parameters that need to be stored in each object.
	// To achieve this, dynamical polymorphism is required
	{
			template <typename T>
			struct args_tuple {};

			template <typename ... Args>
			struct args_tuple<std::function<void(Args...)>>
			{
				using type = std::tuple<Args...>;
			};

			template <std::size_t ... Ns>
			void call(std::index_sequence<Ns...>)
			// call function `f` with arguments `args`
			{
				f(get<Ns>(args)...);
			}

			using func_type = decltype(std::function{std::declval<Func>()});
			using args_type = args_tuple<func_type>::type;

			func_type f;
			args_type args;

		public:

			template <typename F, typename ... Args>
			task_wrapper(F&& f, Args&& ... args)
			// constructor:
			//	`f` is the function object associated to the task.
			//	`args...` are the arguments of the function object.
				: f(std::forward<F>(f)), args{std::forward<Args>(args)...}
			{}

			void perform() override
			// perform the tasks f(args...)
			{
				call(std::make_index_sequence<std::tuple_size_v<args_type>>{});
			}

			operator bool() override
			// check if `f` is callable
			{
				return bool(f);
			}
	};

	class task
	{
		std::unique_ptr<task_base> pt;

		public:

			task() : pt(nullptr) {}
			// default constructor

			template <typename F, typename ... Args>
			task(F&& f, Args&& ... args)
			// add a task by storing a function object `f` and its arguments `args`
				: pt(std::make_unique<task_wrapper<F>>(std::forward<F>(f), std::forward<Args>(args)...))
			{}

			void perform()
			// perform the task
			// undefined behaviour if `pt` is nullptr
			// throws std::bad_function_call if `pt` is a std::function object with no target callable
			{
				pt->perform();
			}

			operator bool()
			// check if pt is callable
			{
				return pt && bool(*pt);
			}
	};

	class notification_queue
	{
		std::size_t work_counter{};
		std::deque<task> tasks;
		std::mutex mutex;
		std::condition_variable main;
		std::condition_variable_any cond;
		// using std::condition_variable_any for std::stop_token support in C++20

		public:

			bool pop(task& t, std::stop_token stok)
			// extract a task `t` to perform, and wait if there is none
			// return false if a stop has been requested.
			// `stok` is the stop token associated to a std::jthread.
			{
				std::unique_lock lock(mutex);
				cond.wait(lock, stok, [this]{ return !tasks.empty(); });
				// during waiting, mutex is unlocked, then it is locked afterwards/during condition checks
				if (stok.stop_requested())
					return false;
				++work_counter;
				t = std::move(tasks.front());
				tasks.pop_front();
				return true;
			}

			bool try_pop(task& t)
			// try to extract a task `t` to perform without waiting
			// return false if failing to extract a task (the mutex is already locked in another thread)
			// or if there is no task left
			{
				std::unique_lock lock(mutex, std::try_to_lock);
				if (!lock || tasks.empty())
					return false;
				++work_counter;
				t = std::move(tasks.front());
				tasks.pop_front();
				return true;
			}

			void finished_task()
			// to be called when the task extracted with pop has terminated
			{
				{
					std::lock_guard lock(mutex);
					--work_counter;
				}
				main.notify_all();
			}

			template <typename F, typename ... Args>
			void push(F&& f, Args&& ... args)
			// put a new task in the queue
			{
				{
					task t(std::forward<F>(f), std::forward<Args>(args)...);
					std::lock_guard lock(mutex);
					tasks.push_back(std::move(t));
				}
				cond.notify_one();
			}

			template <typename F, typename ... Args>
			bool try_push(F&& f, Args&& ... args)
			// try to put a new task in the queue without blocking the current thread
			// return false if failing
			{
				std::unique_lock lock(mutex, std::try_to_lock);
				if (!lock)
					return false;
				task t(std::forward<F>(f), std::forward<Args>(args)...);
				tasks.push_back(std::move(t));
				lock.unlock();
				cond.notify_one();
				return true;
			}

			bool empty()
			// check if the task queue is empty
			{
				std::lock_guard lock(mutex);
				return tasks.empty();
			}

			bool busy()
			// check if there is any task to be completed
			{
				std::lock_guard lock(mutex);
				return !tasks.empty() || work_counter > 0;
			}

			void wait()
			// wait for all tasks to be completed
			{
				std::unique_lock lock(mutex);
				main.wait(lock, [this] { return tasks.empty() && work_counter == 0; } );
			}
	};

	class thread_pool
	// A simple multi-queue thread pool, based on "Better Code: Concurrency", Sean Parent, 2017
	// using a thread pool rather than direct std::thread calls lead to less overhead
	// since threads are reused for many tasks rather than being constructed and destroyed.
	// std::async may implement a thread pool but this is not guaranteed for all compilers.
	{
		struct worker
		{
			notification_queue queue;
			std::jthread thread;
			// jthreads are threads that, on destruction, automatically request a stop and join;
			// 'thread' is declared last so that it will be destroyed first, by standard
		};

		std::atomic<std::size_t> index;
		std::deque<worker> workers;
		// using std::deque, which does not reallocate when adding new elements

		void thread_loop(std::stop_token stok, std::size_t i)
		// `stok` is the stop token associated to the current thread
		// `i` is the index of the current thread
		{
			while (true)
			{
				task t;
				std::size_t chosen_ind = i;
				// simple task stealing algorithm
				for (unsigned j = 0; j < size(); ++j)
					if (workers[(i+j) % size()].queue.try_pop(t))
					{
						chosen_ind = (i+j) % size();
						break;
					}
				if (!t && !workers[i].queue.pop(t, stok))
					break;

				t.perform();
				workers[chosen_ind].queue.finished_task();
			}
		}

		public:

			thread_pool(std::size_t num_threads = std::thread::hardware_concurrency())
			// constructor:
			// `num_threads` is the number of threads to launch and use
			{
				resize(num_threads);
			}

			std::size_t size() const noexcept
			// return the number of threads (i.e. the size of the thread pool)
			{
				return workers.size();
			}

			void resize(std::size_t num_threads)
			// change the number of threads
			{
				wait();
				std::size_t old_num_threads = size();
				workers.resize(num_threads);
				for (std::size_t i = old_num_threads; i < num_threads; ++i)
					workers[i].thread = std::jthread([this, i](std::stop_token stok) { thread_loop(stok, i); });
			}

			template <typename F, typename ... Args>
			void enqueue(F&& f, Args&& ... args)
			// put a new task in a queue to be performed by one thread
			// `f` is a function object (either a function pointer or a lambda)
			// `args...` are the arguments (if any)
			{
				auto i = index++; // index is atomic so that enqueue can be called by any thread (it is thread-safe)
				// simple scheduling mechanism
				for (unsigned j = 0; j < size(); ++j)
					if (workers[(i+j) % size()].queue.try_push(std::forward<F>(f), std::forward<Args>(args)...))
						return;
				workers[i % size()].queue.push(std::forward<F>(f), std::forward<Args>(args)...);
			}

			bool busy()
			// check if there is any task to be completed
			{
				for (auto& w : workers)
					if (w.queue.busy())
						return true;
				return false;
			}

			void wait()
			// wait for all tasks to be completed
			// needed for synchronization
			{
				for (auto& w : workers)
					w.queue.wait();
			}
	};

} // namespace utils

#endif // UTILS_THREAD_POOL_H

























