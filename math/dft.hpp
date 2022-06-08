//  Discrete (fast) Fourier transform
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

#ifndef MATH_DFT_H
#define MATH_DFT_H

#include <array>
#include <algorithm> // copy, reverse
#include <iterator> // random_access_iterator
#include <complex>
#include <numbers> // pi_v

#include "../thread_pool.hpp"

namespace math
{
	constexpr bool access_bit(std::size_t num, unsigned pos)
	{
		return (num >> pos) & 1;
	}

	constexpr std::size_t bit_reverse(std::size_t num, std::size_t n_bits)
	// returns the bit-reverse of num, given the number of bits n_bits
	// if n_bits is 8:
	// 000 -> 000 (0 -> 0)
	// 001 -> 100 (1 -> 4)
	// 010 -> 010 (2 -> 2)
	// 011 -> 110 (3 -> 6)
	// etc...
	{
		using std::size_t;
		size_t res = 0;
		for (size_t i = 0; i < n_bits; ++i)
			res |= access_bit(num, i) << (n_bits-1-i);
		return res;
	}

	constexpr std::size_t num_bits(std::size_t num)
	// least amount of bits needed to represent num
	{
		std::size_t res = 0;
		while (num >>= 1)
			++res;
		return res;
	}

	template <std::forward_iterator InIt, std::random_access_iterator OutIt>
	constexpr void bit_reversal_permutation(InIt first, InIt last, OutIt d_first)
	// performs a bit-reversal permutation
	// requires OutIt to be a random access iterator
	{
		using std::size_t;
		size_t n = std::distance(first, last);
		if ((n&(n-1)) != 0)
			throw("number of elements is not a power of 2!");
		size_t n_bits = num_bits(n);
		for (size_t i = 0; first != last; (void)++first, (void)++i)
			*(d_first + bit_reverse(i, n_bits)) = *first;
	}

	template <typename T>
	struct dft // discrete Fourier transform
	{
		dft() : tmp(1)
		{}

		template <std::random_access_iterator It>
		requires requires (It it)
			{
				{*it} -> std::convertible_to<std::complex<T>>;
			}
		constexpr void fft(It first, It last, std::vector<std::complex<T>>::iterator temp)
		// fast Fourier transform: radix-2 (in-place) algorithm
		// it takes O(n log n) time
		// requires It to be a random access iterator to std::complex<T>
		{
			using std::size_t;
			using std::sin;
			using std::cos;
			size_t n = std::distance(first, last);
			bit_reversal_permutation(first, last, temp);
			std::copy(temp, temp + n, first);
			if ((n&(n-1)) != 0)
				throw("n is not a power of 2!");
			for (size_t i = 1; i <= n/2; i *= 2)
			{
				T ki = -std::numbers::pi_v<T> / i;
				std::complex<T> wi = std::complex<T>(cos(ki), sin(ki));
				for (size_t j = 0; j < n; j += i*2)
				{
					std::complex<T> w = 1;
					for (size_t k = j; k < j+i; ++k)
					{
						std::complex<T> prevk = *(first + k), wvipk = w*(*(first + (i+k)));
						*(first + k) = prevk + wvipk;
						*(first + (i+k)) = prevk - wvipk;
						w *= wi;
					}
				}
			}
		}

		template <typename Vec>
		constexpr Vec fft(Vec x)
		{
			tmp[0].resize(x.size());
			fft(begin(x), end(x), tmp[0].begin());
			return x;
		}

		template <std::random_access_iterator It>
		constexpr void ifft(It first, It last, std::vector<std::complex<T>>::iterator temp)
		// inverse fast Fourier transform
		// requires It to be a random access iterator
		{
			std::reverse(first+1, last);
			fft(first, last, temp);
			T n = last - first;
			for (; first != last; ++first)
				*first /= n;
		}

		template <typename Vec>
		constexpr Vec ifft(Vec x)
		{
			tmp[0].resize(x.size());
			ifft(begin(x), end(x), tmp[0].begin());
			return x;
		}

		template <std::size_t N, typename Vec, bool b_inverse = false>
		requires (N >= 2)
		void fftn(Vec& v, const std::array<std::size_t, N>& n, std::size_t num_threads = 1)
		// N-dimensional fast Fourier transform
		// 1-dim FFTs are performed for each direction, exploiting the separability property of N-dim DFT
		// The 1-dim FFTs are easily parallelized
		{
			using std::size_t;
			tmp.resize(num_threads);
			auto eval_last = [this, num_threads, n, &v](size_t idx, size_t)
				{
					size_t n_ij = 1;
					for (size_t i = 0; i < N-1; ++i)
						n_ij *= n[i];
					size_t block = (n_ij-1)/num_threads + 1;
					for (size_t i = idx*block; (i < n_ij) && (i < (idx+1)*block); ++i)
						if constexpr (!b_inverse)
							fft(begin(v) + i*n[N-1], begin(v) + (i+1)*n[N-1], tmp[idx].begin());
						else
							ifft(begin(v) + i*n[N-1], begin(v) + (i+1)*n[N-1], tmp[idx].begin());
				};
			auto eval = [this, num_threads, n, &v](size_t idx, size_t pos)
				{
					size_t n_ij = 1, npos = 1;
					for (size_t i = 0; i < N; ++i)
						if (i != pos)
							n_ij *= n[i];
					auto it = tmp[idx].begin();
					auto itn = it + n[pos];
					size_t block = (n_ij-1)/num_threads + 1;
					for (size_t j = pos+1; j < N; ++j)
						npos *= n[j];
					std::array<T, N-1> ind{};
					for (size_t i = idx*block; (i < n_ij) && (i < (idx+1)*block); ++i)
					{
						size_t ipart = i;
						for (size_t j = N-1; j > pos; --j)
						{
							ind[j-1] = ipart % n[j];
							ipart /= n[j];
						}
						for (int j = pos-1; j >= 0; --j)
						{
							ind[j] = ipart % n[j];
							ipart /= n[j];
						}
						size_t beg = 0;
						for (size_t j = 0; j < pos; ++j)
						{
							beg *= n[j];
							beg += ind[j];
						}
						beg *= n[pos];
						for (size_t j = pos+1; j < N; ++j)
						{
							beg *= n[j];
							beg += ind[j-1];
						}
						for (size_t k = 0; k < n[pos]; ++k)
							*(it + k) = *(begin(v) + (beg + k*npos));
						if constexpr (!b_inverse)
							fft(it, itn, itn);
						else
							ifft(it, itn, itn);
						for (size_t k = 0; k < n[pos]; ++k)
							*(begin(v) + (beg + k*npos)) = *(it + k);
					}
				};
			if (num_threads == 1)
			{
				tmp[0].resize(n[N-1]);
				eval_last(0, N-1);
				for (int pos = N-2; pos >= 0; --pos)
				{
					tmp[0].resize(n[pos]*2);
					eval(0, pos);
				}
			}
			else
			{
				tp.resize(num_threads);
				for (size_t i = 0; i < num_threads; ++i)
					tmp[i].resize(n[N-1]);
				for (size_t i = 0; i < num_threads; ++i)
					tp.enqueue(eval_last, i, N-1);
				tp.wait();
				for (int pos = N-2; pos >= 0; --pos)
				{
					for (size_t i = 0; i < num_threads; ++i)
						tmp[i].resize(n[pos]*2);
					for (size_t i = 0; i < num_threads; ++i)
						tp.enqueue(eval, i, pos);
					tp.wait();
				}
			}
		}

		template <std::size_t N, typename Vec>
		constexpr void ifftn(Vec& v, const std::array<std::size_t, N>& n, std::size_t num_threads = 1)
		{
			fftn<N, Vec, true>(v, n, num_threads);
		}

		private:

			std::vector<std::vector<std::complex<T>>> tmp;
			// a vector of std::complex<T> for each thread, to minimize chances of false sharing
			thread_pool<std::size_t, std::size_t> tp;
	};

}

#endif // MATH_DFT_H




















