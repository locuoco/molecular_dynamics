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
#include <algorithm> // copy, reverse, iter_swap
#include <iterator> // random_access_iterator
#include <complex>
#include <numbers> // numbers::pi_v
#include <concepts> // convertible_to
#include <stdexcept> // invalid_argument

#include "../utils/thread_pool.hpp"

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
		std::size_t res = 1;
		while (num >>= 1)
			++res;
		return res;
	}

	template <std::random_access_iterator It>
	constexpr void bit_reversal_permutation(It first, It last)
	// performs a bit-reversal permutation in-place
	// requires It to be a random access iterator
	{
		using std::size_t;
		size_t n = std::distance(first, last);
		if ((n&(n-1)) != 0 || !n)
			throw std::invalid_argument("number of elements is not a power of 2!");
		size_t n_bits = num_bits(n-1);
		for (size_t i = 0; i < n; ++i)
		{
			size_t j = bit_reverse(i, n_bits);
			if (i < j) // bit reversal is symmetric!
				std::iter_swap(first + i, first + j);
		}
	}

	template <typename T>
	struct dft // discrete Fourier transform
	{
		constexpr void compute_twiddle_factors(std::size_t n)
		// twiddle_factor[j][k] = exp(-2 pi i k/(2^j)), where i = imaginary unit
		// precomputed twiddle factors are needed to achieve optimal O(sqrt(log(N))) RMS rounding error
		// without sacrificing computation speed (since sin, cos are relatively slow)
		{
			using std::size_t;
			using std::sin;
			using std::cos;
			if ((n&(n-1)) != 0 || !n)
				throw std::invalid_argument(std::string("number of elements is not a power of 2: ") + std::to_string(n));
			size_t n_bits = num_bits(n);
			if (n_bits > twiddle_factors.size())
			{
				twiddle_factors.resize(n_bits);
				size_t n_elems = 1;
				for (size_t i = 0; i < n_bits; ++i)
				{
					if (twiddle_factors[i].size() != n_elems)
					{
						twiddle_factors[i].resize(n_elems);
						T half_n = n_elems/2;
						twiddle_factors[i][0] = 1;
						if (n_elems > 1)
							twiddle_factors[i][n_elems/2] = -1;
						if (n_elems >= 4)
							twiddle_factors[i][n_elems/4] = std::complex<T>(0, -1);
						for (size_t j = 1; j < n_elems/4; j += 2)
						{
							T k = -std::numbers::pi_v<T> * (j/half_n);
							twiddle_factors[i][j] = std::complex<T>(cos(k), sin(k));
						}
						for (size_t j = 2; j < n_elems/4; j += 2)
							twiddle_factors[i][j] = twiddle_factors[i-1][j >> 1];
						for (size_t j = n_elems/4+1; j < n_elems/2; ++j)
						{
							std::complex<T> t(twiddle_factors[i][n_elems/2-j]);
							twiddle_factors[i][j] = std::complex<T>(-t.real(), t.imag());
						}
						for (size_t j = n_elems/2+1; j < n_elems; ++j)
							twiddle_factors[i][j] = conj(twiddle_factors[i][n_elems-j]);
					}
					n_elems *= 2;
				}
			}
		}

		template <std::random_access_iterator It>
		requires requires (It it, std::complex<T> z)
			{
				{z} -> std::convertible_to<decltype(*it)>;
			}
		constexpr void fft(It first, It last)
		// fast Fourier transform: radix-2 (in-place) algorithm
		// it takes O(n log n) time
		// requires It to be a random access iterator to std::complex<T>
		{
			using std::size_t;
			using std::ptrdiff_t;
			size_t n = std::distance(first, last);
			compute_twiddle_factors(n);
			bit_reversal_permutation(first, last);
			size_t m = 1;
			for (size_t i = 1; i <= n/2; ++m, i *= 2)
				for (size_t j = 0; j < n; j += i*2)
				{
					auto twiddle = twiddle_factors[m].begin();
					for (size_t k = j; k < j+i; ++k)
					{
						std::complex<T> x0 = *(first + k), x1 = (*twiddle++)*(*(first + (i+k)));
						*(first + k) = x0 + x1;
						*(first + (i+k)) = x0 - x1;
					}
				}
		}

		template <std::random_access_iterator It>
		constexpr void ifft(It first, It last)
		// inverse fast Fourier transform
		// requires It to be a random access iterator to std::complex<T>
		{
			std::reverse(first+1, last);
			fft(first, last);
			T n = last - first;
			for (; first != last; ++first)
				*first /= n;
		}

		template <std::random_access_iterator It>
		constexpr void rfft_dit(It first, It last, std::ptrdiff_t sign)
		{
			using std::size_t;
			size_t half_n = std::distance(first, last);
			++first;
			--last;
			size_t n_bits = num_bits(half_n);
			decltype(twiddle_factors[0].begin()) twiddle;
			compute_twiddle_factors(half_n*2);
			if (sign > 0)
				twiddle = twiddle_factors[n_bits].end() - half_n/2;
			else
				twiddle = twiddle_factors[n_bits].begin() + half_n/2;
			for (; first != last; ++first, --last)
			{
				twiddle -= sign;
				std::complex<T> zk = *first/T(2), znk = *last/T(2);
				std::complex<T> zkc = conj(zk), znkc = conj(znk);
				std::complex<T> xe = zk + znkc, xo = zk - znkc;
				std::complex<T> xen = znk + zkc, xon = znk - zkc;

				*first = xe + xo * (*twiddle);
				*last = xen + xon * conj(*twiddle);
			}
			first -> imag(-first -> imag());
		}

		template <std::random_access_iterator It>
		constexpr void rfft(It first, It last)
		// real fast Fourier transform
		// requires It to be a random access iterator to to std::complex<T>
		// the even elements must be put in real parts and odd elements in imaginary parts
		{
			fft(first, last);
			T xe = first -> real();
			T xo = first -> imag();
			first -> real(xe + xo);
			first -> imag(xe - xo);
			// ^ the first and the last elements are on the first complex number,
			// since they are both real
			rfft_dit(first, last, -1);
		}

		template <std::random_access_iterator It>
		constexpr void irfft(It first, It last)
		// inverse real fast Fourier transform
		// requires It to be a random access iterator to to std::complex<T>
		// the first and the last elements are on the first complex number
		{
			T xe = first -> real();
			T xo = first -> imag();
			first -> real((xe + xo)/2);
			first -> imag((xe - xo)/2);
			rfft_dit(first, last, 1);
			ifft(first, last);
		}

		template <bool b_inverse, bool b_real, std::random_access_iterator It>
		constexpr void bfft(It first, It last)
		{
			if constexpr (!b_real)
			{
				if constexpr (!b_inverse)
					fft(first, last);
				else
					ifft(first, last);
			}
			else
			{
				if constexpr (!b_inverse)
					rfft(first, last);
				else
					irfft(first, last);
			}
		}

		template <typename Vec>
		constexpr Vec fft(Vec x)
		{
			fft(begin(x), end(x));
			return x;
		}

		template <typename Vec>
		constexpr Vec ifft(Vec x)
		{
			ifft(begin(x), end(x));
			return x;
		}

		template <typename Vec>
		constexpr Vec rfft(Vec x)
		{
			rfft(begin(x), end(x));
			return x;
		}

		template <typename Vec>
		constexpr Vec irfft(Vec x)
		{
			irfft(begin(x), end(x));
			return x;
		}

		template <bool b_inverse, bool b_real, typename Vec>
		constexpr Vec bfft(Vec x)
		{
			bfft<b_inverse, b_real>(begin(x), end(x));
			return x;
		}

		template <bool b_inverse, bool b_real, std::size_t N, typename Vec>
		requires (N >= 2)
		void bfftn(Vec& v, const std::array<std::size_t, N>& n, utils::thread_pool& tp)
		// N-dimensional fast Fourier transform (performed in-place)
		// 1-dim FFTs are performed for each direction, exploiting the separability property of N-dim DFT
		// The 1-dim FFTs are easily parallelized
		{
			using std::size_t;
			auto num_threads = tp.size();
			tmp.resize(num_threads);
			auto eval_last = [this, num_threads, n, &v](size_t idx)
				{
					size_t n_ij = 1;
					for (size_t i = 0; i < N-1; ++i)
						n_ij *= n[i];
					size_t block = (n_ij-1)/num_threads + 1;
					for (size_t i = idx*block; (i < n_ij) && (i < (idx+1)*block); ++i)
						bfft<b_inverse, b_real>(begin(v) + i*n[N-1], begin(v) + (i+1)*n[N-1]);
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
						bfft<b_inverse, false>(it, itn);
						for (size_t k = 0; k < n[pos]; ++k)
							*(begin(v) + (beg + k*npos)) = *(it + k);
					}
				};
			if (num_threads == 1)
			{
				if constexpr (!b_inverse)
					eval_last(0);
				for (int pos = N-2; pos >= 0; --pos)
				{
					tmp[0].resize(n[pos]);
					eval(0, pos);
				}
				if constexpr (b_inverse)
					eval_last(0);
			}
			else
			{
				if constexpr (!b_inverse)
				{
					if constexpr (b_real)
						compute_twiddle_factors(2*n[N-1]);
					else
						compute_twiddle_factors(n[N-1]);
					for (size_t i = 0; i < num_threads; ++i)
						tp.enqueue(eval_last, i);
					tp.wait();
				}
				for (int pos = N-2; pos >= 0; --pos)
				{
					compute_twiddle_factors(n[pos]);
					for (size_t i = 0; i < num_threads; ++i)
						tmp[i].resize(n[pos]);
					for (size_t i = 0; i < num_threads; ++i)
						tp.enqueue(eval, i, pos);
					tp.wait();
				}
				if constexpr (b_inverse)
				{
					if constexpr (b_real)
						compute_twiddle_factors(2*n[N-1]);
					else
						compute_twiddle_factors(n[N-1]);
					for (size_t i = 0; i < num_threads; ++i)
						tp.enqueue(eval_last, i);
					tp.wait();
				}
			}
		}

		template <std::size_t N, typename Vec>
		constexpr void fftn(Vec& v, const std::array<std::size_t, N>& n, utils::thread_pool& tp)
		{
			bfftn<false, false, N, Vec>(v, n, tp);
		}

		template <std::size_t N, typename Vec>
		constexpr void ifftn(Vec& v, const std::array<std::size_t, N>& n, utils::thread_pool& tp)
		{
			bfftn<true, false, N, Vec>(v, n, tp);
		}

		template <std::size_t N, typename Vec>
		constexpr void rfftn(Vec& v, const std::array<std::size_t, N>& n, utils::thread_pool& tp)
		// adjacent real values must be put inside a single complex number along the last dimension
		{
			bfftn<false, true, N, Vec>(v, n, tp);
		}

		template <std::size_t N, typename Vec>
		constexpr void irfftn(Vec& v, const std::array<std::size_t, N>& n, utils::thread_pool& tp)
		{
			bfftn<true, true, N, Vec>(v, n, tp);
		}

		private:

			std::vector<std::vector<std::complex<T>>> tmp;
			// a vector of std::complex<T> for each thread, to minimize chances of false sharing
			std::vector<std::vector<std::complex<T>>> twiddle_factors;
	};

}

#endif // MATH_DFT_H




















