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
#include <concepts> // convertible_to, random_access_iterator, forward_iterator
#include <stdexcept> // invalid_argument

#include "helper.hpp" // powm1
#include "../utils/thread_pool.hpp"

namespace math
{
	template <std::forward_iterator It>
	constexpr void impulse(It first, It last, std::size_t pos)
	// Create a unit impulse signal in the range [first, last),
	// i.e.: 1 at position `pos` and 0 elsewhere.
	// If pos >= std::distance(first, last), then the resulting signal will be 0 everywhere.
	// Note: `It` is required to be a forward (and output) iterator
	{
		for (std::size_t i = 0; first != last; (void)++i, ++first)
			*first = (i == pos);
	}

	constexpr bool access_bit(std::size_t num, unsigned pos)
	// access bit at position `pos` in number `num`
	{
		return (num >> pos) & 1;
	}

	constexpr std::size_t bit_reverse(std::size_t num, std::size_t n_bits)
	// return the bit-reverse of `num`, given the number of bits `n_bits`
	// if `n_bits` is 3:
	// 000 -> 000 (0 -> 0)
	// 001 -> 100 (1 -> 4)
	// 010 -> 010 (2 -> 2)
	// 011 -> 110 (3 -> 6)
	// 100 -> 001 (4 -> 1)
	// 101 -> 101 (5 -> 5)
	// 110 -> 011 (6 -> 3)
	// 111 -> 111 (7 -> 7)
	{
		using std::size_t;
		size_t res = 0;
		for (size_t i = 0; i < n_bits; ++i)
			res |= access_bit(num, i) << (n_bits-1-i);
		return res;
	}

	constexpr std::size_t num_bits(std::size_t num)
	// return the least amount of bits needed to represent `num`
	{
		std::size_t res = 1;
		while (num >>= 1)
			++res;
		return res;
	}

	template <std::random_access_iterator It>
	constexpr void bit_reversal_permutation(It first, It last)
	// perform an in-place bit-reversal permutation of the range [first, last).
	// Throw a `std::invalid_argument` if the number of elements in the range
	// is not a power of 2.
	// example of usage for a container `c`:
	//	bit_reversal_permutation(begin(c), end(c));
	// Note: `It` is required to be a random access iterator
	{
		using std::size_t;
		size_t n = std::distance(first, last);
		if (!math::is_pow_of_2(n))
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
	struct dft
	// discrete Fourier transform (used in PPPM method)
	// contains several methods for the calculation of the discrete Fourier transform
	// thanks to FFT algorithms
	{
		constexpr void compute_twiddle_factors(std::size_t n)
		// compute all twiddle factors for a sequence of power-of-2 length `n` or less.
		// Throw a `std::invalid_argument` if `n` is not a power of 2.
		// twiddle_factor[j][k] = exp(-2 pi i k/(2^j)), where i = imaginary unit.
		// Precomputed twiddle factors are needed to achieve optimal O(sqrt(log(N))) RMS rounding error
		// without sacrificing computation speed (since sin, cos are relatively slow).
		{
			using std::size_t;
			using std::sin;
			using std::cos;
			if (!math::is_pow_of_2(n))
				throw std::invalid_argument(std::string("number of elements is not a power of 2: ") + std::to_string(n));
			size_t n_bits = num_bits(n);
			// check whether twiddle factors are already computed
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
						// the odd elements j of the first quarter of the twiddle factors
						// are computed using e^(ik) = cos(k) + i sin(k) with k = -2 pi j / n
						for (size_t j = 1; j < n_elems/4; j += 2)
						{
							T k = -std::numbers::pi_v<T> * (j/half_n);
							twiddle_factors[i][j] = std::complex<T>(cos(k), sin(k));
						}
						// the even elements j of the first quarter of the twiddle factors
						// are computed from twiddle factors calculated for n/2 (i.e. the previous level)
						for (size_t j = 2; j < n_elems/4; j += 2)
							twiddle_factors[i][j] = twiddle_factors[i-1][j >> 1];
						// all other twiddle factors are calculated using properties of sine and cosine
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

		constexpr std::complex<T> twiddle_factor(std::size_t k, std::size_t n)
		// Return exp(-2 pi i k/n), where i = imaginary unit.
		// Subroutine `compute_twiddle_factors` will throw a `std::invalid_argument` if the
		// number of elements in the range is not a power of 2.
		{
			compute_twiddle_factors(n);
			return twiddle_factors[num_bits(n-1)][k&(n-1)];
		}

		template <std::random_access_iterator It>
		requires requires (It it, std::complex<T> z)
			{
				{z} -> std::convertible_to<decltype(*it)>;
			}
		constexpr void fft(It first, It last)
		// compute a 1-dimensional fast Fourier transform of the range [first, last)
		// using a radix-2 in-place algorithm.
		// `It` is required to be a random access iterator to a type which `std::complex<T>`
		// is convertible to.
		// Subroutine `compute_twiddle_factors` will throw a `std::invalid_argument` if the
		// number of elements in the range is not a power of 2.
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
		// compute the 1-dimensional inverse fast Fourier transform of the range [first, last).
		// The algorithm is in place.
		// `It` is required to be a random access iterator to a type which `std::complex<T>`
		// is convertible to.
		// Subroutine `fft` will throw a `std::invalid_argument` if the number of elements
		// in the range is not a power of 2.
		{
			std::reverse(first+1, last);
			fft(first, last);
			T n = last - first;
			for (; first != last; ++first)
				*first /= n;
		}

		template <std::random_access_iterator It>
		constexpr void rfft(It first, It last)
		// compute a 1-dimensional fast Fourier transform of the subrange [first, last-1)
		// from real input reinterpret-casted to a complex type. A padding of an additional
		// element is required, as the algorithm is *in place* (the output replaces the input).
		// `It` is required to be a random access iterator to a type which `std::complex<T>`
		// is convertible to.
		// While the input is interpreted as real numbers, the output should be interpreted
		// as complex numbers. Since for real input x the DFT has the (even-symmetric) property:
		//	DFT(x)_k = conj(DFT(x)_(n-k))
		// where conj is the complex conjugate and n is the number of real elements the input, then only
		// n/2+1 numbers will need to be stored.
		// Subroutine `fft` will throw a `std::invalid_argument` if the number of elements
		// in the subrange [first, last-1) is not a power of 2.
		{
			fft(first, last-1);
			T xe = first->real();
			T xo = first->imag();
			*first = xe + xo;
			*(last-1) = xe - xo;
			rfft_correct(first, last-1, -1);
		}

		template <std::random_access_iterator It>
		constexpr void irfft(It first, It last)
		// compute a 1-dimensional inverse fast Fourier transform of the range [first, last).
		// The algorithm is *in place* (the output replaces the input).
		// `It` is required to be a random access iterator to a type which `std::complex<T>`
		// is convertible to.
		// This method assumes that the input has the same format as the output of `rfft`.
		// In particular, the first and the last (n+1-th) elements are assumed to be real (imaginary
		// component assumed to be 0).
		// Although the output has a complex type (the same as the input), it should be reinterpreted
		// as real (this can be done easily using `reinterpret_cast`).
		// The last element(s) will be set to 0, as it is a padding used only in Fourier space.
		// Subroutine `rfft_correct` will throw a `std::invalid_argument` if the number of elements
		// in the subrange [first, last-1) is not a power of 2.
		{
			T xe = first->real();
			T xo = (last-1)->real();
			first->real((xe + xo)/2);
			first->imag((xe - xo)/2);
			*(last-1) = T(0);
			rfft_correct(first, last-1, 1);
			ifft(first, last-1);
		}

		template <bool b_inverse, bool b_real, std::random_access_iterator It>
		constexpr void bfft(It first, It last)
		// a wrapper function which performs a FFT of the range [first, last):
		//	if `b_inverse` is false and `b_real` is false, then `fft` is called.
		//	if `b_inverse` is false and `b_real` is true, then `rfft` is called.
		//	if `b_inverse` is true and `b_real` is false, then `ifft` is called.
		//	if `b_inverse` is true and `b_real` is true, then `irfft` is called.
		// The algorithm is *in place* (the output replaces the input).
		// `It` is required to be a random access iterator to a type which `std::complex<T>`
		// is convertible to.
		// A `std::invalid_argument` will be thrown if the number of elements in the
		// range is not supported (power of 2 for fft, ifft, and power of 2 + 1 for rfft
		// and irfft).
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
		// compute the 1-d forward FFT of a container `x`.
		// It is a wrapper of `fft(It, It)`. Functions `begin` and `end` are required to
		// be overloaded for type `Vec` and return a random access iterator to a type
		// which `std::complex<T>` is convertible to.
		// The algorithm is out of place (the input is not replaced).
		// Subroutine `fft` will throw a `std::invalid_argument` if the number of elements
		// is not a power of 2.
		{
			using std::begin;
			using std::end;
			fft(begin(x), end(x));
			return x;
		}

		template <typename Vec>
		constexpr Vec ifft(Vec x)
		// compute the 1-d inverse FFT of a container `x`.
		// It is a wrapper of `ifft(It, It)`. Functions `begin` and `end` are required to
		// be overloaded for type `Vec` and return a random access iterator to a type
		// which `std::complex<T>` is convertible to.
		// The algorithm is out of place (the input is not replaced).
		// Subroutine `ifft` will throw a `std::invalid_argument` if the number of elements
		// is not a power of 2.
		{
			using std::begin;
			using std::end;
			ifft(begin(x), end(x));
			return x;
		}

		template <typename Vec>
		constexpr Vec rfft(Vec x)
		// compute the 1-d forward FFT of a container `x` (with input interpreted as real).
		// It is a wrapper of `rfft(It, It)`. Functions `begin` and `end` are required to
		// be overloaded for type `Vec` and return a random access iterator to a type
		// which `std::complex<T>` is convertible to.
		// The algorithm is out of place (the input is not replaced).
		// Subroutine `rfft` will throw a `std::invalid_argument` if the number of complex elements
		// is not a power of 2 + 1.
		{
			using std::begin;
			using std::end;
			rfft(begin(x), end(x));
			return x;
		}

		template <typename Vec>
		constexpr Vec irfft(Vec x)
		// compute the 1-d inverse FFT of a container `x` (with output interpreted as real).
		// It is a wrapper of `irfft(It, It)`. Functions `begin` and `end` are required to
		// be overloaded for type `Vec` and return a random access iterator to a type
		// which `std::complex<T>` is convertible to.
		// The algorithm is out of place (the input is not replaced).
		// Subroutine `ifft` will throw a `std::invalid_argument` if the number of complex elements
		// is not a power of 2 + 1.
		{
			using std::begin;
			using std::end;
			irfft(begin(x), end(x));
			return x;
		}

		template <bool b_inverse, bool b_real, typename Vec>
		constexpr Vec bfft(Vec x)
		// compute the 1-d FFT of a container `x`, depending on the non-type template arguments.
		// It is a wrapper of `bfft(It, It)`. Functions `begin` and `end` are required to
		// be overloaded for type `Vec` and return a random access iterator to a type
		// which `std::complex<T>` is convertible to.
		// The algorithm is out of place (the input is not replaced).
		// Subroutine `bfft` will throw a `std::invalid_argument` if the number of elements
		// is not supported (power of 2 for fft, ifft, and power of 2 + 1 for rfft and irfft).
		{
			using std::begin;
			using std::end;
			bfft<b_inverse, b_real>(begin(x), end(x));
			return x;
		}

		template <bool b_inverse, bool b_real, std::size_t N, typename Vec>
		requires (N >= 2)
		void bfftn(Vec& v, const std::array<std::size_t, N>& n, utils::thread_pool& tp)
		// N-dimensional fast Fourier transform (performed in-place)
		// 1-dim FFTs are performed for each direction, exploiting the separability property of N-dim DFT.
		// The 1-dim FFTs are easily parallelized using a thread pool `tp`.
		// `v` is a random-access container (like a std::valarray or std::deque) for which the FFT needs to be computed.
		// `n` specifies how the shape of the container should be interpreted.
		// An inverse transform is computed if `b_inverse` is true and the real variants are used
		// if `b_real` is true.
		// A `std::invalid_argument` will be thrown if the shape of the array is not supported (see other methods
		// for more info).
		{
			using std::size_t;
			using std::begin;
			using std::end;
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
						compute_twiddle_factors(2*(n[N-1]-1));
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
						compute_twiddle_factors(2*(n[N-1]-1));
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
		// N-dimensional forward fast Fourier transform (performed in-place).
		// Wrapper of `bfftn` with b_inverse=false and b_real=false (see also `fft`).
		{
			bfftn<false, false, N, Vec>(v, n, tp);
		}

		template <std::size_t N, typename Vec>
		constexpr void ifftn(Vec& v, const std::array<std::size_t, N>& n, utils::thread_pool& tp)
		// N-dimensional inverse fast Fourier transform (performed in-place).
		// Wrapper of `bfftn` with b_inverse=true and b_real=false (see also `ifft`).
		{
			bfftn<true, false, N, Vec>(v, n, tp);
		}

		template <std::size_t N, typename Vec>
		constexpr void rfftn(Vec& v, const std::array<std::size_t, N>& n, utils::thread_pool& tp)
		// N-dimensional forward fast Fourier transform (performed in-place).
		// Wrapper of `bfftn` with b_inverse=false and b_real=true (see also `rfft`).
		{
			bfftn<false, true, N, Vec>(v, n, tp);
		}

		template <std::size_t N, typename Vec>
		constexpr void irfftn(Vec& v, const std::array<std::size_t, N>& n, utils::thread_pool& tp)
		// N-dimensional inverse fast Fourier transform (performed in-place).
		// Wrapper of `bfftn` with b_inverse=true and b_real=true (see also `irfft`).
		{
			bfftn<true, true, N, Vec>(v, n, tp);
		}

		template <std::random_access_iterator It>
		requires requires (It it, std::complex<T> z)
			{
				{z} -> std::convertible_to<decltype(*it)>;
			}
		constexpr void dft_impulse(It first, It last, std::size_t pos)
		// Create the DFT of a unit impulse signal in the range [first, last), i.e.:
		//	 DFT(x)_l = e^(-ik) with k = 2 * (l*pos) * pi / n, and i = imaginary unit
		// `pos` is the index corresponding to the position of impulse before DFT.
		// `It` is required to be a random access iterator to a type which `std::complex<T>`
		// is convertible to.
		// If pos >= std::distance(first, last), a `std::invalid_argument` is thrown.
		// Subroutine `compute_twiddle_factors` will throw a `std::invalid_argument` if the
		// number of elements inside the range is not a power of 2.
		{
			size_t n = std::distance(first, last);
			if (pos >= n)
				throw std::invalid_argument("`pos` argument must be smaller than the size of the range [first, last).");
			compute_twiddle_factors(n);
			size_t n_bits = num_bits(n-1);
			for (std::size_t i = 0; i < n; ++i)
				*(first + i) = twiddle_factors[n_bits][(i*pos)&(n-1)];
		}

		template <std::random_access_iterator It>
		requires requires (It it, std::complex<T> z)
			{
				{z} -> std::convertible_to<decltype(*it)>;
			}
		constexpr void dft_impulse_real(It first, It last, std::size_t pos)
		// Create the DFT of a unit impulse signal in the subrange [first, last-1), i.e.:
		//	 DFT(x)_l = e^(-ik) with k = 2 * (l*pos) * pi / n, and i = imaginary unit
		// `pos` is the index corresponding to the position of impulse before DFT.
		// `It` is required to be a random access iterator to a type which `std::complex<T>`
		// is convertible to.
		// A padding of an additional element is required. As an input it will be ignored, but
		// as an output it will contain the (n/2+1)-th element of the DFT.
		// Differently from `dft_impulse`, the format assumed is the same as the one used
		// for transforms of real input signals.
		// If pos >= 2*std::distance(first, last-1), a `std::invalid_argument` is thrown.
		// Subroutine `compute_twiddle_factors` will throw a `std::invalid_argument` if the
		// number of elements inside the subrange [first, last-1) is not a power of 2.
		{
			size_t n = 2*std::distance(first, last-1);
			if (pos >= n)
				throw std::invalid_argument("`pos` argument must be smaller than 2x the size of the subrange [first, last-1).");
			compute_twiddle_factors(n);
			size_t n_bits = num_bits(n-1);
			// first element must be 1
			// (n/2+1)-th element must be either 1 or -1 depending on j's parity
			*first = T(1);
			*(last-1) = T(math::powm1(pos));
			for (std::size_t i = 1; i < n/2; ++i)
				*(first + i) = twiddle_factors[n_bits][(i*pos)&(n-1)];
		}

		private:

			template <std::random_access_iterator It>
			constexpr void rfft_correct(It first, It last, std::ptrdiff_t sign)
			// used in `rfft` and `irfft`:
			// these methods use the complex-valued `fft` and `ifft` routines but then
			// a correction is required since the input (for forward transform) or the
			// output (for inverse transform) must be interpreted as real.
			// Subroutine `compute_twiddle_factors` will throw a `std::invalid_argument` if the
			// number of elements inside the range is not a power of 2.
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
				for (; first != last; (void)++first, --last)
				{
					twiddle -= sign;
					std::complex<T> zk = *first/T(2), znk = *last/T(2);
					std::complex<T> zkc = conj(zk), znkc = conj(znk);
					std::complex<T> xe = zk + znkc, xo = zk - znkc;
					std::complex<T> xen = znk + zkc, xon = znk - zkc;

					*first = xe + xo * (*twiddle);
					*last = xen + xon * conj(*twiddle);
				}
				first->imag(-first->imag());
			}

			std::vector<std::vector<std::complex<T>>> tmp;
			// a vector of std::complex<T> for each thread, to minimize chances of false sharing

			std::vector<std::vector<std::complex<T>>> twiddle_factors;
	};

}

#endif // MATH_DFT_H




















