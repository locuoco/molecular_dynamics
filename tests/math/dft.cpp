//  Discrete (fast) Fourier transform tests
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
#include <random> // mt19937, uniform_real_distribution
#include <algorithm> // generate, is_permutation
#include <cassert>
#include <valarray>
#include <array>
#include <numeric> // iota
#include <chrono>

/*

g++ dft.cpp -o dft -std=c++20 -Wall -Wextra -pedantic -Ofast -pthread -fmax-errors=1

*/

#include "../../math/dft.hpp"

void test_bit_reversal()
// test basic properties of bit reversal of an incrementing sequence
{
	std::size_t n = 1 << 3;
	std::vector<int> seq(n), res;
	std::iota(seq.begin(), seq.end(), 0);

	res = seq;

	math::bit_reversal_permutation(res.begin(), res.end());

	for (const auto& elem : res)
		std::cout << elem << ' ';
	std::cout << '\n';

	bool ok = true;
	for (std::size_t i = 0; i < n/2; ++i)
		ok &= (res[i+n/2] == res[i]+1);
	assert(ok);
	assert(std::is_permutation(res.begin(), res.end(), seq.begin()));
}

utils::thread_pool tp;
std::mt19937 mersenne_twister;

template <typename T>
std::uniform_real_distribution<T> dist(-1, 1);

double error_threshold = 1e-15;

template <typename T>
auto gen()
{
	return std::complex<T>(dist<T>(mersenne_twister), dist<T>(mersenne_twister));
}

template <bool b_real = false>
void test_fft_impulse()
// test DFT of unit impulse (Kronecker delta)
{
	using std::sqrt;
	using std::sin;
	using std::cos;
	using std::size_t;
	using std::ptrdiff_t;
	math::dft<double> dft;
	size_t nmin = 1 << 10, nmax = 1ull << 20;
	for (size_t n = nmin; n <= nmax; n <<= 1)
	{
		std::vector<std::complex<double>> x(n);
		double *x_real = reinterpret_cast<double*>(x.data());
		size_t j;

		if constexpr (!b_real)
		{
			std::uniform_int_distribution<size_t> idist(0, n-1);
			j = idist(mersenne_twister);
			for (size_t i = 0; i < n; ++i)
				x[i] = (i == j);
		}
		else
		{
			std::uniform_int_distribution<size_t> idist(0, 2*n-1);
			j = idist(mersenne_twister);
			for (size_t i = 0; i < 2*n; ++i)
				x_real[i] = (i == j);
		}

		x = dft.bfft<false, b_real>(x);

		double rmse = 0;
		if constexpr (!b_real)
		{
			for (size_t i = 0; i < n; ++i)
			{
				double k = -2*double((i*j) % n)*std::numbers::pi/n;
				double diff = norm(x[i]-std::complex<double>(cos(k), sin(k)));
				rmse += diff;
			}
			rmse /= n;
		}
		else
		{
			double diff = x[0].real()-1;
			rmse += diff*diff;
			diff = x[0].imag()-(1-2*(ptrdiff_t(j)%2));
			rmse += diff*diff;
			for (size_t i = 1; i < n; ++i)
			{
				double k = -double((i*j) % (2*n))*std::numbers::pi/n;
				diff = norm(x[i]-std::complex<double>(cos(k), sin(k)));
				rmse += diff;
			}
			rmse /= n+1;
		}
		rmse = sqrt(rmse);
		std::cout << "n = " << n << " -- RMSE = " << rmse << std::endl;

		assert(rmse < error_threshold);
	}
}

template <std::size_t N = 1, bool b_real = false>
void test_fft_ifft()
// test that composing FFT and inverse FFT one obtains identity up to rounding errors
{
	using std::sqrt;
	using std::size_t;
	math::dft<double> dft;
	size_t nN = 1, nmin, nmax;
	if constexpr (N == 1)
		nmin = 1 << 10, nmax = 1ull << 20;
	else
		nmin = 1 << 2, nmax = 1ull << 25;
	dft.compute_twiddle_factors(nmax >> N);
	for (size_t n = nmin; (nN << N) <= nmax; n <<= 1)
	{
		std::array<size_t, N> nlist;
		nN = 1;
		for (size_t i = 0; i < N; ++i)
		{
			nlist[i] = n;
			nN *= n;
		}
		if constexpr (N > 1)
		{
			nlist[N-1] /= 2;
			nlist[N-2] *= 2;
		}
		std::valarray<std::complex<double>> xref(nN), x;

		std::generate(begin(xref), end(xref), gen<double>);

		if constexpr (N == 1)
		{
			x = dft.bfft<false, b_real>(xref);
			x = dft.bfft<true, b_real>(x);
		}
		else
		{
			x = xref;
			dft.bfftn<false, b_real, N>(x, nlist, tp);
			dft.bfftn<true, b_real, N>(x, nlist, tp);
		}

		auto sqr = [](std::complex<double> z)
			{
				return std::complex<double>(norm(z), 0);
			};

		double rmse = sqrt((x-xref).apply(sqr).sum().real()/nN);
		if constexpr (N == 1)
			std::cout << "n = ";
		else
			std::cout << "n^" << N << " = ";
		std::cout << nN << " -- RMSE = " << rmse << std::endl;

		assert(rmse < error_threshold);
	}
}

template <std::size_t N = 1>
void test_fft_perf()
// test FFT performance
{
	using std::size_t;
	math::dft<double> dft;
	size_t n_loops = 10;
	size_t nN = 1, nmin, nmax = (1ull << 22);
	if constexpr (N == 1)
		nmin = 1 << 4;
	else
		nmin = 1 << 2;
	dft.compute_twiddle_factors(nmax >> N);
	for (size_t n = nmin; (nN << N) <= nmax; n <<= 1)
	{
		std::array<size_t, N> nlist;
		nN = 1;
		for (size_t i = 0; i < N; ++i)
		{
			nlist[i] = n;
			nN *= n;
		}

		std::vector<std::complex<double>> x(nN);
		decltype(std::chrono::steady_clock::now()) start, finish;
		double timing = 0;

		for (size_t i = 0; i < n_loops+1; ++i)
		{
			std::generate(begin(x), end(x), gen<double>);

			start = std::chrono::steady_clock::now();
			if constexpr (N == 1)
				x = dft.fft(x);
			else
				dft.fftn<N>(x, nlist, tp);
			finish = std::chrono::steady_clock::now();

			if (i)
				timing += std::chrono::duration<double>(finish-start).count();
		}

		timing /= n_loops;
		if constexpr (N == 1)
			std::cout << "n = "; 
		else
			std::cout << "n^" << N << " = ";
		std::cout << nN << " -- t = " << timing << "s -- x_0 = " << x[0] << std::endl;
	}
}

int main()
{
	std::cout << "=== Bit-reversal test ===\n";
	test_bit_reversal();
	std::cout << "=== FFT impulse test ===\n";
	test_fft_impulse();
	test_fft_impulse<true>();
	std::cout << "=== FFT-IFFT test ===\n";
	test_fft_ifft();
	test_fft_ifft<1, true>();
	test_fft_ifft<2>();
	test_fft_ifft<2, true>();
	test_fft_ifft<3>();
	test_fft_ifft<3, true>();
	test_fft_ifft<4>();
	test_fft_ifft<4, true>();
	std::cout << "=== FFT performance test ===\n";
	test_fft_perf();
	test_fft_perf<2>();
	test_fft_perf<3>();
	test_fft_perf<4>();

	/*using std::size_t;
	math::dft<double> dft;
	size_t n = 8;
	std::valarray<std::complex<double>> x(n*n);

	std::generate(begin(x), end(x), []{ return std::complex(dist<double>(mersenne_twister), 0.); });

	dft.fftn<2>(x, {n, n}, tp);

	for (size_t i = 0; i < n; ++i)
	{
		for (size_t j = 0; j < n; ++j)
			std::cout << x[i*n+j] << ", ";
		std::cout << '\n';
	}*/

	return 0;
}
























