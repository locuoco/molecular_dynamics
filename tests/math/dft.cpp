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
#include <random> // mt19937_64
#include <algorithm> // generate, is_permutation
#include <cassert>
#include <valarray>
#include <array>
#include <numeric> // iota
#include <chrono>
#include <cmath> // sin, cos, sqrt

/*

Compilation:
g++ dft.cpp -o dft -std=c++20 -Wall -Wextra -pedantic -Ofast -pthread -fmax-errors=1

*/

#include "../../math/dft.hpp"

utils::thread_pool tp;
constexpr double error_threshold = 1e-15;

std::mt19937_64 mersenne_twister(1234); // seed set to 1234
// mersenne_twister will behave in the same way for all compilers/runs

auto gen()
// generate a pseudo-random complex number x + iy with x,y in [-1, 1]
{
	double two63 = 1ull << 63; // 2 ** 63
	return std::complex(mersenne_twister()/two63-1, mersenne_twister()/two63-1);
}

auto sqr(std::complex<double> z)
// calculate the norm squared of complex number `z`.
// Return a complex number with 0 imaginary part.
{
	return std::complex<double>(norm(z), 0);
}

void test_bit_reversal_incr(std::size_t n)
// test basic property of bit reversal of an incrementing sequence:
//	B(x_(i+n/2)) == B(x_i + 1)  for 0 <= i < n/2  (if n is a power of 2)
// where x is an incrementing sequence and B is the bit reversal permutation.
// The test passes if this property is satisfied.
// `n` is the length of the sequence to test, which must be a power of 2.
{
	std::vector<int> seq(n);
	// iota creates an incrementing (starting from 0)
	std::iota(begin(seq), end(seq), 0);

	math::bit_reversal_permutation(begin(seq), end(seq));

	bool ok = true;
	for (std::size_t i = 0; i < n/2; ++i)
		ok &= (seq[i+n/2] == seq[i]+1);
	assert(ok);
}

void test_bit_reversal_permutation(std::size_t n)
// Starting from an incrementing sequence, the test passes if
// its bit reversal permutation is indeed a permutation.
// `n` is the length of the sequence to test, which must be a power of 2.
{
	std::vector<int> seq(n), res;
	std::iota(begin(seq), end(seq), 0);

	res = seq;

	math::bit_reversal_permutation(begin(res), end(res));

	assert(std::is_permutation(begin(res), end(res), begin(seq)));
}

void test_fft_impulse(std::size_t n, std::size_t j)
// test DFT of unit impulse (Kronecker delta)
// `n` is the length of the sequence to test, which must be a power of 2.
// `j` is the element that is 1 (all others are 0).
// Note that j < n must be true for the test to be meaningful.
// The test passes if the Fourier transform is given by (up to rounding errors):
// DFT(x)_l = e^(ik) with k = 2 * (l*j) * pi / n
{
	using std::sqrt;
	using std::sin;
	using std::cos;
	using std::size_t;

	math::dft<double> dft;
	std::vector<std::complex<double>> x(n);

	for (size_t i = 0; i < n; ++i)
		x[i] = (i == j);

	dft.fft(begin(x), end(x));

	double rmse = 0;
	for (size_t i = 0; i < n; ++i)
	{
		double k = -2*double((i*j) % n)*std::numbers::pi/n;
		double diff = norm(x[i]-std::complex<double>(cos(k), sin(k)));
		rmse += diff;
	}
	rmse = sqrt(rmse/n);

	assert(rmse < error_threshold);
}

void test_fft_impulse_real(std::size_t n, std::size_t j)
// test DFT of unit impulse (Kronecker delta) for real input
// `n` is the length of the sequence to test, which must be a power of 2 greater than 1.
// `j` is the element that is 1 (all others are 0).
// Note that j < n must be true for the test to be meaningful.
// Note also that n real numbers can be reinterpreted as n/2 complex ones.
// The test passes if the Fourier transform is given by (up to rounding errors):
// DFT(x)_l = e^(ik) with k = 2 * (l*j) * pi / n, and i = imaginary unit
{
	using std::sqrt;
	using std::sin;
	using std::cos;
	using std::size_t;

	math::dft<double> dft;
	std::vector<std::complex<double>> x(n/2);
	double *x_real = reinterpret_cast<double*>(x.data());

	for (size_t i = 0; i < n; ++i)
		x_real[i] = (i == j);

	dft.rfft(begin(x), end(x));

	double rmse = 0;
	// first element must be 1
	double diff = x[0].real()-1;
	rmse += diff*diff;
	// (n/2+1)-th element must be either 1 or -1 depending on j's parity
	diff = x[0].imag()-(1-2*(std::ptrdiff_t(j)%2));
	rmse += diff*diff;
	for (size_t i = 1; i < n/2; ++i)
	{
		double k = -2*double((i*j) % n)*std::numbers::pi/n;
		diff = norm(x[i]-std::complex<double>(cos(k), sin(k)));
		rmse += diff;
	}
	rmse = sqrt(rmse/(n/2+1));

	assert(rmse < error_threshold);
}

void test_fft_ifft(std::size_t n)
// test that composing FFT and inverse FFT one obtains identity.
// `n` is the length of the sequence to test, which must be a power of 2.
// Starting from a vector x, the test passes if IFFT(FFT(x)) == x
// up to rounding errors.
{
	using std::sqrt;

	math::dft<double> dft;
	std::valarray<std::complex<double>> xref(n), x;

	std::generate(begin(xref), end(xref), gen);

	x = dft.fft(xref);
	x = dft.ifft(x);

	double rmse = sqrt((x-xref).apply(sqr).sum().real()/n);

	assert(rmse < error_threshold);
}

void test_rfft_irfft(std::size_t n)
// test that composing FFT and inverse FFT one obtains identity.
// `n` is the length of the sequence to test, which must be a power of 2.
// Starting from a vector x, the test passes if IFFT(FFT(x)) == x
// up to rounding errors.
{
	using std::sqrt;

	math::dft<double> dft;
	std::valarray<std::complex<double>> xref(n), x;

	std::generate(begin(xref), end(xref), gen);

	x = dft.rfft(xref);
	x = dft.irfft(x);

	double rmse = sqrt((x-xref).apply(sqr).sum().real()/n);

	assert(rmse < error_threshold);
}

void test_fft_ifft3(const std::array<std::size_t, 3>& n_shape)
// test that composing FFT and inverse FFT one obtains identity.
// `n_shape` is the shape of the array to test, whose entries must be powers of 2.
// Starting from an array x, the test passes if IFFT(FFT(x)) == x
// up to rounding errors.
{
	using std::sqrt;
	using std::size_t;

	size_t n = n_shape[0]*n_shape[1]*n_shape[2];

	math::dft<double> dft;
	std::valarray<std::complex<double>> xref(n), x;

	std::generate(begin(xref), end(xref), gen);

	x = xref;
	dft.fftn<3>(x, n_shape, tp);
	dft.ifftn<3>(x, n_shape, tp);

	double rmse = sqrt((x-xref).apply(sqr).sum().real()/n);

	assert(rmse < error_threshold);
}

void test_rfft_irfft3(const std::array<std::size_t, 3>& n_shape)
// test that composing FFT and inverse FFT one obtains identity.
// `n_shape` is the shape of the array to test, whose entries must be powers of 2.
// Starting from an array x, the test passes if IFFT(FFT(x)) == x
// up to rounding errors.
{
	using std::sqrt;
	using std::size_t;

	size_t n = n_shape[0]*n_shape[1]*n_shape[2];

	math::dft<double> dft;
	std::valarray<std::complex<double>> xref(n), x;

	std::generate(begin(xref), end(xref), gen);

	x = xref;
	dft.rfftn<3>(x, n_shape, tp);
	dft.irfftn<3>(x, n_shape, tp);

	double rmse = sqrt((x-xref).apply(sqr).sum().real()/n);

	assert(rmse < error_threshold);
}

int main()
{
	test_bit_reversal_incr(16);
	test_bit_reversal_permutation(8);
	test_fft_impulse(1<<20, 2468);
	test_fft_impulse(1<<20, 3579);
	test_fft_impulse_real(1<<20, 4680);
	test_fft_impulse_real(1<<20, 1357);
	test_fft_ifft(1<<20);
	test_rfft_irfft(1<<20);
	test_fft_ifft3({16, 64, 32});
	test_rfft_irfft3({16, 64, 32});

	std::cout << "All tests passed successfully!" << std::endl;

	return 0;
}
























