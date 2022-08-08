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
#include <numbers> // numbers::pi

/*

Compilation:
g++ dft.cpp -o dft -std=c++20 -Wall -Wextra -pedantic -Ofast -pthread -fmax-errors=1

*/

#include "../../math/dft.hpp"
#include "../../math/helper.hpp" // rms

utils::thread_pool tp;
constexpr double error_threshold = 1e-15;

std::mt19937_64 mersenne_twister;
// mersenne_twister will behave in the same way for all compilers/runs

auto gen()
// generate a pseudo-random complex number x + iy with x,y in [-1, 1]
{
	double two63 = 1ull << 63; // 2 ** 63
	return std::complex(mersenne_twister()/two63-1, mersenne_twister()/two63-1);
}

void test_gen()
// test that the values generated with `gen` have real and imaginary
// parts between -1 and 1
{
	using std::abs;
	mersenne_twister.seed(1234);
	std::valarray<std::complex<double>> x(1000);

	std::ranges::generate(x, gen);

	for (const auto& elem : x)
		assert(abs(elem.real()) <= 1 && abs(elem.imag()) <= 1);
}

void test_bit_reversal_incr(std::size_t n)
// test basic property of bit reversal of an incrementing sequence:
//	B(x_(i+n/2)) == B(x_i + 1)  for 0 <= i < n/2  (if n is a power of 2)
// where x is an incrementing sequence and B is the bit reversal permutation.
// Another property, that the first n/2 elements are even and the last n/2
// are odd, is also used.
// The test passes if these property is satisfied.
// `n` is the length of the sequence to test, which must be a power of 2
// greater than 1.
{
	std::valarray<int> x(n);
	// iota creates an incrementing (starting from 0)
	std::iota(begin(x), end(x), 0);

	math::bit_reversal_permutation(begin(x), end(x));

	// test that the number of even elements are exactly n/2
	assert(as_const(x)[x%2 == 0].size() == n/2);
	// `as_const` forces the compiler to choose the const overload of operator[]
	// which return a valarray instead of a mask_array (which has no size()
	// method).

	assert((x[x%2 == 1] == as_const(x)[x%2 == 0]+1).min());
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

	assert(std::ranges::is_permutation(res, seq));
}

void test_twiddle_factors(std::size_t k, std::size_t n)
// Test that dft.twiddle_factor(k, n) == exp(-2 pi i k/n), where i = imaginary unit.
// `n` must be a power of 2.
{
	math::dft<double> dft;
	std::complex exponent(0., -2*std::numbers::pi*k/n);

	assert(abs(dft.twiddle_factor(k, n) - exp(exponent)) < error_threshold);
}

void test_fft_impulse(std::size_t n, std::size_t j)
// test DFT of unit impulse (Kronecker delta)
// `n` is the length of the sequence to test, which must be a power of 2.
// `j` is the element that is 1 (all others are 0).
// Note that j < n must be true for the test to be meaningful.
// The test passes if the Fourier transform is given by (up to rounding errors):
// 	DFT(x)_l = e^(-ik) with k = 2 * (l*j) * pi / n
{
	math::dft<double> dft;
	std::valarray<std::complex<double>> x(n), xref(n);

	math::impulse(begin(x), end(x), j);
	dft.fft(begin(x), end(x));

	dft.dft_impulse(begin(xref), end(xref), j);

	assert(math::rms(x-xref) < error_threshold);
}

void test_fft_impulse_real(std::size_t n, std::size_t j)
// test DFT of unit impulse (Kronecker delta) for real input
// `n` is the length of the sequence to test, which must be a power of 2 greater than 1.
// `j` is the element that is 1 (all others are 0).
// Note that j < n must be true for the test to be meaningful.
// Note also that n real numbers can be reinterpreted as n/2 complex ones.
// The test passes if the Fourier transform is given by (up to rounding errors):
// 	DFT(x)_l = e^(-ik) with k = 2 * (l*j) * pi / n, and i = imaginary unit
{
	math::dft<double> dft;
	std::valarray<std::complex<double>> x(n/2+1), xref(n/2+1);
	double *x_real = reinterpret_cast<double*>(&x[0]);

	math::impulse(x_real, x_real+n, j);

	dft.rfft(begin(x), end(x));
	dft.dft_impulse_real(begin(xref), end(xref), j);

	assert(math::rms(x-xref) < error_threshold);
}

void test_fft_linearity(std::size_t n)
// test linearity property of FFT up to rounding errors:
//	FFT(a x + b y) == a FFT(x) + b FFT(y)
// `n` is the length of the sequence to test, which must be a power of 2.
{
	mersenne_twister.seed(1234);
	math::dft<double> dft;
	std::valarray<std::complex<double>> x(n), y(n);
	std::complex<double> a = gen(), b = gen();

	std::ranges::generate(x, gen);
	std::ranges::generate(y, gen);

	std::valarray<std::complex<double>> z = a*x + b*y;

	x = dft.fft(x);
	y = dft.fft(y);
	z = dft.fft(z);

	assert(math::rms(z - (a*x + b*y)) < 1e-12);
}

void test_rfft_linearity(std::size_t n)
// test linearity property of FFT up to rounding errors:
//	FFT(a x + b y) == a FFT(x) + b FFT(y)
// `n` is the length of the sequence to test, which must be a power of 2.
// Real-signal input variant.
{
	mersenne_twister.seed(1234);
	math::dft<double> dft;
	std::valarray<std::complex<double>> x(n+1), y(n+1);
	std::complex<double> c = gen();
	std::complex<double> a = c.real(), b = c.imag();
	// a and b must are set to real to avoid mixing between input values
	// (which, although of type complex, must be interpreted as real).

	std::ranges::generate(x, gen);
	std::ranges::generate(y, gen);

	std::valarray<std::complex<double>> z = a*x + b*y;

	// the last element (padding) will be ignored
	x = dft.rfft(x);
	y = dft.rfft(y);
	z = dft.rfft(z);

	assert(math::rms(z - (a*x + b*y)) < 1e-12);
}

void test_fft_ifft(std::size_t n)
// test that composing FFT and inverse FFT one obtains identity.
// `n` is the length of the sequence to test, which must be a power of 2.
// Starting from a vector x, the test passes if IFFT(FFT(x)) == x
// up to rounding errors.
{
	mersenne_twister.seed(1234);
	math::dft<double> dft;
	std::valarray<std::complex<double>> xref(n), x;

	std::ranges::generate(xref, gen);

	x = dft.fft(xref);
	x = dft.ifft(x);

	assert(math::rms(x-xref) < error_threshold);
}

void test_rfft_irfft(std::size_t n)
// test that composing FFT and inverse FFT one obtains identity.
// `n` is the length of the sequence to test, which must be a power of 2.
// Starting from a vector x, the test passes if IFFT(FFT(x)) == x
// up to rounding errors.
// Real-signal input variant.
{
	mersenne_twister.seed(1234);
	math::dft<double> dft;
	std::valarray<std::complex<double>> xref(n+1), x;

	std::ranges::generate(xref, gen);
	xref[n] = 0;

	x = dft.rfft(xref);
	x = dft.irfft(x);

	assert(math::rms(x-xref) < error_threshold);
}

void test_fft_ifft3(const std::array<std::size_t, 3>& n_shape)
// test that composing FFT and inverse FFT one obtains identity.
// `n_shape` is the shape of the 3-array to test, whose entries must be powers of 2.
// Starting from an array x, the test passes if IFFT(FFT(x)) == x
// up to rounding errors.
{
	mersenne_twister.seed(1234);
	std::size_t n = n_shape[0]*n_shape[1]*n_shape[2];

	math::dft<double> dft;
	std::valarray<std::complex<double>> xref(n), x;

	std::ranges::generate(xref, gen);

	x = xref;
	dft.fftn<3>(x, n_shape, tp);
	dft.ifftn<3>(x, n_shape, tp);

	assert(math::rms(x-xref) < error_threshold);
}

void test_rfft_irfft3(const std::array<std::size_t, 3>& n_shape)
// test that composing FFT and inverse FFT one obtains identity.
// `n_shape` is the shape of the 3-array to test, whose entries must be powers of 2
// (except the last dimension, which must be a power of 2 + 1).
// Starting from an array x, the test passes if IFFT(FFT(x)) == x
// up to rounding errors.
// Real-signal input variant.
{
	mersenne_twister.seed(1234);
	std::size_t n = n_shape[0]*n_shape[1]*n_shape[2];

	math::dft<double> dft;
	std::valarray<std::complex<double>> xref(n), x;

	std::ranges::generate(xref, gen);
	// setting padding elements to 0
	for (std::size_t i = n_shape[2]-1; i < n; i += n_shape[2])
		xref[i] = 0;

	x = xref;
	dft.rfftn<3>(x, n_shape, tp);
	dft.irfftn<3>(x, n_shape, tp);

	assert(math::rms(x-xref) < error_threshold);
}

void test_fft_linearity3(const std::array<std::size_t, 3>& n_shape)
// test linearity property of FFT up to rounding errors:
//	FFT(a x + b y) == a FFT(x) + b FFT(y)
// `n_shape` is the shape of the 3-array to test, whose entries must be powers of 2.
{
	mersenne_twister.seed(1234);
	std::size_t n = n_shape[0]*n_shape[1]*n_shape[2];

	math::dft<double> dft;
	std::valarray<std::complex<double>> x(n), y(n);
	std::complex<double> a = gen(), b = gen();

	std::ranges::generate(x, gen);
	std::ranges::generate(y, gen);

	std::valarray<std::complex<double>> z = a*x + b*y;

	dft.fftn<3>(x, n_shape, tp);
	dft.fftn<3>(y, n_shape, tp);
	dft.fftn<3>(z, n_shape, tp);

	assert(math::rms(z - (a*x + b*y)) < 1e-12);
}

void test_rfft_linearity3(const std::array<std::size_t, 3>& n_shape)
// test linearity property of FFT up to rounding errors:
//	FFT(a x + b y) == a FFT(x) + b FFT(y)
// `n_shape` is the shape of the 3-array to test, whose entries must be powers of 2
// (except the last dimension, which must be a power of 2 + 1).
// Real-signal input variant.
{
	mersenne_twister.seed(1234);
	std::size_t n = n_shape[0]*n_shape[1]*n_shape[2];

	math::dft<double> dft;
	std::valarray<std::complex<double>> x(n), y(n);
	std::complex<double> c = gen();
	std::complex<double> a = c.real(), b = c.imag();
	// a and b must are set to real to avoid mixing between input values
	// (which, although of type complex, must be interpreted as real).

	std::ranges::generate(x, gen);
	std::ranges::generate(y, gen);

	std::valarray<std::complex<double>> z = a*x + b*y;

	// the last element along the last dimension (padding) will be ignored
	dft.rfftn<3>(x, n_shape, tp);
	dft.rfftn<3>(y, n_shape, tp);
	dft.rfftn<3>(z, n_shape, tp);

	assert(math::rms(z - (a*x + b*y)) < 1e-12);
}

int main()
{
	test_gen();
	test_bit_reversal_incr(16);
	test_bit_reversal_permutation(8);
	test_twiddle_factors(1, 2);
	test_twiddle_factors((1<<20)+1, 1<<20);
	test_fft_impulse(1<<20, 2468);
	test_fft_impulse(1<<20, 3579);
	test_fft_impulse_real(1<<20, 4680);
	test_fft_impulse_real(1<<20, 1357);
	test_fft_linearity(1<<20);
	test_rfft_linearity(1<<20);
	test_fft_ifft(1<<20);
	test_rfft_irfft(1<<20);
	test_fft_ifft3({16, 64, 32});
	test_fft_ifft3({128, 128, 128});
	test_rfft_irfft3({16, 64, 33});
	test_rfft_irfft3({128, 128, 65});
	test_fft_linearity3({128, 128, 128});
	test_rfft_linearity3({128, 128, 65});

	std::cout << "All tests passed successfully!" << std::endl;

	return 0;
}
























