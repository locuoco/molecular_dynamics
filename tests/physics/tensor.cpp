//  Tensor tests
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
#include <type_traits> // is_standard_layout_v, is_trivial_v
#include <numeric> // iota
#include <cassert>

/*

Compilation:
g++ tensor.cpp -o tensor -std=c++20 -Wall -Wextra -pedantic -Ofast -fmax-errors=1

*/

#include "../../physics/tensor.hpp"

void test_tensor_pod()
// test (at compile time) that physics::tensor<double, 3> is a POD (plain-old data).
// Since the assertion is checked at compile time and the function is already
// instantiated (it is not templated), the function does not need to be called.
{
	using type = physics::tensor<double, 1, 2, 3, 4>;
	static_assert(std::is_standard_layout_v<type> && std::is_trivial_v<type>, "tensor is not a POD!");
	// std::is_pod_v is deprecated in C++20 (it is the same as is_standard_layout_v && is_trivial_v)
}

void test_tensor_fill_constructor()
// test tensor fill constructor (using 42 as argument).
// The test is passed if all elements of the constructed tensor are equal to 42.
// Since the assertion is checked at compile time and the function is already
// instantiated (it is not templated), the function does not need to be called.
{
	constexpr physics::tensor<long, 3, 3, 3> t(42);
	static_assert(t.all(42));
}

void test_matrix_product()
// a 3x4 matrix filled with 2's multiplied by 4x3 filled with 3's should give
// a 3x3 matrix filled with 24's (i.e. 2*3*4).
// Since the assertion is checked at compile time and the function is already
// instantiated (it is not templated), the function does not need to be called.
{
	constexpr physics::mat<float, 3, 4> m2(2);
	constexpr physics::mat<float, 4, 3> m3(3);
	constexpr physics::mat3f m24(24);
	static_assert(m2 % m3 == m24);
	// operations with integers are exact even with floating points
}

void test_outer_product()
// unit test for outer product
// Since the assertion is checked at compile time and the function is already
// instantiated (it is not templated), the function does not need to be called.
{
	constexpr physics::vec3d v1{1,2,3}, v2{4,5,6};
	constexpr physics::mat3d m{
		 4, 5, 6,
		 8,10,12,
		12,15,18,
	};
	static_assert(outer(v1, v2) == m); // argument-dependent lookup
}

void test_dot_product()
// unit test for dot product
// Since the assertion is checked at compile time and the function is already
// instantiated (it is not templated), the function does not need to be called.
{
	constexpr physics::vec3d v1{1,2,3}, v2{4,5,6};
	static_assert(dot(v1, v2) == 32); // argument-dependent lookup
}

void test_trace()
// test that the trace of an identity matrix is equal to the number of rows.
// Since the assertion is checked at compile time and the function is already
// instantiated (it is not templated), the function does not need to be called.
{
	constexpr std::size_t N = 8;
	// class template argument deduction (since C++17)
	constexpr physics::mat m(physics::make_identity<double, N>());
	static_assert(trace(m) == N); // argument-dependent lookup
}

void test_cross_product()
// unit test for cross product
// Since the assertion is checked at compile time and the function is already
// instantiated (it is not templated), the function does not need to be called.
{
	constexpr physics::vec3d v1{1,2,3}, v2{4,5,6}, res{-3, 6, -3};
	static_assert(cross(v1, v2) == res); // argument-dependent lookup
}

void test_transpose()
// Check that a matrix and its transpose verify the `is_transpose` condition.
// Since the assertion is checked at compile time and the function is already
// instantiated (it is not templated), the function does not need to be called.
{
	constexpr auto init = []
		{
			physics::mat4ld m;
			// the following is valid because physics::tensor inherits from std::array,
			// and creates an increasing sequence in row-major ordering.
			std::iota(begin(m), end(m), 0);
			return m;
		};

	constexpr auto m(init()), mt(transpose(m));

	static_assert(is_transpose(m, mt));
}

int main()
{
	std::cout << "All tests have been passed successfully!" << std::endl;

	return 0;
}
























