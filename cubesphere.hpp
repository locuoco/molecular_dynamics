#ifndef CUBESHPERE_H
#define CUBESHPERE_H

#include <cmath>

#include <GL/glew.h>
#include <glm/ext.hpp>

#include "math.hpp"

template <typename T>
T squareAngle(T l, T x, T y)
// ((x = l/2 or x = -l/2) and -l/2 <= y <= l/2) or ((y = l/2 or y = -l/2) and -l/2 <= x <= l/2)
// l = max(abs(x), abs(y))
{
	T res;
	if (y >= 0 && y <= x)
		res = y;
	else if (y > x && y >= -x)
		res = l - x;
	else if (y >= x && y < -x)
		res = 2*l - y;
	else if (y < x && y <= -x)
		res = 3*l + x;
	else // if (y < 0 && y > -x)
		res = 4*l + y;
	return res;
}

template <typename T>
T squareAngle(T l, const glm::tvec2<T>& v)
{
	return squareAngle(l, v.x, v.y);
}

template <typename T>
T squareAngle(T x, T y)
{
	T l = std::max(std::abs(x), std::abs(y));
	return squareAngle(l, x, y);
}

template <typename T>
T squareAngle(const glm::tvec2<T>& v)
{
	return squareAngle(v.x, v.y);
}

template <typename T>
inline glm::tvec2<T> transformSquareToCircle(T l, const glm::tvec2<T>& v)
{
	using namespace math;
	if (l == 0)
		return glm::tvec2<T>(0);
	T angle = half_pi<T>() * squareAngle(l, v) / l;
	T r = l / 2;
	return glm::tvec2<T>(r * std::cos(angle), r * std::sin(angle));
}

template <typename T>
glm::tvec2<T> transformSquareToCircle(const glm::tvec2<T>& v)
{
	T l = 2*std::max(std::abs(v.x), std::abs(v.y));
	return transformSquareToCircle(l, v);
}

template <typename T>
T cubicInclinationAngle(T l0, T l, T, T, T z)
{
	using namespace math;
	if (z >= l/2 * one_minus_eps<T>())
		return l0/2;
	else if (z <= -l/2 * one_minus_eps<T>())
		return 2*l - l0/2;
	else
		return l - z;
}

template <typename T>
T cubicInclinationAngle(T l0, T l, const glm::tvec3<T> &v)
{
	return cubicInclinationAngle(l0, l, v.x, v.y, v.z);
}

template <typename T>
T cubicInclinationAngle(T l, T x, T y, T z)
{
	return cubicInclinationAngle(2*std::max(std::abs(x), std::abs(y)), l, x, y, z);
}

template <typename T>
T cubicInclinationAngle(T l, const glm::tvec3<T> &v)
{
	return cubicInclinationAngle(l, v.x, v.y, v.z);
}

template <typename T>
T cubicInclinationAngle(T x, T y, T z)
{
	T l0 = 2*std::max(std::abs(x), std::abs(y));
	T l = std::max(l0, 2*std::abs(z));
	return cubicInclinationAngle(l0, l, x, y, z);
}

template <typename T>
T cubicInclinationAngle(const glm::tvec3<T> &v)
{
	return cubicInclinationAngle(v.x, v.y, v.z);
}

template <typename T>
glm::tvec3<T> transformCubeToSphere(T l0, T l, const glm::tvec3<T> &v)
{
	using namespace math;
	if (l == 0)
		return glm::vec3(0);
	T inclinationAngle = half_pi<T>() * cubicInclinationAngle(l0, l, v) / l;
	inclinationAngle = (inclinationAngle < 0) ? 0 :
					   (inclinationAngle > pi<T>()) ? pi<T>() : inclinationAngle;
	T azimuthAngle = (l0 == 0) ? 0 : half_pi<T>() * squareAngle(l0, v.x, v.y) / l0;
	T r = l / 2, sinInclination = std::sin(inclinationAngle);
	return glm::vec3(
		r * sinInclination * std::cos(azimuthAngle),
		r * sinInclination * std::sin(azimuthAngle),
		r * std::cos(inclinationAngle)
	);
}

template <typename T>
glm::tvec3<T> transformCubeToSphere(T l, const glm::tvec3<T> &v)
{
	T l0 = 2*std::max(std::abs(v.x), std::abs(v.y));
	return transformCubeToSphere(l0, l, v);
}

template <typename T>
glm::tvec3<T> transformCubeToSphere(const glm::tvec3<T> &v)
{
	T l0 = 2*std::max(std::abs(v.x), std::abs(v.y));
	T l = std::max(l0, 2*std::abs(v.z));
	return transformCubeToSphere(l0, l, v);
}

template <typename T>
glm::tvec2<T> transformSquareToCircle2(const glm::tvec2<T>& v)
// equal areas
{
	return glm::tvec2<T>(
		v.x * std::sqrt((1 + math::sqrt2_<T>() * v.y) * (1 - math::sqrt2_<T>() * v.y)),
		v.y * std::sqrt((1 + math::sqrt2_<T>() * v.x) * (1 - math::sqrt2_<T>() * v.x))
	);
}

template <typename T>
glm::tvec3<T> transformCubeToSphere2(const glm::tvec3<T>& v)
{
	T x2 = v.x * v.x, y2 = v.y * v.y, z2 = v.z * v.z;
	return glm::tvec3<T>(
		v.x * std::sqrt(1 - (y2 + z2) / 2 + y2 * z2 / 3),
		v.y * std::sqrt(1 - (x2 + z2) / 2 + x2 * z2 / 3),
		v.z * std::sqrt(1 - (x2 + y2) / 2 + x2 * y2 / 3)
	);
}

template <typename T>
glm::tvec3<T> transformCubeToSphereProj(T l, const glm::tvec3<T>& v)
{
	return (l/2) * glm::normalize(v);
}
template <typename T>
glm::tvec3<T> transformCubeToSphereProj(const glm::tvec3<T>& v)
{
	T l = 2*std::max(std::max(std::abs(v.x), std::abs(v.y)), std::abs(v.z));
	return transformCubeToSphereProj(l, v);
}

template <typename T>
glm::tvec3<T> transformSphereToCubeProj(T l, const glm::tvec3<T>& v)
{
	T l0 = 2*std::max(std::max(std::abs(v.x), std::abs(v.y)), std::abs(v.z));
	if (l0 == 0)
		return glm::tvec3<T>(0);
	T rescale = l / l0;
	return rescale * v;
}
template <typename T>
glm::tvec3<T> transformSphereToCubeProj(const glm::tvec3<T>& v)
{
	T l = 2*glm::length(v);
	return transformSphereToCubeProj(l, v);
}

template <typename T>
glm::tvec3<T> transformCubeToSphereTan(T l, const glm::tvec3<T>& v)
// equal angles
{
	if (l == 0)
		return glm::tvec3<T>(0);
	return transformCubeToSphereProj(l, tan(v * (math::half_pi<T>() / l)));
}
template <typename T>
glm::tvec3<T> transformCubeToSphereTan(const glm::tvec3<T>& v)
{
	T l = 2*std::max(std::max(std::abs(v.x), std::abs(v.y)), std::abs(v.z));
	return transformCubeToSphereTan(l, v);
}

template <typename T>
glm::tvec3<T> transformSphereToCubeTan(T l, const glm::tvec3<T>& v)
{
	return (l / math::half_pi<T>()) * atan(transformSphereToCubeProj(l, v));
}
template <typename T>
glm::tvec3<T> transformSphereToCubeTan(const glm::tvec3<T>& v)
{
	T l = 2*glm::length(v);
	return (l / math::half_pi<T>()) * atan(transformSphereToCubeProj(l, v));
}

inline GLuint genSphere(GLfloat *verts, const glm::vec3 &center, GLfloat radius, GLuint div)
{
	GLuint num = genCubeGrid(verts, center, 2*radius, div);

	if (verts != nullptr)
		for (unsigned int i = 0; i < num/3; ++i)
		{
			glm::vec3 v(verts[i*3 + 0], verts[i*3 + 1], verts[i*3 + 2]);
			v = transformCubeToSphereTan(v);
			verts[i*3 + 0] = v.x;
			verts[i*3 + 1] = v.y;
			verts[i*3 + 2] = v.z;
		}
	return num;
}

inline void genSphereNormals(GLfloat *normals, const GLfloat *verts, const glm::vec3 &center, GLuint n)
{
	for (GLuint i = 0; i < n / 3; ++i)
	{
		glm::vec3 v(verts[i*3 + 0], verts[i*3 + 1], verts[i*3 + 2]);

		glm::vec3 normal = glm::normalize(v - center);

		normals[i*3 + 0] = normal.x;
		normals[i*3 + 1] = normal.y;
		normals[i*3 + 2] = normal.z;
	}
}

#endif




















