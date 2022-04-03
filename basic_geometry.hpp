#ifndef BASIC_GEOMETRY_H
#define BASIC_GEOMETRY_H

#include <cmath>

#include <GL/glew.h>
#include <glm/ext.hpp>

#include "math.hpp"

struct vertex
{
	GLfloat position[3];
	GLfloat normal[3];
	GLfloat color[3];

	vertex() = default;

	vertex(const glm::vec3 &p, const glm::vec3 &n, const glm::vec3 &c)
	{
		position[0] = p.x;
		position[1] = p.y;
		position[2] = p.z;
		normal[0] = n.x;
		normal[1] = n.y;
		normal[2] = n.z;
		color[0] = c.x;
		color[1] = c.y;
		color[2] = c.z;
	}

	static void format(GLint vbo)
	{
		glEnableVertexAttribArray(0);
		glEnableVertexAttribArray(1);
		glEnableVertexAttribArray(2);

		glVertexAttribFormat(0, sizeof(position)/sizeof(GLfloat), GL_FLOAT, GL_FALSE, 0);
		glVertexAttribFormat(1, sizeof(normal)/sizeof(GLfloat), GL_FLOAT, GL_FALSE, 0);
		glVertexAttribFormat(2, sizeof(color)/sizeof(GLfloat), GL_FLOAT, GL_FALSE, 0);

		glVertexBindingDivisor(0, 0);
		glVertexBindingDivisor(1, 0);
		glVertexBindingDivisor(2, 0);

		glVertexAttribBinding(0, 0);
		glVertexAttribBinding(1, 1);
		glVertexAttribBinding(2, 2);

		glBindVertexBuffer(0, vbo, offsetof(vertex, position), sizeof(vertex));
		glBindVertexBuffer(1, vbo, offsetof(vertex, normal), sizeof(vertex));
		glBindVertexBuffer(2, vbo, offsetof(vertex, color), sizeof(vertex));
	}

	static void unformat()
	{
		glDisableVertexAttribArray(0);
		glDisableVertexAttribArray(1);
		glDisableVertexAttribArray(2);
	}
};

inline GLuint genTriangle2D(GLfloat *verts, const glm::vec2 &v0, const glm::vec2 &v1, const glm::vec2 &v2)
// returns number of GLfloats needed in verts
// preserve triangle orientation
{
	if (verts != nullptr)
	{
		verts[0] = v0.x;
		verts[1] = v0.y;
		verts[2] = v1.x;
		verts[3] = v1.y;
		verts[4] = v2.x;
		verts[5] = v2.y;
	}
	return 3*2; // 3 vertices, 2 coordinates per vertex
}

inline GLuint genTriangle3D(GLfloat *verts, const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &v2)
// returns number of GLfloats needed in verts
// preserve triangle orientation
{
	if (verts != nullptr)
	{
		verts[0] = v0.x;
		verts[1] = v0.y;
		verts[2] = v0.z;
		verts[3] = v1.x;
		verts[4] = v1.y;
		verts[5] = v1.z;
		verts[6] = v2.x;
		verts[7] = v2.y;
		verts[8] = v2.z;
	}
	return 3*3; // 3 vertices, 3 coordinates per vertex
}

inline GLuint genQuad2D(GLfloat *verts, const glm::vec2 &v0, const glm::vec2 &v1, const glm::vec2 &v2, const glm::vec2 &v3)
// returns number of GLfloats needed in verts
// preserve orientation
// vertices must be given in clockwise or counter-clockwise order
/*
v3		  v2
+---------*
|      /  |
|    /    |
|  /      |
*---------+
v0        v1
*/
{
	GLuint num = genTriangle2D(verts, v0, v1, v2);
	if (verts != nullptr)
		verts += num;
	num += genTriangle2D(verts, v0, v2, v3);
	return num;
}

inline GLuint genQuad3D(GLfloat *verts, const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &v2, const glm::vec3 &v3)
// returns number of GLfloats needed in verts
// preserve orientation
// vertices must be given in clockwise or counter-clockwise order
{
	GLuint num = genTriangle3D(verts, v0, v1, v2);
	if (verts != nullptr)
		verts += num;
	num += genTriangle3D(verts, v0, v2, v3);
	return num;
}

inline GLuint genQuadGrid2D(GLfloat *verts, GLuint div, const glm::vec2 &v0, const glm::vec2 &v1, const glm::vec2 &v2, const glm::vec2 &v3)
{
	GLuint num = 0;
	if (verts != nullptr)
	{
		GLfloat *beg = verts;
		GLfloat div2_ = 1.f / (div*div);
		for (GLuint i = 0; i < div; ++i)
			for (GLuint j = 0; j < div; ++j)
			{
				num += genQuad2D(
					verts,
					(GLfloat((div-i)*(div-j)) * v0 + GLfloat(i*(div-j)) * v1 + GLfloat((div-i)*j) * v3 + GLfloat(i*j) * v2) * div2_,
					(GLfloat((div-i-1)*(div-j)) * v0 + GLfloat((i+1)*(div-j)) * v1 + GLfloat((div-i-1)*j) * v3 + GLfloat((i+1)*j) * v2) * div2_,
					(GLfloat((div-i-1)*(div-j-1)) * v0 + GLfloat((i+1)*(div-j-1)) * v1 + GLfloat((div-i-1)*(j+1)) * v3 + GLfloat((i+1)*(j+1)) * v2) * div2_,
					(GLfloat((div-i)*(div-j-1)) * v0 + GLfloat(i*(div-j-1)) * v1 + GLfloat((div-i)*(j+1)) * v3 + GLfloat(i*(j+1)) * v2) * div2_
				);
				verts = beg + num;
			}
	}
	else
		num = div * div * genQuad2D(nullptr, v0, v1, v2, v3);
	return num;
}

inline GLuint genQuadGrid3D(GLfloat *verts, GLuint div, const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &v2, const glm::vec3 &v3)
{
	GLuint num = 0;
	if (verts != nullptr)
	{
		GLfloat *beg = verts;
		GLfloat div2_ = 1.f / (div*div);
		for (GLuint i = 0; i < div; ++i)
			for (GLuint j = 0; j < div; ++j)
			{
				num += genQuad3D(
					verts,
					(GLfloat((div-i)*(div-j)) * v0 + GLfloat(i*(div-j)) * v1 + GLfloat((div-i)*j) * v3 + GLfloat(i*j) * v2) * div2_,
					(GLfloat((div-i-1)*(div-j)) * v0 + GLfloat((i+1)*(div-j)) * v1 + GLfloat((div-i-1)*j) * v3 + GLfloat((i+1)*j) * v2) * div2_,
					(GLfloat((div-i-1)*(div-j-1)) * v0 + GLfloat((i+1)*(div-j-1)) * v1 + GLfloat((div-i-1))*(j+1) * v3 + GLfloat((i+1)*(j+1)) * v2) * div2_,
					(GLfloat((div-i)*(div-j-1)) * v0 + GLfloat(i*(div-j-1)) * v1 + GLfloat((div-i)*(j+1)) * v3 + GLfloat(i*(j+1)) * v2) * div2_
				);
				verts = beg + num;
			}
	}
	else
		num = div * div * genQuad3D(nullptr, v0, v1, v2, v3);
	return num;
}

inline GLuint genCube(GLfloat *verts, const glm::vec3 &center, GLfloat side)
// returns number of GLfloats needed in verts
// vertices are in counter-clockwise order
// just make a reflection through any plane to have a clockwise order
{
	GLuint num = 0;
	if (verts != nullptr)
	{
		GLfloat s = side/2;
		GLfloat *beg = verts;
		glm::vec3 v[8]
		{
			center + glm::vec3(-s, -s, s), // 0
			center + glm::vec3(s, -s, s), // 1
			center + glm::vec3(s, s, s), // 2
			center + glm::vec3(-s, s, s), // 3
			center + glm::vec3(s, -s, -s), // 4
			center + glm::vec3(-s, -s, -s), // 5
			center + glm::vec3(-s, s, -s), // 6
			center + glm::vec3(s, s, -s) // 7
		};
		// front
		num += genQuad3D(verts, v[0], v[1], v[2], v[3]);
		verts = beg + num;
		// right
		num += genQuad3D(verts, v[1], v[4], v[7], v[2]);
		verts = beg + num;
		// back
		num += genQuad3D(verts, v[4], v[5], v[6], v[7]);
		verts = beg + num;
		// left
		num += genQuad3D(verts, v[5], v[0], v[3], v[6]);
		verts = beg + num;
		// top
		num += genQuad3D(verts, v[3], v[2], v[7], v[6]);
		verts = beg + num;
		// bottom
		num += genQuad3D(verts, v[5], v[4], v[1], v[0]);
	}
	else
		num = 6 * genQuad3D(nullptr, center, center, center, center);
	
	return num;
}

inline GLuint genCubeGrid(GLfloat *verts, const glm::vec3 &center, GLfloat side, GLuint div)
// returns number of GLfloats needed in verts
// vertices are in counter-clockwise order
// just make a reflection through any plane to have a clockwise order
{
	GLuint num = 0;
	if (verts != nullptr)
	{
		GLfloat s = side/2;
		GLfloat *beg = verts;
		glm::vec3 v[8]
		{
			center + glm::vec3(-s, -s, s), // 0
			center + glm::vec3(s, -s, s), // 1
			center + glm::vec3(s, s, s), // 2
			center + glm::vec3(-s, s, s), // 3
			center + glm::vec3(s, -s, -s), // 4
			center + glm::vec3(-s, -s, -s), // 5
			center + glm::vec3(-s, s, -s), // 6
			center + glm::vec3(s, s, -s) // 7
		};
		// front
		num += genQuadGrid3D(verts, div, v[0], v[1], v[2], v[3]);
		verts = beg + num;
		// right
		num += genQuadGrid3D(verts, div, v[1], v[4], v[7], v[2]);
		verts = beg + num;
		// back
		num += genQuadGrid3D(verts, div, v[4], v[5], v[6], v[7]);
		verts = beg + num;
		// left
		num += genQuadGrid3D(verts, div, v[5], v[0], v[3], v[6]);
		verts = beg + num;
		// top
		num += genQuadGrid3D(verts, div, v[3], v[2], v[7], v[6]);
		verts = beg + num;
		// bottom
		num += genQuadGrid3D(verts, div, v[5], v[4], v[1], v[0]);
	}
	else
		num = 6 * genQuadGrid3D(nullptr, div, center, center, center, center);
	
	return num;
}

inline void genNormals(GLfloat *normals, const GLfloat *verts, GLuint n)
{
	for (GLuint i = 0; i < n / 9; ++i)
	{
		glm::vec3 v0(verts[i*9 + 0], verts[i*9 + 1], verts[i*9 + 2]);
		glm::vec3 v1(verts[i*9 + 3], verts[i*9 + 4], verts[i*9 + 5]);
		glm::vec3 v2(verts[i*9 + 6], verts[i*9 + 7], verts[i*9 + 8]);

		glm::vec3 normal = glm::normalize(glm::cross(v1 - v0, v2 - v0));

		normals[i*9 + 0] = normal.x;
		normals[i*9 + 1] = normal.y;
		normals[i*9 + 2] = normal.z;
		normals[i*9 + 3] = normal.x;
		normals[i*9 + 4] = normal.y;
		normals[i*9 + 5] = normal.z;
		normals[i*9 + 6] = normal.x;
		normals[i*9 + 7] = normal.y;
		normals[i*9 + 8] = normal.z;
	}
}

#endif




















