//  Fragment shader for main 3D graphics (imposters raytracing, lighting,
//    tone maps, gamma correction etc...)
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

#version 330 core

// Input from vertex shader (interpolated)
in vec3 fragPos;
flat in vec3 fragCenter;
flat in vec3 fragCol;
flat in float fragSize;

// Values that stay constant for all pixels
uniform mat4 Proj;
uniform vec3 lightPos;
uniform float ambient = 0.01; //0.001
uniform float specularStrength = 1;
uniform float shininess = 16;
uniform float gamma = 2.2 + 1./30;

const float pi = 3.14159265;
const float sphRadius = 1;

// ACES tone map

vec3 aces(vec3 x)
{
	mat3 m1 = mat3(
		0.59719, 0.07600, 0.02840,
		0.35458, 0.90834, 0.13383,
		0.04823, 0.01566, 0.83777
	);
	mat3 m2 = mat3(
		 1.60475, -0.10208, -0.00327,
		-0.53108,  1.10813, -0.07276,
		-0.07367, -0.00605,  1.07602
	);
	vec3 v = m1 * x;
	vec3 a = v * (v + 0.0245786) - 0.000090537;
	vec3 b = v * (0.983729 * v + 0.4329510) + 0.238081;
	return clamp(m2 * (a / b), 0., 1.);
}

vec3 aces_approx(vec3 x)
{
	const float a = 2.51, b = 0.03, c = 2.43, d = 0.59, e = 0.14;
	return clamp((x * (a*x + b)) / (x * (c*x + d) + e), 0., 1.);
}

// PBR

float NormalDistributionFunction(vec3 n, vec3 h, float a)
// n: surface normal
// h: halfway vector
// a: roughness
{
	float a2 = a*a;
	float ndoth = max(dot(n, h), 0);
	float div = ndoth*ndoth * (a2 - 1) + 1;
	return a2 / (pi * div*div);
}

float GeometrySchlickBeckmann(float ndotv, float k)
// ndotv: dot product between normal vector and another direction
// k: roughness parameter (related to a)
{
	return ndotv / (ndotv*(1-k) + k);
}

float GeometryFunction(float ndotv, float ndotl, float a)
// ndotv: dot product between normal vector and view direction
// ndotl: dot product between normal vector and light direction
// a: roughness
{
	float a1 = a + 1;
	float k = a1*a1 / 8;
	return GeometrySchlickBeckmann(ndotv, k) * GeometrySchlickBeckmann(ndotl, k);
}

vec3 FresnelEquation(inout vec3 h, inout vec3 v, inout vec3 F0)
// h: halfway vector
// v: view direction vector
// F0: base reflectivity
{
	float hdotv = max(dot(h, v), 0);
	return F0 + (1 - F0) * pow(1 - hdotv, 5);
}

// Output (screen)
layout(location = 0) out vec4 color;

void main()
{
	// r2 is the square distance of sphere from the viewer
	float r2 = dot(fragCenter, fragCenter);
	// R2 is the square radius of the sphere to render
	float R2 = fragSize * fragSize;
	if (r2 < R2)
		discard; // the surface is behind the near plane
	// the ray direction is taken by normalizing the fragment position
	// relative to the viewer
	vec3 rayDir = normalize(fragPos);
	// calculate the position of the sphere point
	float B = dot(rayDir, fragCenter);
	float C = r2 - R2;
	float Discriminant = B*B - C;
	if (Discriminant < 0)
		discard; // the fragment of the square is outside the sphere
	float t = B - sqrt(Discriminant); // point of the sphere nearest to the point of view
	vec3 realPos = rayDir * t;
	vec3 realNorm = realPos - fragCenter;

	vec4 clipPos = Proj * vec4(realPos, 1);
	gl_FragDepth = clipPos.z / clipPos.w; // explicitly set the depth buffer

	vec3 n = normalize(realNorm);
	vec3 v = -rayDir;

	vec3 lightx[2], lightcol[2];
	lightx[0] = lightPos;
	lightx[1] = -lightPos;
	lightcol[0] = 5*vec3(1, 1, 1);
	lightcol[1] = 5*vec3(.1, .1, .1);

	float metallic = .1;
	float roughness = .35;

	vec3 lo = vec3(0);
	vec3 albedo = fragCol;

	for (int i = 0; i < 2; ++i)
	{
		vec3 l = normalize(lightx[i]);
		vec3 h = normalize(v + l);
		vec3 radiance = lightcol[i];

		float ndotv = max(dot(n, v), 0);
		float ndotl = max(dot(n, l), 0);

		vec3 F0 = vec3(0.04);
		F0 = mix(F0, albedo, metallic);
		vec3 F = FresnelEquation(h, v, F0);

		float D = NormalDistributionFunction(n, h, roughness);
		float G = GeometryFunction(ndotv, ndotl, roughness);

		vec3 specular = D * F * G / (4 * ndotv * ndotl + 0.0001);

		vec3 kD = vec3(1) - F;
		kD *= 1 - metallic;

		lo += (kD * albedo / pi + specular) * radiance * ndotl;
	}
	
	// tone mapping + gamma correction
	color = vec4(pow(aces(lo), vec3(1/gamma)), 1);
}




















