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
uniform float specularStrength = 0.5;
uniform float shininess = 32;
uniform float gamma = 2.2 + 1./30;

const float pi = 3.14159265;
const float sphRadius = 1;

// ACES tone map

vec3 aces(vec3 x)
{
	const float a = 2.51, b = 0.03, c = 2.43, d = 0.59, e = 0.14;
	return clamp((x * (a*x + b)) / (x * (c*x + d) + e), 0., 1.);
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

	// calculate directional light (ambient, diffuse, specular)
	vec3 norm = normalize(realNorm);
	vec3 lightDir = normalize(lightPos);
	float cosangle = dot(norm, lightDir);
	float diffuse = max(cosangle, 0);
	float antidiffuse = -min(cosangle, 0);

	vec3 viewDir = -rayDir;
	vec3 halfwayDir = normalize(lightDir + viewDir);
	float cosangle2 = dot(norm, halfwayDir);
	float spec = max(cosangle2, 0);
	float antispec = -min(cosangle2, 0);
	float energyFactor = (8 + shininess) / (8 * pi);
	float specular = energyFactor * pow(spec, shininess) * specularStrength;
	float antispecular = energyFactor * pow(antispec, shininess) * specularStrength;

	float light = ambient + diffuse + specular + .1*antidiffuse + .1*antispecular;
	
	// tone mapping + gamma correction
	color = vec4(pow(aces(fragCol * light), vec3(1/gamma)), 1);
}




















