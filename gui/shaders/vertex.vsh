//  Vertex shader for main 3D graphics (imposters rotation)
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

// Input vertex, same code will be called for each vertex
layout(location = 0) in vec2 pos;

layout(location = 10) in vec4 inst_pos;

// Values that stay constant for all vertices
uniform mat4 MV;
uniform mat4 Proj;
uniform float gamma = 2.2 + 1./30;

// Output (to fragment shader)
out vec3 fragPos;
flat out vec3 fragCenter;
flat out vec3 fragCol;
flat out float fragSize;

const vec4 atom_pars[19] = vec4[](
	vec4(1.00, 1.00, 1.00, 1.00), // default (white)
	vec4(1.00, 1.00, 1.00, 1.20), // hydrogen (white)
	vec4(.851, 1.00, 1.00, 1.40), // helium
	vec4(.800, .502, 1.00, 1.82), // lithium
	vec4(.761, 1.00, .000, 1.53), // beryllium
	vec4(1.00, .710, .710, 1.92), // boron
	vec4(.565, .565, .565, 1.70), // carbon (dark gray)
	vec4(.188, .313, .973, 1.55), // nitrogen (blue)
	vec4(1.00, .051, .051, 1.52), // oxygen (red)
	vec4(.565, .878, .313, 1.47), // fluorine
	vec4(.702, .890, .961, 1.54), // neon
	vec4(.671, .361, .949, 2.27), // sodium
	vec4(.541, 1.00, .000, 1.73), // magnesium
	vec4(.749, .651, .651, 1.84), // aluminium
	vec4(.941, .784, .627, 2.10), // silicon
	vec4(1.00, .502, .000, 1.80), // phosphorus (orange)
	vec4(1.00, 1.00, .188, 1.80), // sulfur (yellow)
	vec4(.121, .941, .121, 1.75), // chlorine (green)
	vec4(.502, .820, .890, 1.88)  // argon
);

void main()
{
	vec4 pos0 = MV * (vec4(inst_pos.xyz, 1)); // position of the center of the sphere (square) in view space
	fragCenter = pos0.xyz;
	int atom_id = int(inst_pos.w);
	if (atom_id < 19)
	{
		fragCol = atom_pars[atom_id].xyz;
		fragSize = atom_pars[atom_id].w;
	}
	else
	{
		// use default color and size
		fragCol = atom_pars[0].xyz;
		fragSize = atom_pars[0].w;
	}
	fragCol = pow(fragCol, vec3(gamma)); // gamma correction
	fragSize *= 0.5;
	float boxScal = fragSize*max(1., 1 / (.1 + .1*length(fragCenter)));
	vec3 zAxis = normalize(-fragCenter);
	vec3 xAxis = normalize(cross(vec3(0, 1, 0), zAxis));
	vec3 yAxis = cross(zAxis, xAxis);
	mat3 lookAtCamera = mat3(xAxis, yAxis, zAxis);
	pos0 += vec4(lookAtCamera*(vec3(pos*boxScal, 0)), 0); // position of the square impostor vertex in view space

	fragPos = pos0.xyz;

	gl_Position = Proj * pos0;
}























