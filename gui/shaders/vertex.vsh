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

// Jmol atomic CPK color scheme
const vec4 atom_pars[56] = vec4[](
	vec4(1.00, .000, .000, 0.00), //  0 default (red point)
	vec4(1.00, 1.00, 1.00, 1.20), //  1 hydrogen (white)
	vec4(.851, 1.00, 1.00, 1.40), //  2 helium
	vec4(.800, .502, 1.00, 1.82), //  3 lithium
	vec4(.761, 1.00, .000, 1.53), //  4 beryllium
	vec4(1.00, .710, .710, 1.92), //  5 boron
	vec4(.565, .565, .565, 1.70), //  6 carbon (dark gray)
	vec4(.188, .313, .973, 1.55), //  7 nitrogen (blue)
	vec4(1.00, .051, .051, 1.52), //  8 oxygen (red)
	vec4(.565, .878, .313, 1.47), //  9 fluorine
	vec4(.702, .890, .961, 1.54), // 10 neon
	vec4(.671, .361, .949, 2.27), // 11 sodium
	vec4(.541, 1.00, .000, 1.73), // 12 magnesium
	vec4(.749, .651, .651, 1.84), // 13 aluminium
	vec4(.941, .784, .627, 2.10), // 14 silicon
	vec4(1.00, .502, .000, 1.80), // 15 phosphorus (orange)
	vec4(1.00, 1.00, .188, 1.80), // 16 sulfur (yellow)
	vec4(.121, .941, .121, 1.75), // 17 chlorine (green)
	vec4(.502, .820, .890, 1.88), // 18 argon
	vec4(.561, .251, .831, 2.75), // 19 potassium
	vec4(.239, 1.00, .000, 2.31), // 20 calcium
	vec4(.902, .902, .902, 2.11), // 21 scandium
	vec4(.749, .761, .780, 2.04), // 22 titanium (interpolated radius)
	vec4(.651, .651, .671, 1.97), // 23 vanadium (interpolated radius)
	vec4(.541, .600, .780, 1.90), // 24 chromium (interpolated radius)
	vec4(.612, .478, .780, 1.84), // 25 manganese (interpolated radius)
	vec4(.878, .400, .200, 1.77), // 26 iron (interpolated radius)
	vec4(.941, .565, .627, 1.70), // 27 cobalt (interpolated radius)
	vec4(.314, .816, .314, 1.63), // 28 nichel
	vec4(.784, .502, .200, 1.40), // 29 copper
	vec4(.490, .502, .690, 1.39), // 30 zinc
	vec4(.761, .561, .561, 1.87), // 31 gallium
	vec4(.400, .561, .561, 2.11), // 32 germanium
	vec4(.741, .502, .890, 1.85), // 33 arsenic
	vec4(1.00, .631, .000, 1.90), // 34 selenium
	vec4(.651, .160, .160, 1.85), // 35 bromine
	vec4(.361, .722, .820, 2.02), // 36 krypton
	vec4(.439, .180, .690, 3.03), // 37 rubidium
	vec4(.000, 1.00, .000, 2.49), // 38 strontium
	vec4(.580, 1.00, 1.00, 2.38), // 39 yttrium (interpolated radius)
	vec4(.580, .878, .878, 2.28), // 40 zirconium (interpolated radius)
	vec4(.451, .761, .788, 2.17), // 41 niobium (interpolated radius)
	vec4(.329, .710, .710, 2.06), // 42 molybdenum (interpolated radius)
	vec4(.231, .620, .620, 1.95), // 43 technetium (interpolated radius)
	vec4(.141, .561, .561, 1.85), // 44 ruthenium (interpolated radius)
	vec4(.039, .490, .549, 1.74), // 45 rhodium (interpolated radius)
	vec4(.000, .412, .522, 1.63), // 46 palladium
	vec4(.753, .753, .753, 1.72), // 47 silver
	vec4(1.00, .850, .561, 1.58), // 48 cadmium
	vec4(.651, .459, .451, 1.93), // 49 indium
	vec4(.400, .502, .502, 2.17), // 50 tin
	vec4(.620, .388, .710, 2.06), // 51 antimony
	vec4(.831, .478, .000, 2.06), // 52 tellurium
	vec4(.580, .000, .580, 1.98), // 53 iodine
	vec4(.259, .620, .690, 2.16), // 54 xenon
	vec4(.341, .090, .561, 3.43)  // 55 caesium
);

void main()
{
	vec4 pos0 = MV * (vec4(inst_pos.xyz, 1));
	fragCenter = pos0.xyz; // position of the center of the sphere (square) in view space
	int atom_id = int(inst_pos.w);
	if (atom_id < 56)
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
	fragSize *= 0.75;
	float boxScal = fragSize*max(1., 1 / (.1 + .1*length(fragCenter)));
	vec3 zAxis = normalize(-fragCenter);
	vec3 xAxis = normalize(cross(vec3(0, 1, 0), zAxis));
	vec3 yAxis = cross(zAxis, xAxis);
	mat3 lookAtCamera = mat3(xAxis, yAxis, zAxis);
	pos0 += vec4(lookAtCamera*(vec3(pos*boxScal, 0)), 0);

	fragPos = pos0.xyz; // position of the square impostor vertex in view space

	gl_Position = Proj * pos0;
}























