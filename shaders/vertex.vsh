// VERTEX SHADER
#version 330 core

// Input vertex, same code will be called for each vertex
layout(location = 0) in vec3 pos;

layout(location = 10) in vec4 inst_pos;

// Values that stay constant for all vertices
uniform mat4 MV;
uniform mat4 Proj;

// Output (to fragment shader)
out vec3 fragPos;
flat out vec3 fragCenter;
flat out vec3 fragCol;
flat out float fragSize;

void main()
{
	vec4 pos0 = MV * (vec4(inst_pos.xyz, 1)); // position of the center of the sphere (square) in view space
	fragCenter = pos0.xyz;
	int atom_id = int(inst_pos.w);
	switch (atom_id)
	{
		case 1: // hydrogen
			fragCol = vec3(1, 1, 1); // white
			fragSize = 1.20; // van der Waals radius
			break;
		case 6: // carbon
			fragCol = vec3(.1, .1, .1); // dark gray
			fragSize = 1.70;
			break;
		case 7: // nitrogen
			fragCol = vec3(0, 0, 1); // blue
			fragSize = 1.55;
			break;
		case 8: // oxygen
			fragCol = vec3(1, 0, 0); // red
			fragSize = 1.52;
			break;
		case 15: // phosphorus
			fragCol = vec3(1, .3, 0); // orange
			fragSize = 1.80;
			break;
		case 16: // sulfur
			fragCol = vec3(1, 1, 0); // yellow
			fragSize = 1.80;
			break;
		case 17: // chlorine
			fragCol = vec3(0, 1, 0); // green
			fragSize = 1.75;
			break;
		default:
			fragCol = vec3(1, 1, 1);
			fragSize = 1.;
			break;
	}
	fragSize *= 0.5;
	float boxScal = fragSize*max(1., 1 / (.1 + .1*length(fragCenter)));
	vec3 zAxis = normalize(-fragCenter);
	vec3 xAxis = normalize(cross(vec3(0, 1, 0), zAxis));
	vec3 yAxis = cross(zAxis, xAxis);
	mat3 lookAtCamera = mat3(xAxis, yAxis, zAxis);
	pos0 = vec4(lookAtCamera*(pos*boxScal), 0) + pos0; // position of the square impostor vertex in view space

	fragPos = pos0.xyz;

	gl_Position = Proj * pos0;
}























