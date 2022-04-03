// VERTEX SHADER (FONT)
#version 330 core

// Input (buffer)
layout(location = 0) in vec4 coord;
uniform mat4 Pmat;

// Output (fragment shader)
out vec2 texcoord;

void main()
{
	gl_Position = Pmat * vec4(coord.xy, 0, 1);
	texcoord = coord.zw;
}