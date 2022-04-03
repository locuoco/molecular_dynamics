// VERTEX SHADER (POST)
#version 330 core

// Input (buffer)
layout(location = 0) in vec2 coord;

// Output (fragment shader)
out vec2 texCoord;

void main()
{
	vec2 pos = coord * 2 - 1;
	gl_Position = vec4(pos, 0, 1);

	texCoord = coord;
}





























