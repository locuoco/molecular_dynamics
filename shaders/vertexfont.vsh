//  Vertex shader for text/font rendering
//  Copyright (C) 2018-2022 Alessandro Lo Cuoco (alessandro.locuoco@gmail.com)

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