//  Fragment shader for text/font rendering
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

// Input (vertex shader & uniform)
in vec2 texcoord;
uniform sampler2D tex;
uniform vec4 col;

// Output (schermo)
out vec4 color;

void main()
{
	// outputs the pixel color and transparency
	// texture `tex` contains only the alpha channel of the character to print
	// `col` is the color of the text
	color = vec4(1, 1, 1, texture(tex, texcoord).r) * col;
}