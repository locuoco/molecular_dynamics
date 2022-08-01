//  Shaders loaders
//  Copyright (C) 2021-2022 Alessandro Lo Cuoco (alessandro.locuoco@gmail.com)

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

#ifndef GUI_SHADER_H
#define GUI_SHADER_H

#include <fstream>
#include <sstream>
#include <string>
#include <exception> // runtime_error

#include <GL/glew.h>

std::string load_code(const char* shader_path)
// Load shader code from file.
// `shader_path` is the path of the file.
// Return the shader code.
// Throw a `std::runtime_error` if the file does not exist.
{
	std::string shader_code;
	std::ifstream shstream(shader_path, std::ios::in);
	if (shstream)
	{
		std::stringstream sstr;
		sstr << shstream.rdbuf();
		shader_code = sstr.str();
		shstream.close();
	}
	else
		throw std::runtime_error(std::string("Error: Cannot open ") + shader_path + ".\n"
			"Check that the root folder of the repository has been set as the working "
			"directory or added to the PATH environment variable.");
	return shader_code;
}

void compile_code(GLuint shader, const std::string &code)
// Compile shader code.
// `shader` is the shader ID (returned by `glCreateShader`).
// `code` is the shader code extracted with `load_code`.
// Throw a `std::runtime_error` if there is a compilation error.
{
	GLint st;
	int n;

	const char *c_code = code.c_str();
	glShaderSource(shader, 1, &c_code, nullptr);
	glCompileShader(shader);

	glGetShaderiv(shader, GL_COMPILE_STATUS, &st);
	glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &n);
	if (n > 0)
	{
		char* mess = new char[n + 1];
		glGetShaderInfoLog(shader, n + 1, nullptr, mess);
		std::clog << mess << std::endl;
		delete[] mess;
	}
	if (!st)
		throw std::runtime_error("Error during compilation!");
}

void program_log(GLuint prog)
// Print program link log and check program status.
// `prog` is the program ID (returned by `glCreateProgram`).
// Throw a `std::runtime_error` if there is a program linking error.
{
	GLint st;
	int n;

	glGetProgramiv(prog, GL_LINK_STATUS, &st);
	glGetProgramiv(prog, GL_INFO_LOG_LENGTH, &n);
	if (n > 0)
	{
		char* mess = new char[n + 1];
		glGetProgramInfoLog(prog, n + 1, nullptr, mess);
		std::clog << mess << std::endl;
		delete[] mess;
	}
	if (!st)
		throw std::runtime_error("Error with linking program :(");
}

GLuint load_shader(const char* vs_path, const char* fs_path)
// Load vertex and fragment shaders.
// `vs_path` is the path of the vertex shader.
// `fs_path` is the path of the fragment shader.
// Return the program ID.
{
	GLuint vs = glCreateShader(GL_VERTEX_SHADER);
	GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);

	compile_code(vs, load_code(vs_path));
	compile_code(fs, load_code(fs_path));

	GLuint prog = glCreateProgram();
	glAttachShader(prog, vs);
	glAttachShader(prog, fs);
	glLinkProgram(prog);

	program_log(prog);

	glDetachShader(prog, vs);
	glDetachShader(prog, fs);

	glDeleteShader(vs);
	glDeleteShader(fs);

	return prog;
}

GLuint load_shader_tessellation(const char* vs_path, const char* tcs_path, const char* tes_path, const char* fs_path)
// Load vertex, tesselation and fragment shaders.
// `vs_path` is the path of the vertex shader.
// `tcs_path` is the path of the tesselation control shader.
// `tes_path` is the path of the tesselation evaluation shader.
// `fs_path` is the path of the fragment shader.
// Return the program ID.
{
	GLuint vs = glCreateShader(GL_VERTEX_SHADER);
	GLuint tcs = glCreateShader(GL_TESS_CONTROL_SHADER);
	GLuint tes = glCreateShader(GL_TESS_EVALUATION_SHADER);
	GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);

	compile_code(vs, load_code(vs_path));
	compile_code(tcs, load_code(tcs_path));
	compile_code(tes, load_code(tes_path));
	compile_code(fs, load_code(fs_path));

	GLuint prog = glCreateProgram();
	glAttachShader(prog, vs);
	glAttachShader(prog, tcs);
	glAttachShader(prog, tes);
	glAttachShader(prog, fs);
	glLinkProgram(prog);

	program_log(prog);

	glDetachShader(prog, vs);
	glDetachShader(prog, tcs);
	glDetachShader(prog, tes);
	glDetachShader(prog, fs);

	glDeleteShader(vs);
	glDeleteShader(tcs);
	glDeleteShader(tes);
	glDeleteShader(fs);

	return prog;
}

#endif // GUI_SHADER_H












