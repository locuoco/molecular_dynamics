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
			"Check that the working directory is set to the root folder of the repository.");
	return shader_code;
}

void compile_code(GLuint shader, const std::string &code)
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
		throw std::runtime_error(std::string("Error in compiling:\n") + code + '\n');
}

void program_log(GLuint prog)
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
		throw std::runtime_error("Error with linking program :(\n");
}

GLuint load_shader(const char* vs_path, const char* fs_path)
{
	GLuint vs = glCreateShader(GL_VERTEX_SHADER);
	GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);

	std::string vs_code = load_code(vs_path);
	std::string fs_code = load_code(fs_path);

	compile_code(vs, vs_code);
	compile_code(fs, fs_code);

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
{
	GLuint vs = glCreateShader(GL_VERTEX_SHADER);
	GLuint tcs = glCreateShader(GL_TESS_CONTROL_SHADER);
	GLuint tes = glCreateShader(GL_TESS_EVALUATION_SHADER);
	GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);

	std::string vs_code = load_code(vs_path);
	std::string tcs_code = load_code(tcs_path);
	std::string tes_code = load_code(tes_path);
	std::string fs_code = load_code(fs_path);

	compile_code(vs, vs_code);
	compile_code(tcs, tcs_code);
	compile_code(tes, tes_code);
	compile_code(fs, fs_code);

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












