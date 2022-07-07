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

#ifndef SHADER_H
#define SHADER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <GL/glew.h>

std::string loadCode(const char* shaderPath)
{
	std::string shaderCode;
	std::ifstream shstream(shaderPath, std::ios::in);
	if (shstream)
	{
		std::stringstream sstr;
		sstr << shstream.rdbuf();
		shaderCode = sstr.str();
		shstream.close();
	}
	else
	{
		std::cerr << "Cannot open " << shaderPath << '.' << std::endl;
		throw;
	}
	return shaderCode;
}

void compileCode(GLuint shader, const std::string &code)
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
}

void programLog(GLuint prog)
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
}

GLuint loadShader(const char* vsPath, const char* fsPath)
{
	GLuint vs = glCreateShader(GL_VERTEX_SHADER);
	GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);

	std::string vsCode = loadCode(vsPath);
	std::string fsCode = loadCode(fsPath);

	compileCode(vs, vsCode);
	compileCode(fs, fsCode);

	GLuint prog = glCreateProgram();
	glAttachShader(prog, vs);
	glAttachShader(prog, fs);
	glLinkProgram(prog);

	programLog(prog);

	glDetachShader(prog, vs);
	glDetachShader(prog, fs);

	glDeleteShader(vs);
	glDeleteShader(fs);

	return prog;
}

GLuint loadShaderTessellation(const char* vsPath, const char* tcsPath, const char* tesPath, const char* fsPath)
{
	GLuint vs = glCreateShader(GL_VERTEX_SHADER);
	GLuint tcs = glCreateShader(GL_TESS_CONTROL_SHADER);
	GLuint tes = glCreateShader(GL_TESS_EVALUATION_SHADER);
	GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);

	std::string vsCode = loadCode(vsPath);
	std::string tcsCode = loadCode(tcsPath);
	std::string tesCode = loadCode(tesPath);
	std::string fsCode = loadCode(fsPath);

	compileCode(vs, vsCode);
	compileCode(tcs, tcsCode);
	compileCode(tes, tesCode);
	compileCode(fs, fsCode);

	GLuint prog = glCreateProgram();
	glAttachShader(prog, vs);
	glAttachShader(prog, tcs);
	glAttachShader(prog, tes);
	glAttachShader(prog, fs);
	glLinkProgram(prog);

	programLog(prog);

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

#endif












