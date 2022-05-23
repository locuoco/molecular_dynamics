//  Graphics for molecular dynamics
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

#ifndef GRAPHICS_H
#define GRAPHICS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <thread>
#include <vector>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <glm/ext.hpp>

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

#include "shader.hpp"
#include "controls.hpp"
#include "physics/molecule.hpp"

#define TEST(x)
#include "Font.hpp"

#ifdef __MINGW32__
#define EXPORT __declspec(dllexport)
#elif defined(_MSC_VER)
#define EXPORT _declspec(dllexport)
#elif defined(__GNUC__)
#define EXPORT __attribute__((visibility("default")))
#else
#define EXPORT
#endif

// enable optimus, which chooses dedicated GPU if present
extern "C"
{
    EXPORT unsigned long NvOptimusEnablement = 1;
    EXPORT int AmdPowerXpressRequestHighPerformance = 1;
}

class graphics
{
	public:

		graphics() : dt(0), frame(0), FPS(0)
		{
			glfwSetErrorCallback(&ErrorCallback);
			// OpenGL graphic interface initialization (GLFW + GLEW)
			if (!glfwInit())
			{
				std::cerr << "Error: Cannot initialize GLFW." << std::endl;
				throw;
			}
			//glfwWindowHint(GLFW_SAMPLES, 4); // for multisampling
			glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
			glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
			glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
			glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
			glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);

			vidmode = glfwGetVideoMode(glfwGetPrimaryMonitor());

			glfwGetMonitorContentScale(glfwGetPrimaryMonitor(), &xscale, &yscale);
			w = 1280*xscale, h = 720*yscale;

			window = glfwCreateWindow(w, h, "Molecular Dynamics", nullptr, nullptr);
			// glfwGetPrimaryMonitor() for fullscreen (as 4th argument)

			if (!window)
			{
				std::cerr << "GLFW Error: Cannot create window." << std::endl;
				throw;
			}

			glfwMakeContextCurrent(window);
			glfwSwapInterval(0);

			std::clog << "GPU: " << glGetString(GL_RENDERER) << std::endl;

			if (Init(w, h) == -1)
				throw;

			glfwPollEvents();

			glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);

			glfwSetScrollCallback(window, updateScroll);

			lastTime = glfwGetTime();
		}

		graphics(const graphics&) = delete;

		graphics& operator=(const graphics&) = delete;

		~graphics()
		{
			delete font;

			glDeleteProgram(progID);

			glDeleteFramebuffers(1, &fb);
			glDeleteTextures(1, &texScene);
			glDeleteTextures(1, &texBlueNoise);
			glDeleteRenderbuffers(1, &rb);
			glDeleteBuffers(1, &vbPost);
			glDeleteBuffers(1, &vb_atom);
			glDeleteVertexArrays(1, &va);

			glfwTerminate();
		}

		bool should_close()
		{
			return glfwWindowShouldClose(window) || glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS;
		}

		template <typename MolSys>
		void draw(const MolSys& molsys)
		{
			using std::remainder;
			atomPosType.resize(4*molsys.n);

			for (unsigned i = 0; i < molsys.n; ++i)
			{
				atomPosType[i*4+0] = remainder(molsys.x[i][0], molsys.side);
				atomPosType[i*4+1] = remainder(molsys.x[i][1], molsys.side);
				atomPosType[i*4+2] = remainder(molsys.x[i][2], molsys.side);
				atomPosType[i*4+3] = physics::atom_number[int(molsys.id[i])];
			}

			glfwPollEvents();

			glBindFramebuffer(GL_FRAMEBUFFER, fb);

			glDrawBuffer(GL_COLOR_ATTACHMENT0);

			glClearColor(0, 0, 0, 1);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			glViewport(0, 0, w, h);

			glUseProgram(progID);

			//glEnable(GL_MULTISAMPLE);
			glEnable(GL_DEPTH_TEST);
			glDepthFunc(GL_GEQUAL);
			glClearDepth(0);

			glClipControl(GL_LOWER_LEFT, GL_ZERO_TO_ONE);

			glEnable(GL_CULL_FACE);
			glFrontFace(GL_CCW);
			//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

			glDisable(GL_BLEND);

			updateControls(window, w, h, (float)dt);

			glm::dmat4 Model(1);
			glm::dmat4 dMV = controls::View * Model;
			glm::mat4 MV = dMV;
			glm::mat4 fProj = controls::Proj;

			glUniformMatrix4fv(mvID, 1, GL_FALSE, &MV[0][0]);
			glUniformMatrix4fv(projID, 1, GL_FALSE, &fProj[0][0]);

			//glm::mat3 NormalMat(transpose(inverse(dMV)));

			//glUniformMatrix3fv(normalMatID, 1, GL_FALSE, &NormalMat[0][0]);

			glm::dvec4 lightPosWorld = glm::dvec4(
				400*std::sin(2*lastTime/100),
				400*std::sin(std::sqrt(2.)*lastTime/100),
				400*std::cos(2*lastTime/100),
				1
			);
			glm::vec3 lightPos = controls::View * lightPosWorld;

			glUniform3fv(lightPosID, 1, &lightPos[0]);
			glUniform1f(gammaID, controls::gamma);

			glEnableVertexAttribArray(0);
			glBindBuffer(GL_ARRAY_BUFFER, vb_atom);
			glVertexAttribPointer(
				0,			// location id in vertex shader
				3,			// size
				GL_FLOAT,	// type
				GL_FALSE,	// normalized?
				0,			// stride
				nullptr		// array buffer offset
			);

			if (molsys.n > 0)
			{
				glEnableVertexAttribArray(10);
				glBindBuffer(GL_ARRAY_BUFFER, pb_inst);
				glBufferData(GL_ARRAY_BUFFER, physics::molecular_system<float>::max_atoms*4*sizeof(float), nullptr, GL_STREAM_DRAW);
				glBufferSubData(GL_ARRAY_BUFFER, 0, molsys.n*4*sizeof(float), atomPosType.data());
				glVertexAttribPointer(
					10,			// location id in vertex shader
					4,			// size
					GL_FLOAT,	// type
					GL_FALSE,	// normalized?
					0,			// stride
					nullptr		// array buffer offset
				);

				// glVertexAttribDivisor(m, n)
				// n = 0 => reuse the entire vertex buffer m for each object (mesh)
				// otherwise, use a different vertex element m every n object (instancing buffer, n=1)
				glVertexAttribDivisor(0, 0);
				glVertexAttribDivisor(10, 1);

				glDrawArraysInstanced(
					GL_TRIANGLE_STRIP,	// type of primitive
					0,					// offset
					mesh_size/3,		// number of vertices
					molsys.n
				);
			}

			glDisableVertexAttribArray(0);
			glDisableVertexAttribArray(10);

			// POST-PROCESSING

			// FXAA (fast approximate anti-aliasing)

			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

			glBindFramebuffer(GL_FRAMEBUFFER, 0);
			glClearColor(1, 1, 1, 1);
			glClear(GL_COLOR_BUFFER_BIT);

			glDisable(GL_BLEND);
			glDisable(GL_DEPTH_TEST);

			glUseProgram(progFXAAID);

			glUniform1i(sceneID, 0);
			glUniform1i(blueNoiseID, 1);

			glUniform3f(rgbMaxID, (1 << (vidmode -> redBits - 2)) - 1,
								  (1 << (vidmode -> greenBits - 2)) - 1,
								  (1 << (vidmode -> blueBits - 2)) - 1); // (1 << (vidmode -> redBits - 2)) - 1
			glUniform1ui(frameID, frame);

			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, texScene);
			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, texBlueNoise);

			renderQuad();

			// TEXT (render on same framebuffer)

			glEnable(GL_BLEND);
			glBlendEquation(GL_FUNC_ADD);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			glDisable(GL_DEPTH_TEST);

			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

			font -> Begin();

			font -> Color(0, 1, 0, 1);
			font -> Draw(
				std::string("FPS: ") + std::to_string(FPS) + '\n' + std::to_string(controls::pos.x) + ' '
																  + std::to_string(controls::pos.y) + ' '
																  + std::to_string(controls::pos.z),
				12*xscale, 24*yscale, 0.5f*xscale
			);
			font -> Draw(
				std::string("FoV: ") + std::to_string(controls::FoV) + "\ngamma: " + std::to_string(controls::gamma),
				w-160*xscale, 36*yscale, 0.5f*xscale
			);
			font -> Draw(
				std::string("T = ") + std::to_string(molsys.temperature()) +
				std::string(" K\nP = ") + std::to_string(molsys.pressure_fixedT()/physics::atm<double>) +
				std::string(" atm\nV = ") + std::to_string(molsys.volume()/1'000) +
				std::string(" nm^3\nE = ") + std::to_string(molsys.total_energy()) +
				std::string(" kcal/mol\nN = ") + std::to_string(molsys.n) +
				std::string("\nrho = ") + std::to_string(molsys.density()/physics::kg_per_m3<double>) +
				std::string(" kg/m^3"),
				w-160*xscale, h-36*yscale, 0.5f*xscale
			);
			font -> End();
			
			glfwSwapBuffers(window);

			GLenum err;
			if ((err = glGetError()) != GL_NO_ERROR)
			{
				std::cerr << "Error with code: " << err << std::endl;
				throw;
			}

			double nextTime = glfwGetTime();
			dt = nextTime - lastTime;
			if (FPS == 0)
				FPS = 1 / dt;
			else
				FPS = (FPS * 9 + 1 / dt) / 10;
			dt = std::min(dt, 1.);
			++frame;
			lastTime = nextTime;
		}

	private:

		GLFWwindow* window;
		const GLFWvidmode *vidmode;
		Font *font;
		GLuint progID, progFXAAID; // programs
		GLuint mvID, projID, normalMatID, lightPosID, gammaID; // uniforms of progID
		GLuint rgbMaxID, frameID, sceneID, blueNoiseID; // uniforms of progFXAAID
		GLuint fb; // framebuffers
		GLuint texScene, texBlueNoise; // textures
		GLuint va, rb; // vertex arrays and render buffers
		GLuint vb_atom, // atom buffers
			   pb_inst, // (position) instancing buffer
			   vbPost; // post-processing vertex buffer (square over the screen)
				// vertex buffer objects
		GLint fbIntFormat = GL_RGBA16F; // internal format (half-precision floating point)
		GLenum fbFormat = GL_RGBA;
		unsigned int mesh_size;
		float xscale, yscale;
		int w, h;
		double lastTime, dt;
		unsigned int frame;
		int FPS;
		std::vector<float> atomPosType;

		static void DbgMessage(GLenum, GLenum, GLuint, GLenum, GLsizei, const GLchar *message, const void*)
		{
			std::cerr << message << std::endl;
		}

		int Init(const int w, const int h)
		{
			glewExperimental = GL_TRUE;
			GLenum st = glewInit();

			if (st != GLEW_OK)
			{
				std::cerr << "GLEW Error: " << glewGetErrorString(st) << std::endl;
				return -1;
			}

			glEnable(GL_DEBUG_OUTPUT);
			glPatchParameteri(GL_PATCH_VERTICES, 3);

			glDebugMessageCallback(DbgMessage, nullptr);

			std::clog << "GLEW version: " << glewGetString(GLEW_VERSION) << std::endl;

			glGenVertexArrays(1, &va);
			glBindVertexArray(va);

			mesh_size = 4*3;
			float verts[]
			{
				-1, 1, 0,
				-1, -1, 0,
				1, 1, 0,
				1, -1, 0,
			}; // square

			glGenBuffers(1, &vb_atom);
			glBindBuffer(GL_ARRAY_BUFFER, vb_atom);
			glBufferData(GL_ARRAY_BUFFER, mesh_size*sizeof(float), verts, GL_STATIC_DRAW);

			glGenBuffers(1, &pb_inst);
			glBindBuffer(GL_ARRAY_BUFFER, pb_inst);
			glBufferData(GL_ARRAY_BUFFER, physics::molecular_system<float>::max_atoms*4*sizeof(float), nullptr, GL_STREAM_DRAW);

			GLfloat square[12]
			{
				0, 1,
				0, 0,
				1, 0,
				0, 1,
				1, 0,
				1, 1,
			};
			glGenBuffers(1, &vbPost);
			glBindBuffer(GL_ARRAY_BUFFER, vbPost);
			glBufferData(GL_ARRAY_BUFFER, 12*sizeof(GLfloat), square, GL_STATIC_DRAW);

			glGenFramebuffers(1, &fb);
			glBindFramebuffer(GL_FRAMEBUFFER, fb);

			glGenTextures(1, &texScene);
			glBindTexture(GL_TEXTURE_2D, texScene);
			glTexImage2D(GL_TEXTURE_2D, 0, fbIntFormat, w, h, 0, fbFormat, GL_FLOAT, nullptr);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR); // without bilinear interpolation FXAA does not work properly!
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
			glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texScene, 0);
			glBindTexture(GL_TEXTURE_2D, 0);

			glGenRenderbuffers(1, &rb);
			glBindRenderbuffer(GL_RENDERBUFFER, rb);
			glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT32F, w, h);
			glBindRenderbuffer(GL_RENDERBUFFER, 0);

			glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rb);

			if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
				std::cerr << "Error: Framebuffer is not complete!" << std::endl;

			glBindFramebuffer(GL_FRAMEBUFFER, 0);

			glGenTextures(1, &texBlueNoise);
			glBindTexture(GL_TEXTURE_2D, texBlueNoise);
			int bnw, bnh, bnc;
			unsigned char *blueNoiseData = stbi_load("LDR_RGB1_0.png", &bnw, &bnh, &bnc, 3);
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, bnw, bnh, 0, GL_RGB, GL_UNSIGNED_BYTE, blueNoiseData);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

			glBindTexture(GL_TEXTURE_2D, 0);

			std::free(blueNoiseData); // was allocated with malloc() and not new[]

			progFXAAID = loadShader("shaders/vertexpost.vsh", "shaders/fragmentfxaa.fsh");

			progID = loadShader("shaders/vertex.vsh", "shaders/fragment.fsh");

			if (!progID || !progFXAAID)
				return -1;

			font = new Font("LiberationSans-Regular.ttf", 24, w, h);

			if (!font -> good())
				return -1;

			mvID = glGetUniformLocation(progID, "MV");
			projID = glGetUniformLocation(progID, "Proj");
			normalMatID = glGetUniformLocation(progID, "NormalMat");
			lightPosID = glGetUniformLocation(progID, "lightPos");
			gammaID = glGetUniformLocation(progID, "gamma");

			rgbMaxID = glGetUniformLocation(progFXAAID, "rgbMax");
			frameID = glGetUniformLocation(progFXAAID, "frame");
			sceneID = glGetUniformLocation(progFXAAID, "screenTex");
			blueNoiseID = glGetUniformLocation(progFXAAID, "blueNoiseTex");

			return 0;
		}

		static void ErrorCallback(int, const char* description)
		{
			std::cout << description << std::endl;
		}

		void renderQuad()
		{
			glEnableVertexAttribArray(0);
			glBindBuffer(GL_ARRAY_BUFFER, vbPost);
			glVertexAttribPointer(
				0,			// location id in vertex shader
				2,			// size
				GL_FLOAT,	// type
				GL_FALSE,	// normalized?
				0,			// stride
				nullptr		// array buffer offset
			);

			glDrawArrays(GL_TRIANGLES, 0, 6);

			glDisableVertexAttribArray(0);
		}

};

#endif // GRAPHICS_H



















