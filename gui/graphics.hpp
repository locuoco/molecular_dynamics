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

#ifndef GUI_GRAPHICS_H
#define GUI_GRAPHICS_H

#include <iostream> // cerr, clog, endl
#include <vector>
#include <cmath> // remainder, sin, cos
#include <numbers> // numbers::sqrt2
#include <string>
#include <exception> // runtime_error

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "shader.hpp" // load_shader
#include "controls.hpp"
#include "../physics/molecule.hpp"

#define USE_PHYSICS
#include "font.hpp"

// define EXPORT for MinGW, MSVC and GCC compilers
// other compilers are not supported for this
#ifdef __MINGW32__
#define EXPORT __declspec(dllexport)
#elif defined(_MSC_VER)
#define EXPORT _declspec(dllexport)
#elif defined(__GNUC__)
#define EXPORT __attribute__((visibility("default")))
#else
#define EXPORT // nothing
#endif

// enable optimus, which chooses a dedicated GPU if present.
// Should work with MinGW, MSVC and GCC compilers
extern "C"
{
    EXPORT unsigned long NvOptimusEnablement = 1;
    EXPORT int AmdPowerXpressRequestHighPerformance = 1;
}

#define CHECK_GL_ERROR() checkGlError(__FILE__, __LINE__)

inline void checkGlError(const char *file, int line)
// check if there are any errors, throwing a `std::runtime_error`
// with a message containing the file and the line where it occurred.
// `file` is the path and filename where the check is done.
// `line` is the line number where the check is done.
// use the macro `CHECK_GL_ERROR()` to automatically generate the
// `file` and `line` parameters.
{
	if (GLenum err = glGetError(); err != GL_NO_ERROR)
	{
		std::cerr << "Error with code: " << err << '\n';
		throw std::runtime_error(std::string("Error message: ") + reinterpret_cast<const char*>(gluErrorString(err))
			+ ' ' + file + ' ' + std::to_string(line));
	}
}

class graphics
{
	public:

		graphics() : dt(0), frame(0), fps(0)
		// (default) constructor
		// it throws `std::runtime_error` if GLFW could not be initialized, if the window could
		// not be created or if init_resources() throws
		{
			glfwSetErrorCallback(&ErrorCallback);
			// OpenGL graphic interface initialization (GLFW + GLEW)
			if (!glfwInit())
				throw std::runtime_error("Error: Cannot initialize GLFW.");
			// This class requires OpenGL 3.3 (March 11th 2010)
			glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
			glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
			// Requesting core profile, i.e. modern OpenGL API rather than
			// the compatibility (deprecated) one
			glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
			// Forward compatibility is needed for OpenGL 3.2+ support in MacOS (see GLFW docs)
			glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
			// Make the window not resizable
			glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);

			// vidmode is used to know the bit-depth for each color channel of the monitor being used
			vidmode = glfwGetVideoMode(glfwGetPrimaryMonitor());

			// full-HD and better monitors may have more pixels per inch squared,
			// so the window needs to be rescaled accordingly
			glfwGetMonitorContentScale(glfwGetPrimaryMonitor(), &xscale, &yscale);
			w = 1280*xscale, h = 720*yscale;

			// create the window. `nullptr` for the 4th parameter for windowed mode.
			// Otherwise, glfwGetPrimaryMonitor(), which returns the pointer to the window
			// object associated to the monitor, will enable full-screen mode.
			window = glfwCreateWindow(w, h, "Molecular Dynamics", nullptr, nullptr);

			if (!window)
				throw std::runtime_error("GLFW Error: Cannot create window.");

			// glfwMakeContextCurrent can be used to switch between contexts,
			// for just one window this function must be called only once at the beginning
			glfwMakeContextCurrent(window);
			// VSync is disabled (set glfwSwapInterval to 1 to enable)
			glfwSwapInterval(0);

			std::clog << "GPU: " << glGetString(GL_RENDERER) << std::endl;

			camera_controls = new controls(window, w, h);

			init_resources();

			last_time = glfwGetTime();
		}

		graphics(const graphics&) = delete;
		graphics& operator=(const graphics&) = delete;

		~graphics()
		// destructor
		{
			delete text;
			delete camera_controls;

			glDeleteProgram(prog_id);
			glDeleteProgram(prog_post_id);

			// delete all buffers (free GPU memory)
			glDeleteFramebuffers(1, &fb);
			glDeleteTextures(1, &tex_scene);
			glDeleteTextures(1, &tex_blue_noise);
			glDeleteRenderbuffers(1, &rb);
			glDeleteBuffers(1, &vb_post);
			glDeleteBuffers(1, &pb_inst);
			glDeleteBuffers(1, &vb_atom);
			glDeleteVertexArrays(1, &va_post);
			glDeleteVertexArrays(1, &va_atom);

			glfwTerminate();
		}

		bool should_close()
		// return true if the quit button is selected or if `ESC` has been pressed
		{
			return glfwWindowShouldClose(window) || glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS;
		}

		template <typename MolSys, typename CharT = char>
		void draw(MolSys& molsys, const std::basic_string<CharT>& custom_text = "")
		// draw all the atoms contained in `molsys` system.
		// `custom_text` is an (optional) text to print on screen.
		// `Molsys` must be a `molecular_system` type.
		// Throw a `std::runtime_error` if `CHECK_GL_ERROR` finds an error
		{
			using std::remainder;
			molsys.fetch();
			atom_attribs.resize(4*molsys.n);
			for (unsigned i = 0; i < molsys.n; ++i)
			{
				// the first three components are the position of the atom
				// x_i in [-side_i/2, side_i/2]^3 where `side_i` is one side of the simulation box
				atom_attribs[i*4+0] = remainder(molsys.x[i][0], molsys.side[0]);
				atom_attribs[i*4+1] = remainder(molsys.x[i][1], molsys.side[1]);
				atom_attribs[i*4+2] = remainder(molsys.x[i][2], molsys.side[2]);
				// the fourth component is the identity of the atom (needed
				// to choose color and size of the sphere)
				atom_attribs[i*4+3] = physics::atom_number[int(molsys.id[i])];
			}

			// update controls
			camera_controls->update(dt);

			// drawing subroutines
			draw_background_and_spheres(molsys.n);
			draw_post_processing();
			draw_text(molsys, custom_text);

			// swap the default framebuffer with the window color buffer
			// (because GLFW uses a double buffer by default for performance purposes)
			glfwSwapBuffers(window);

			CHECK_GL_ERROR();

			estimate_fps();
		}

	controls *camera_controls = nullptr;

	private:

		GLFWwindow* window = nullptr;
		const GLFWvidmode *vidmode = nullptr;
		font *text = nullptr;
		GLuint prog_id, prog_post_id; // programs (3D sphere rendering and post-processing)
		GLuint mv_id, proj_id, normal_mat_id, light_pos_id, gamma_id; // uniforms of prog_id
		GLuint rgb_max_id, frame_id, scene_id, blue_noise_id; // uniforms of prog_post_id
		GLuint fb, rb; // framebuffers and render buffers
		GLuint tex_scene, tex_blue_noise; // textures
		GLuint va_atom, va_post; // vertex arrays
		GLuint vb_atom, // atom buffer object
			   pb_inst, // (position) instancing buffer bject
			   vb_post; // post-processing vertex buffer object (square over the screen)
		unsigned int mesh_size;
		float xscale, yscale;
		int w, h;
		double last_time, dt;
		unsigned int frame;
		int fps;
		std::vector<float> atom_attribs;

		static void DbgMessage(GLenum, GLenum, GLuint, GLenum, GLsizei, const GLchar *message, const void*)
		// callback used for OpenGL debugging.
		// `message` is the message received from OpenGL debugger.
		// All other arguments are ignored.
		{
			std::clog << message << std::endl;
		}

		static void ErrorCallback(int, const char* description)
		// callback used for GLFW errors.
		// `description` is the message received from GLFW debugger.
		{
			std::cerr << description << std::endl;
		}

		void init_resources()
		// initialize and load all initial resources on GPU
		// throws `std::runtime_error` if GLEW could not be initialized, if shaders or font
		// could not be loaded or if `CHECK_GL_ERROR` finds an error
		{
			// setting glewExperimental is needed for GLEW 1.13 or earlier
			glewExperimental = GL_TRUE;
			// initialize GLEW
			if (GLuint st = glewInit(); st != GLEW_OK)
				throw std::runtime_error(std::string("GLEW Error: ") + reinterpret_cast<const char*>(glewGetErrorString(st)));

			// enable OpenGL debug message callback
			// requires OpenGL 4.3 to work properly
			glEnable(GL_DEBUG_OUTPUT);
			glDebugMessageCallback(DbgMessage, nullptr);

			CHECK_GL_ERROR();

			std::clog << "GLEW version: " << glewGetString(GLEW_VERSION) << std::endl;

			// initialize buffers
			init_atom_vertex_buffers();
			init_post_vertex_buffers();
			init_framebuffer();
			init_blue_noise();

			// load vertex and fragment shaders for post-processing (FXAA+dithering) into the `prog_post_id` program
			prog_post_id = load_shader("gui/shaders/vertexpost.vsh", "gui/shaders/fragmentpost.fsh");

			// load vertex and fragment shaders for colored-illuminated spheres into the `prog_id` program
			prog_id = load_shader("gui/shaders/vertex.vsh", "gui/shaders/fragment.fsh");

			if (!prog_id || !prog_post_id)
				throw std::runtime_error("Error in loading shader programs.");

			text = new font("gui/LiberationSans-Regular.ttf", 24, w, h);

			if (!text->good())
				throw std::runtime_error("Error in loading font.");

			CHECK_GL_ERROR();

			// Get the addresses of the following uniforms defined inside vertex/fragment shaders
			mv_id = glGetUniformLocation(prog_id, "MV");
			proj_id = glGetUniformLocation(prog_id, "Proj");
			normal_mat_id = glGetUniformLocation(prog_id, "NormalMat");
			light_pos_id = glGetUniformLocation(prog_id, "lightPos");
			gamma_id = glGetUniformLocation(prog_id, "gamma");

			rgb_max_id = glGetUniformLocation(prog_post_id, "rgbMax");
			frame_id = glGetUniformLocation(prog_post_id, "frame");
			scene_id = glGetUniformLocation(prog_post_id, "screenTex");
			blue_noise_id = glGetUniformLocation(prog_post_id, "blueNoiseTex");

			CHECK_GL_ERROR();
		}

		void init_atom_vertex_buffers()
		// initialize vertex buffers to draw atoms.
		// Throw a `std::runtime_error` if `CHECK_GL_ERROR` finds an error
		{
			// create a vertex array object and bind to it
			// all vertex buffers are created and used within a vertex array
			glGenVertexArrays(1, &va_atom);
			glBindVertexArray(va_atom);

			// vertices of a 2x2 square
			// only 4 vertices needed if the square is drawn with a triangle strip
			// (instead of a triangle list)
			mesh_size = 4*2;
			float verts[]
			{
				-1, 1,
				-1, -1,
				1, 1,
				1, -1,
			};

			// these rows specify that the inputs of vertex shader have
			// the following ids (0 and 10). 0 will be used for the square vertices
			// and 10 for the instance buffer (containing positions and identities
			// of all atoms).
			glEnableVertexAttribArray(0);
			glEnableVertexAttribArray(10);

			// generate vertex buffer for the single atom (it is a square, the sphere
			// will be rendered inside it as an impostor), sending the data to the GPU
			// `GL_STATIC_DRAW` promises that the buffer won't be written again
			glGenBuffers(1, &vb_atom);
			glBindBuffer(GL_ARRAY_BUFFER, vb_atom);
			glBufferData(GL_ARRAY_BUFFER, mesh_size*sizeof(float), verts, GL_STATIC_DRAW);
			// specify `vb_atom` buffer format and attributes
			glVertexAttribPointer(
				0,        // location id in vertex shader
				2,        // number of coordinates (horizontal and vertical)
				GL_FLOAT, // type
				GL_FALSE, // normalized? (used for fixed point types)
				0,        // stride (0 assumes no padding)
				nullptr   // array buffer offset
			);

			// generate instancing buffer. It will store the positions of all atoms, and
			// a square will be drawn at each position. `nullptr` as the third parameter
			// of `glBufferData` means that the buffer is only allocated for now.
			// `GL_DYNAMIC_DRAW` promises that the buffer will be written to frequently
			// (every frame indeed)
			glGenBuffers(1, &pb_inst);
			glBindBuffer(GL_ARRAY_BUFFER, pb_inst);
			glBufferData(GL_ARRAY_BUFFER, physics::max_atoms*4*sizeof(float), nullptr, GL_DYNAMIC_DRAW);
			// specify `pb_inst` buffer format and attributes
			glVertexAttribPointer(
				10,       // location id in vertex shader
				4,        // number of coordinates (3 spatial, 1 for atom id)
				GL_FLOAT, // type
				GL_FALSE, // normalized? (used for fixed point types)
				0,        // stride (0 assumes no padding)
				nullptr   // array buffer offset
			);

			// glVertexAttribDivisor(id, n) is used for instanced draw and will cause the following to happen:
			// if n = 0 => reuse the entire vertex buffer for the specified id for each instance (mesh)
			// otherwise, use a different element for the specified id every `n` instances
			// for n = 1 for an instancing buffer (like `pb_inst`, content copied from `atom_attribs`)
			glVertexAttribDivisor(0, 0);
			glVertexAttribDivisor(10, 1);

			// unbind vertex buffer
			glBindBuffer(GL_ARRAY_BUFFER, 0);
			// unbind vertex array
			glBindVertexArray(0);

			CHECK_GL_ERROR();
		}

		void init_post_vertex_buffers()
		// initialize vertex buffer for post-processing.
		// Throw a `std::runtime_error` if `CHECK_GL_ERROR` finds an error
		{
			// create a vertex array object and bind to it
			// all vertex buffers are created and used within a vertex array
			glGenVertexArrays(1, &va_post);
			glBindVertexArray(va_post);

			// create another square for the screen texture used for post-processing
			// (a triangle list will be used this time, for the sake of variation)
			float square[12]
			{
				0, 1,
				0, 0,
				1, 0,
				0, 1,
				1, 0,
				1, 1,
			};
			// similarly as what has been done in `init_atom_vertex_buffers`
			glGenBuffers(1, &vb_post);
			glBindBuffer(GL_ARRAY_BUFFER, vb_post);
			glBufferData(GL_ARRAY_BUFFER, 12*sizeof(GLfloat), square, GL_STATIC_DRAW);
			// enable vertex array
			glEnableVertexAttribArray(0);
			// specify `vb_post` buffer format and attributes
			glVertexAttribPointer(
				0,        // location id in vertex shader
				2,        // number of coordinates per vertex
				GL_FLOAT, // type
				GL_FALSE, // normalized? (for fixed point types)
				0,        // stride (0 assumes no padding)
				nullptr   // array buffer offset
			);

			// unbind vertex buffer
			glBindBuffer(GL_ARRAY_BUFFER, 0);
			// unbind vertex array
			glBindVertexArray(0);

			CHECK_GL_ERROR();
		}

		void init_framebuffer()
		// create `fb` framebuffer, `rb` renderbuffer (containing a 32-bit float depth buffer)
		// and `tex_scene` texture (onto which the sphere will be drawn to before post-processing)
		// then attach `tex_scene` texture and `rb` renderbuffer to `fb`.
		// Throw a `std::runtime_error` if `CHECK_GL_ERROR` finds an error
		{
			// generate a framebuffer, onto which we can attach a texture and a render buffer
			// and do nasty things (i.e. rendering the scene)...
			glGenFramebuffers(1, &fb);
			// ... and bind to it
			glBindFramebuffer(GL_FRAMEBUFFER, fb);

			// generate a texture and bind to it (it will stores the pixels/colors of the scene)
			glGenTextures(1, &tex_scene);
			glBindTexture(GL_TEXTURE_2D, tex_scene);

			// specify the texture format (and allocates memory on GPU)
			// `GL_RGBA16F` is the texture internal format (one half-precision floating point per channel)
			// its dimensions will be the same as the window (w x h)
			// `GL_RGBA` is the channel format: (red, green, blue, alpha)
			// the alpha channel is unused. `GL_RGBA16F` is chosen since it is more supported than `GL_RGB16F`
			// and requires less space than `GL_RGB32F`.
			// `GL_FLOAT` means that we want a floating point texture rather than integer one.
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, w, h, 0, GL_RGBA, GL_FLOAT, nullptr);
			// min filter set to nearest, we don't need it
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
			// mag filter set to (bi)linear. It means that when sampling a color in positions
			// which do not correspond to a pixel, a bilinear interpolation is performed between neighboring
			// pixels (without it FXAA does not work properly!)
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
			// "clamp to edge" means that when sampling a color outside the texture, the pixel
			// from the closest position is returned. S and T are the two coordinates of the texture
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
			// attach the texture to the current framebuffer (i.e. `fb`) as a color buffer with index 0
			glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex_scene, 0);
			glBindTexture(GL_TEXTURE_2D, 0);

			// create a w x h renderbuffer and use a 32-bit depth buffer.
			// A depth buffer stores the depth of a pixel. Before drawing a fragment over a pixel,
			// its depth will be checked against the current depth. Only fragments closest to the viewer
			// will be conserved, the others will be either overwritten or discarded.
			glGenRenderbuffers(1, &rb);
			glBindRenderbuffer(GL_RENDERBUFFER, rb);
			glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT32F, w, h);
			glBindRenderbuffer(GL_RENDERBUFFER, 0);

			// attach the renderbuffer to the `fb` framebuffer
			glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rb);

			if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
				std::cerr << "Error: Framebuffer is not complete!" << std::endl;

			// bind back to the default framebuffer (i.e., the window)
			glBindFramebuffer(GL_FRAMEBUFFER, 0);

			CHECK_GL_ERROR();
		}

		void init_blue_noise()
		// load blue noise image from file and store it inside a texture on GPU
		// Throw a `std::runtime_error` if `CHECK_GL_ERROR` finds an error
		{
			// create a new texture
			glGenTextures(1, &tex_blue_noise);
			glBindTexture(GL_TEXTURE_2D, tex_blue_noise);
			// load blue noise texture from file (3 is the number of expected channels)
			int bnw, bnh, bnc;
			unsigned char *blueNoiseData = stbi_load("gui/LDR_RGB1_0.png", &bnw, &bnh, &bnc, 3);
			// GL_RGB8 is the internal format, it is the standard 24-bit format for LDR (low dynamic-range) images.
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, bnw, bnh, 0, GL_RGB, GL_UNSIGNED_BYTE, blueNoiseData);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
			// `GL_REPEAT` intuitively means that the colors outside the [0, 1] range are sampled as if
			// the texture is periodic
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

			glBindTexture(GL_TEXTURE_2D, 0);

			stbi_image_free(blueNoiseData);

			CHECK_GL_ERROR();
		}

		void draw_background_and_spheres(std::size_t n)
		// set background to black and draw the spheres on the framebuffer called `fb`
		// (in particular, colors are written inside the framebuffer color attachment
		// represented by the texture `tex_scene`)
		// `n` is the number of particles to draw (given by `molsys.n`).
		{
			// convert view and projection matrices into single-precision floating-point matrices
			physics::mat4f fview = camera_controls->view;
			physics::mat4f fproj = camera_controls->proj;

			// a light is rotating with time (`last_time` is a measure of real time, not the simulation time)
			physics::vec4d lightPosWorld(
				0,
				0,
				400,
				1
			);
			// moving into view coordinates
			physics::vec4f lightPos = camera_controls->view % lightPosWorld;

			// bind to the `fb` framebuffer, onto which we will draw the spheres
			glBindFramebuffer(GL_FRAMEBUFFER, fb);
			// set draw buffer as the first (and only) color attachment of the framebuffer
			glDrawBuffer(GL_COLOR_ATTACHMENT0);

			// enable depth test, which checks the depths of fragments when drawing
			glEnable(GL_DEPTH_TEST);
			// since we are using a reversed depth buffer, the depth function is
			// GL_GEQUAL (the default is GL_LESS)
			glDepthFunc(GL_GEQUAL);

			// depth range set to [-1, 1]
			// the depth is always calculated in the range [-1, 1] and then mapped
			// to the range specified by glDepthRange. Choosing [-1, 1] as the range
			// means no mapping is performed.
			glDepthRange(-1, 1);
			// clear-depth set to 0 (default is 1)
			glClearDepth(0);

			// clear-color set to black
			glClearColor(0, 0, 0, 1);
			// clear color and depth buffers
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			// set the viewport as the whole window/screen
			// one might choose to draw only on a part of the window by setting the parameters accordingly
			glViewport(0, 0, w, h);

			// use the shader program that draws spheres (loaded in init method)
			glUseProgram(prog_id);

			// cull faces whose vertices are in clockwise order with respect to the viewer
			// (actually has no real effect in our case since all spheres are rendered inside squares
			// that are always oriented towards the viewer, so there are no back faces)
			glEnable(GL_CULL_FACE);
			glFrontFace(GL_CCW); // front face is counter-clockwise
			// disabling blending
			glDisable(GL_BLEND);

			if (n == 0)
				return; // no particles to draw

			// set the uniforms (variables defined inside the shaders)
			// GL_TRUE means that the matrices will be transposed, since GLSL
			// uses a column-major syntax
			glUniformMatrix4fv(mv_id, 1, GL_TRUE, &fview(0, 0)); // View matrix
			glUniformMatrix4fv(proj_id, 1, GL_TRUE, &fproj(0, 0)); // Projection matrix

			// set light position
			glUniform3fv(light_pos_id, 1, &lightPos[0]);
			// set gamma correction exponent
			glUniform1f(gamma_id, camera_controls->gamma);

			// use the atom vertex array
			glBindVertexArray(va_atom);

			glBindBuffer(GL_ARRAY_BUFFER, pb_inst);
			// write on part of the `pb_inst` buffer allocated in `init` method
			glBufferSubData(GL_ARRAY_BUFFER, 0, atom_attribs.size()*sizeof(float), atom_attribs.data());

			glDrawArraysInstanced(
				GL_TRIANGLE_STRIP,    // type of primitive
				0,                    // offset
				mesh_size/2,          // number of vertices per instance
				atom_attribs.size()/4 // number of instances
			);
		}

		void draw_post_processing()
		// perform FXAA (fast approximate anti-aliasing) and blue noise dithering
		{
			// bind to default (0) framebuffer (i.e. the screen)
			glBindFramebuffer(GL_FRAMEBUFFER, 0);

			// disable blending
			glDisable(GL_BLEND);
			// disable depth test
			glDisable(GL_DEPTH_TEST);
			// use the shader program for post-processing (loaded in init method) 
			glUseProgram(prog_post_id);

			// set the texture index for the scene (tex_scene)
			glUniform1i(scene_id, 0);
			// set the texture index for the blue noise (tex_blue_noise)
			glUniform1i(blue_noise_id, 1);
			// set the number of quantization levels for each channel
			// used for dithering
			glUniform3f(rgb_max_id, (1 << (vidmode->redBits - 2)) - 1,
			                        (1 << (vidmode->greenBits - 2)) - 1,
			                        (1 << (vidmode->blueBits - 2)) - 1);
			// set frame number
			glUniform1ui(frame_id, frame);

			// bind scene texture to index 0
			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, tex_scene);
			// bind blue noise texture to index 1
			glActiveTexture(GL_TEXTURE1);
			glBindTexture(GL_TEXTURE_2D, tex_blue_noise);

			// use the post-processing vertex array
			glBindVertexArray(va_post);

			// 6 vertices are used for drawing a square (composed of two triangles)
			glDrawArrays(GL_TRIANGLES, 0, 6);
		}

		template <typename Molsys, typename CharT>
		void draw_text(Molsys& molsys, const std::basic_string<CharT>& custom_text)
		// draw text on screen
		{
			// bind to default (0) framebuffer (i.e. the screen)
			glBindFramebuffer(GL_FRAMEBUFFER, 0);

			// enable alpha blending
			// when fragments overlaps, one can choose how the resulting color
			// is calculated. The destination color is the one already in the buffer while
			// the source color is the color that must be added/blended.
			glEnable(GL_BLEND);

			// blended colors are added in this way:
			// c_r = c_s * f_s + c_d * f_d
			// where c_r is the resulting color, c_s is the source color, c_d is the destination
			// color, while f_s is the source blend factor and f_d is the destination blend factor.
			glBlendEquation(GL_FUNC_ADD);

			// source blend factor (f_s) is given by the source alpha channel
			// destination blend factor (f_d) is given by 1-alpha, where alpha is the source alpha channel
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

			// disable depth test
			glDisable(GL_DEPTH_TEST);

			text->begin();
			text->color(0, 1, 0, 1); // opaque green color
			text->draw(
				std::string("FPS: ") + std::to_string(fps) + '\n' + std::to_string(camera_controls->pos[0]) + ' '
				                                                  + std::to_string(camera_controls->pos[1]) + ' '
				                                                  + std::to_string(camera_controls->pos[2]),
				12*xscale, 24*yscale, 0.5f*xscale
			);
			text->draw(
				std::string("FoV: ") + std::to_string(camera_controls->fov) + "\ngamma: " + std::to_string(camera_controls->gamma),
				w-160*xscale, 36*yscale, 0.5f*xscale
			);
			text->draw(
				std::wstring(L"T = ") + std::to_wstring(molsys.temperature()) +
				L" K\nP = " + std::to_wstring(molsys.pressure()/physics::atm<>) +
				L" atm\nV = " + std::to_wstring(molsys.volume()/1'000) +
				L" nm^3\nE = " + std::to_wstring(molsys.internal_energy()) +
				L" kcal/mol\nN = " + std::to_wstring(molsys.n) +
				L"\nρ = " + std::to_wstring(molsys.density()/physics::kg_per_m3<>) +
				L" kg/m^3",
				w-160*xscale, h-36*yscale, 0.5f*xscale
			);
			text->draw(custom_text, 12*xscale, h-36*yscale, 0.5f*xscale);
			text->end();
		}

		void estimate_fps()
		// estimate `fps` (frames per second) and calculate `dt` (real-time step)
		{
			double next_time = glfwGetTime();
			dt = next_time - last_time;
			if (fps == 0)
				fps = 1 / dt;
			else
				fps = (fps * 9 + 1 / dt) / 10; // moving average
			dt = std::min(dt, 1.);
			++frame;
			last_time = next_time;
		}
};

#endif // GUI_GRAPHICS_H



















