//  Text/Font manager, based on Joey de Vries' code
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

#ifndef GUI_FONT_H
#define GUI_FONT_H

#ifdef ENABLE_DBG
#define FONT_DBG_GL(x) GLenum err_ = glGetError(); if (err_) std::cerr << err_ << " in " << x << std::endl
#else
#define FONT_DBG_GL(x)
#endif

#include <iostream> // cerr, endl
#include <string>
#include <map>
#include <utility> // pair

#include <GL/glew.h>
#include <ft2build.h> // Compile with -lfreetype
#include FT_FREETYPE_H

#include "shader.hpp" // load_shader

#ifdef USE_PHYSICS
// if USE_PHYSICS is defined, then `physics/tensor.hpp` must be available
#include "../physics/tensor.hpp"

#else
// if USE_PHYSICS is undefined, then OpenGL Mathematics library must be available
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>
#include <glm/mat4x4.hpp>
#include <glm/trigonometric.hpp>
#include <glm/gtc/matrix_transform.hpp>

#endif

class font
{
#ifdef USE_PHYSICS
		using vec2i = physics::vec2i;
		using vec3f = physics::vec3f;
		using vec4f = physics::vec4f;
		using mat4f = physics::mat4f;
#else
		using vec2i = glm::ivec2;
		using vec3f = glm::vec3;
		using vec4f = glm::vec4;
		using mat4f = glm::mat4;
#endif

		struct character
		{
			GLuint id;     // texture id
			vec2i sz;      // size of character
			vec2i bearing; // character bearing
			GLint advance; // character advance width
		};

		std::map<wchar_t, character> characters;
		mutable vec4f current_color;

		unsigned int w, h; // window/screen width and height
		static constexpr unsigned int tab = 64, margin = 16; // tabulation and margin in pixels
		GLuint va, vb; // vertex array and vertex buffer IDs
		GLuint transf_id, color_id; // Uniform IDs
		GLuint prog; // GPU program ID
		bool init; // flag which says if the `font` object has been constructed correctly
		
		void load_char(FT_Face ff, FT_ULong c)
		// load a character texture
		// `ff` is the face object from which we will choose the character
		// `c` is the (wide) character to load
		// The character metrics and texture id will be inserted inside the `characters` map
		{
			if (FT_Load_Char(ff, c, FT_LOAD_RENDER))
			{
				std::cerr << "Cannot load character " << c << std::endl;
				return;
			}

			FT_GlyphSlot gs = ff->glyph;

			GLuint texture;
			glGenTextures(1, &texture);
			glBindTexture(GL_TEXTURE_2D, texture);
			glTexImage2D(
				GL_TEXTURE_2D,
				0,
				GL_RED,
				gs->bitmap.width,
				gs->bitmap.rows,
				0,
				GL_RED,
				GL_UNSIGNED_BYTE,
				gs->bitmap.buffer
			);

			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

			character ch = {
				texture, 
				vec2i(gs->bitmap.width, gs->bitmap.rows),
				vec2i(gs->bitmap_left, gs->bitmap_top),
				(GLint)gs->advance.x
			};
			characters.insert(std::pair<wchar_t, character>(c, ch));
		}

		bool color_equal(float r, float g, float b, float a) const
		// Check if the (r, g, b, a) color is exactly the same as the current one
		{
			return r == current_color[0] && g == current_color[1] && b == current_color[2] && a == current_color[3];
		}

		bool color_equal(const vec4f& c) const
		// Check if the `c` color (as a vector of floats) is exactly the same as the current one
		{
			return c == current_color;
		}

	public:

		font() : init(false) {}
		// default constructor. Do not initialize font.

		font(const std::string& font_file, const unsigned int font_size, const int wscreen, const int hscreen)
		// initializing constructor:
		//	`font_file` is the path of the font to use.
		//	`font_size` is the font size to load.
		//	`wscreen` and `hscreen` are the screen (window) width and height respectively.
			: current_color(1,1,1,1), w(wscreen), h(hscreen), init(false)
		{
			prog = load_shader("gui/shaders/vertexfont.vsh", "gui/shaders/fragmentfont.fsh");

			if (!prog)
			{
				std::cerr << "Error: cannot load font shaders." << std::endl;
				return;
			}

			FT_Library ft;

			if (FT_Init_FreeType(&ft))
			{
				std::cerr << "Error: Cannot initialize FreeType." << std::endl;
				return;
			}

			FT_Face ff;

			FT_Error err = FT_New_Face(ft, font_file.c_str(), 0, &ff);

			if (!ff)
			{
				switch (err)
				{
					case FT_Err_Unknown_File_Format:
						std::cerr << "FreeType error: Font format not supported." << std::endl;
						break;
					default:
						std::cerr << "FreeType error: Font couldn't be read. "
									 "The font file may be corrupted." << std::endl;
				}
				FT_Done_FreeType(ft);
				return;
			}

			FT_Set_Pixel_Sizes(ff, 0, font_size);

			glGenVertexArrays(1, &va);
			glBindVertexArray(va);

			glGenBuffers(1, &vb);
			glBindBuffer(GL_ARRAY_BUFFER, vb);
			glBufferData(GL_ARRAY_BUFFER, 24*sizeof(GLfloat), nullptr, GL_DYNAMIC_DRAW);

			glEnableVertexAttribArray(0);
			glVertexAttribPointer(
				0,
				4,
				GL_FLOAT,
				GL_FALSE,
				0,
				(void*)0
			);

			glBindBuffer(GL_ARRAY_BUFFER, 0);
			glBindVertexArray(0);

			transf_id = glGetUniformLocation(prog, "Pmat");
			color_id = glGetUniformLocation(prog, "col");

			glUseProgram(prog);
#ifdef USE_PHYSICS
			mat4f m_Pmat = physics::orthographic_projection(0., double(w), 0., double(h));
			glUniformMatrix4fv(transf_id, 1, GL_TRUE, &m_Pmat(0, 0));
#else
			mat4f m_Pmat = glm::ortho(0., double(w), 0., double(h));
			glUniformMatrix4fv(transf_id, 1, GL_FALSE, &m_Pmat[0][0]);
#endif

			glPixelStorei(GL_UNPACK_ALIGNMENT, GL_TRUE);

			for (wchar_t c = 0x0000; c < 0x0080; ++c) // ASCII characters
				load_char(ff, c);

			for (wchar_t c = 0x0391; c < 0x03CA; ++c) // Greek characters
				load_char(ff, c);

			FT_Done_Face(ff);
			FT_Done_FreeType(ft);

			FONT_DBG_GL("font::font(const char*, unsigned int, int, int, GLuint)");

			init = true;
		}

		// copying disabled
		font(const font&) = delete;
		font& operator=(const font&) = delete;

		font(font&& f)
		// move constructor
			: characters(std::move(f.characters)), current_color(f.current_color), w(f.w), h(f.h), vb(f.vb),
			  transf_id(f.transf_id), color_id(f.color_id), prog(f.prog), init(f.init)
		{
			f.vb = 0;
			f.transf_id = 0;
			f.color_id = 0;
			f.prog = 0;
			f.init = false;
		}

		font& operator=(font&& f)
		// move-assign operator
		{
			if (this != &f)
			{
				characters = std::move(f.characters);
				current_color = f.current_color;
				w = f.w;
				h = f.h;
				vb = f.vb;
				transf_id = f.transf_id;
				color_id = f.color_id;
				prog = f.prog;
				init = f.init;
				f.vb = 0;
				f.transf_id = 0;
				f.color_id = 0;
				f.prog = 0;
				f.init = false;
			}
			return *this;
		}

		~font()
		// destructor
		{
			glDeleteProgram(prog);
			glDeleteBuffers(1, &vb);
			glDeleteVertexArrays(1, &va);
		}

		bool good() const
		// check if the object has been constructed and initialized successfully
		{
			return init;
		}

		void color(float r, float g, float b, float a) const
		// Set the current color to (r, g, b, a)
		{
			if (!color_equal(r, g, b, a))
			{
				glUniform4f(color_id, r, g, b, a);
				current_color = vec4f(r, g, b, a);
				FONT_DBG_GL("font::color(float, float, float, float)");
			}
		}

		void color(float r, float g, float b) const
		// Set the current color to (r, g, b, 1)
		{
			color(r,g,b,1.f);
		}

		void color(const vec3f& col) const
		// Set the current color to (col, 1) (1 is the alpha channel)
		{
			color(col[0],col[1],col[2],1.f);
		}

		void color(const vec4f& col) const
		// Set the current color to `col`
		{
			if (!color_equal(col))
			{
				glUniform4fv(color_id, 1, &col[0]);
				current_color = col;
				FONT_DBG_GL("font::color(const vec4f&)");
			}
		}

		void begin() const
		// function to call before drawing text. It binds the resources
		// used by the draw calls.
		{
			glBindVertexArray(va);
			glBindBuffer(GL_ARRAY_BUFFER, vb);
			glUseProgram(prog);
			glActiveTexture(GL_TEXTURE0);

			FONT_DBG_GL("font::begin()");
		}

		void end() const
		// function to call after drawing text. It unbinds the resources
		// used by the draw calls.
		{
			glBindTexture(GL_TEXTURE_2D, 0);
			glBindBuffer(GL_ARRAY_BUFFER, 0);
			glBindVertexArray(0);
		}

		template <typename CharT>
		void draw(const std::basic_string<CharT>& text, float x, float y, const float scale = 1) const
		// draw a string of text.
		// `text` is the text to draw. It may contain ASCII characters, newlines, tabulations, spaces and greek letters
		//	(greek letters require the use of a wide string `std::wstring` as the underlying type of `text` instead of
		//	a standard `std::string`).
		// `x` and `y` are the coordinates of the upper-left position of the first character to render.
		//	(in pixels and wrt to the upper-left of the window/screen).
		// `scale` is a size scaling factor of the font to render. While the size of the loaded font is fixed, it can be
		//	rescaled with this argument at each draw call.
		// Template arguments:
		// `CharT` is the type of the string character. It is deduced as `char` for `std::string`s and as `wchar_t` for
		//	`std::wstring`s. It does not need to be explicitly specified.
		{
			if (!text.length())
				return;

			const float xi = x; // initial horizontal position, needed for newlines
			// x is the current character horizontal position

			auto newline = [this, xi, scale, &x, &y]
			// lambda function to create a newline
				{
					GLint t = characters.at('A').sz[1]; // 'A' is the tallest character
					t += characters.at('g').bearing[1]; // 'g' is the character with greatest bearing
					x = xi;
					y -= t * scale;
				};

			for (auto i = text.begin(); i != text.end(); ++i)
				switch (*i)
				{
					case '\n':
						newline();
						break;

					case '\t':
						// search for a fixed position such that it is past the current horizontal position x
						for (unsigned int i = 0; i < w - margin; i += tab)
						{
							float t = xi + i;
							if (t > x)
							{
								x = t;
								goto ok;
							}
						}
						// if it is not found within the margin, make a new line
						newline();
						ok:
						break;

					case ' ':
						if (x < w - margin)
							x += (characters.at(*i).advance >> 6) * scale;
						else
							newline();
						break;

					default:
						character ch = characters.at(*i);

						float xpos = x + ch.bearing[0] * scale;
						float xl = xpos + ch.sz[0] * scale;
						float yh = y + ch.bearing[1] * scale;
						float ypos = yh - ch.sz[1] * scale;

						float vertices[24] =
						{
							xpos, yh,   0.f, 0.f,
							xpos, ypos, 0.f, 1.f,
							xl,   ypos, 1.f, 1.f,
							xpos, yh,   0.f, 0.f,
							xl,   ypos, 1.f, 1.f,
							xl,   yh,   1.f, 0.f,
						};

						glBindTexture(GL_TEXTURE_2D, ch.id);
						glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices);

						glDrawArrays(GL_TRIANGLES, 0, 6);

						FONT_DBG_GL("font::draw(const std::basic_string<CharT>&, float, float, float)");

						x += (ch.advance >> 6) * scale;
				}
		}
};

#endif // GUI_FONT_H




























