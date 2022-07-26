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
#define DBG_GL(x) GLenum err_ = glGetError(); if (err_) std::cerr << err_ << " in " << x << std::endl
#else
#define DBG_GL(x)
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
			GLuint id;
			vec2i sz;
			vec2i bearing;
			GLint advance;
		};

		std::map<wchar_t, character> characters;
		mutable vec4f prev_color;

		unsigned int w, h;
		static constexpr unsigned int tab = 64, margin = 16;
		GLuint va, vb, transf_id, color_id, prog;
		bool init;
		
		void load_char(FT_Face ff, FT_ULong c)
		{
			if (FT_Load_Char(ff, c, FT_LOAD_RENDER))
			{
				std::cerr << "Cannot load character " << c << std::endl;
				return;
			}

			FT_GlyphSlot gs = ff -> glyph;

			GLuint texture;
			glGenTextures(1, &texture);
			glBindTexture(GL_TEXTURE_2D, texture);
			glTexImage2D(
				GL_TEXTURE_2D,
				0,
				GL_RED,
				gs -> bitmap.width,
				gs -> bitmap.rows,
				0,
				GL_RED,
				GL_UNSIGNED_BYTE,
				gs -> bitmap.buffer
			);

			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

			character ch = {
				texture, 
				vec2i(gs -> bitmap.width, gs -> bitmap.rows),
				vec2i(gs -> bitmap_left, gs -> bitmap_top),
				(GLint)gs -> advance.x
			};
			characters.insert(std::pair<wchar_t, character>(c, ch));
		}

		bool color_equal(float r, float g, float b, float a) const
		{
			return (r == prev_color[0] && g == prev_color[1] && b == prev_color[2] && a == prev_color[3]);
		}

		bool color_equal(const vec4f& c) const
		{
			return c == prev_color;
		}

	public:

		font() : init(false)
		{}

		font(const char* FontFile, const unsigned int font_size, const int wscreen, const int hscreen)
			: prev_color(1,1,1,1), w(wscreen), h(hscreen), init(false)
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

			FT_Error err = FT_New_Face(ft, FontFile, 0, &ff);

			if (!ff)
			{
				switch(err)
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

			for (wchar_t c = 0x0000; c < 0x0080; ++c)
				load_char(ff, c);

			for (wchar_t c = 0x0391; c < 0x03CA; ++c) // Greek characters
				load_char(ff, c);

			FT_Done_Face(ff);
			FT_Done_FreeType(ft);

			DBG_GL("font::font(const char*, unsigned int, int, int, GLuint)");

			init = true;
		}

		font(const font&) = delete;
		font& operator=(const font&) = delete;

		font(font&& f) :
			characters(std::move(f.characters)), prev_color(f.prev_color), w(f.w), h(f.h), vb(f.vb), transf_id(f.transf_id), color_id(f.color_id),
			prog(f.prog), init(f.init)
		{
			f.vb = 0;
			f.transf_id = 0;
			f.color_id = 0;
			f.prog = 0;
			f.init = false;
		}

		font& operator=(font&& f)
		{
			if (this != &f)
			{
				characters = std::move(f.characters);
				prev_color = f.prev_color;
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
		{
			glDeleteProgram(prog);
			glDeleteBuffers(1, &vb);
		}

		bool good() const
		{
			return init;
		}

		void color(float r, float g, float b, float a) const
		{
			if (!color_equal(r, g, b, a))
			{
				glUniform4f(color_id, r, g, b, a);
				prev_color = vec4f(r, g, b, a);
				DBG_GL("font::color(float, float, float, float)");
			}
		}

		void color(float r, float g, float b) const
		{
			color(r,g,b,1.f);
		}

		void color(const vec3f& col) const
		{
			color(col[0],col[1],col[2],1.f);
		}

		void color(const vec4f& col) const
		{
			if (!color_equal(col))
			{
				glUniform4fv(color_id, 1, &col[0]);
				prev_color = col;
				DBG_GL("font::color(const vec4f&)");
			}
		}

		void begin() const
		{
			glBindVertexArray(va);
			glBindBuffer(GL_ARRAY_BUFFER, vb);
			glUseProgram(prog);
			glActiveTexture(GL_TEXTURE0);

			DBG_GL("font::begin()");
		}

		void end() const
		{
			glBindTexture(GL_TEXTURE_2D, 0);
			glBindBuffer(GL_ARRAY_BUFFER, 0);
			glBindVertexArray(0);
		}

		template <typename T>
		void draw(const std::basic_string<T>& text, float x, float y, const float scale) const
		// valid for string and wstring
		{
			if (!text.length())
				return;

			const float xi = x;

			for (auto i = text.begin(); i != text.end(); ++i)
				switch (*i)
				{
					case '\n':
					{
						GLint t = characters.at('A').sz[1];
						t += characters.at('g').bearing[1];
						x = xi;
						y -= t * scale;
						break;
					}
					case '\t':
					{
						for (unsigned int i = 0; i < w - margin; i += tab)
						{
							float t = xi + i;
							if (t > x)
							{
								x = t;
								goto ok;
							}
						}
						{
							GLint t = characters.at('A').sz[1];
							t += characters.at('g').bearing[1];
							x = xi;
							y -= t * scale;
						}
						ok:
						break;
					}
					case ' ':
					{
						if (x < w - margin)
							x += (characters.at(*i).advance >> 6) * scale;
						else
						{
							GLint t = characters.at('A').sz[1];
							t += characters.at('g').bearing[1];
							x = xi;
							y -= t * scale;
						}
						break;
					}
					default:
					{
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

						DBG_GL("font::draw(const std::basic_string<T>&, float, float, float)");

						x += (ch.advance >> 6) * scale;
					}
				}
		}
};

#endif // GUI_FONT_H




























