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

#ifdef ENABLE_TEST
#define TEST(x) GLenum _err_ = glGetError(); if (_err_) std::cerr << _err_ << " in " << x << std::endl
#else
#define TEST(x)
#endif

#include <iostream> // cerr, endl
#include <string>
#include <map>
#include <utility> // pair

#include <GL/glew.h>
#include <ft2build.h> // Compile with -lfreetype
#include FT_FREETYPE_H

#ifdef USE_POINT

#include "../physics/point.hpp"

#else

#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>
#include <glm/mat4x4.hpp>
#include <glm/trigonometric.hpp>
#include <glm/gtc/matrix_transform.hpp>

#endif

class Font
{
#ifdef USE_POINT
		using vec2i = physics::point2i;
		using vec3f = physics::point3f;
		using vec4f = physics::point4f;
		using mat4f = physics::mat4f;
#else
		using vec2i = glm::ivec2;
		using vec3f = glm::vec3;
		using vec4f = glm::vec4;
		using mat4f = glm::mat4;
#endif

		struct Character {
			GLuint ID;
			vec2i Size;
			vec2i Bearing;
			GLint Advance;
		};

		std::map<wchar_t, Character> characters;
		mutable vec4f precCol;

		unsigned int w, h;
		static constexpr unsigned int tab = 64, margin = 16;
		GLuint vb, transf, col, prog;
		char init;
		
		void LoadChar(FT_Face ff, FT_ULong c)
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

			Character ch = {
				texture, 
				vec2i(gs -> bitmap.width, gs -> bitmap.rows),
				vec2i(gs -> bitmap_left, gs -> bitmap_top),
				(GLint)gs -> advance.x
			};
			characters.insert(std::pair<wchar_t, Character>(c, ch));
		}

		bool coleq(float r, float g, float b, float a) const
		{
			if (r == precCol[0] && g == precCol[1] && b == precCol[2] && a == precCol[3])
				return true;
			return false;
		}
		bool coleq(const vec4f& c) const
		{
			if (c == precCol)
				return true;
			return false;
		}

	public:

		Font() : init(false)
		{

		}

		Font(const char* FontFile, const unsigned int font_size, const int wscreen, const int hscreen)
			: precCol(1,1,1,1), w(wscreen), h(hscreen), init(false)
		{
			prog = loadShader("gui/shaders/vertexfont.vsh", "gui/shaders/fragmentfont.fsh");

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

			glGenBuffers(1, &vb);

			glBindBuffer(GL_ARRAY_BUFFER, vb);
			glBufferData(GL_ARRAY_BUFFER, 24*sizeof(GLfloat), nullptr, GL_DYNAMIC_DRAW);

			transf = glGetUniformLocation(prog, "Pmat");
			col = glGetUniformLocation(prog, "col");

			glUseProgram(prog);
#ifdef USE_POINT
			mat4f m_Pmat = physics::orthographic_projection(0., double(w), 0., double(h));
			glUniformMatrix4fv(transf, 1, GL_TRUE, &m_Pmat(0, 0));
#else
			mat4f m_Pmat = glm::ortho(0., double(w), 0., double(h));
			glUniformMatrix4fv(transf, 1, GL_FALSE, &m_Pmat[0][0]);
#endif

			glPixelStorei(GL_UNPACK_ALIGNMENT, GL_TRUE);

			for (wchar_t c = 0x0000; c < 0x0080; ++c)
				LoadChar(ff, c);

			for (wchar_t c = 0x0391; c < 0x03CA; ++c) // Greek characters
				LoadChar(ff, c);

			FT_Done_Face(ff);
			FT_Done_FreeType(ft);

			TEST("Font::Font(const char*, unsigned int, int, int, GLuint)");

			init = true;
		}

		Font(const Font&) = delete;
		Font& operator=(const Font&) = delete;

		Font(Font&& font) :
			characters(std::move(font.characters)), precCol(font.precCol), w(font.w), h(font.h), vb(font.vb), transf(font.transf), col(font.col),
			prog(font.prog), init(font.init)
		{
			font.vb = 0;
			font.transf = 0;
			font.col = 0;
			font.prog = 0;
			font.init = false;
		}

		Font& operator=(Font&& font)
		{
			if (this != &font)
			{
				characters = std::move(font.characters);
				precCol = font.precCol;
				w = font.w;
				h = font.h;
				vb = font.vb;
				transf = font.transf;
				col = font.col;
				prog = font.prog;
				init = font.init;
				font.vb = 0;
				font.transf = 0;
				font.col = 0;
				font.prog = 0;
				font.init = false;
			}
			return *this;
		}

		~Font()
		{
			glDeleteProgram(prog);
			glDeleteBuffers(1, &vb);
		}

		bool good() const
		{
			return (bool)init;
		}

		void Color(float r, float g, float b, float a) const
		{
			if (!coleq(r, g, b, a))
			{
				glUniform4f(col, r, g, b, a);
				precCol = vec4f(r, g, b, a);
				TEST("Font::Color(float, float, float, float)");
			}
		}

		void Color(float r, float g, float b) const
		{
			Color(r,g,b,1.f);
		}

		void Color(const vec3f& color) const
		{
			Color(color[0],color[1],color[2],1.f);
		}

		void Color(const vec4f& color) const
		{
			if (!coleq(color))
			{
				glUniform4fv(col, 1, &color[0]);
				precCol = color;
				TEST("Font::Color(const vec4f&)");
			}
		}

		void Begin() const
		{
			glUseProgram(prog);
			glActiveTexture(GL_TEXTURE0);

			glEnableVertexAttribArray(0);
			glBindBuffer(GL_ARRAY_BUFFER, vb);
			glVertexAttribPointer(
				0,
				4,
				GL_FLOAT,
				GL_FALSE,
				0,
				(void*)0
			);
			TEST("Font::Begin()");
		}

		void End() const
		{
			glBindTexture(GL_TEXTURE_2D, 0);
			glDisableVertexAttribArray(0);
		}

		template<typename T>
		void Draw(const std::basic_string<T>& text, float x, float y, const float scale) const
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
						GLint t = characters.at('A').Size[1];
						t += characters.at('g').Bearing[1];
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
							GLint t = characters.at('A').Size[1];
							t += characters.at('g').Bearing[1];
							x = xi;
							y -= t * scale;
						}
						ok:
						break;
					}
					case ' ':
					{
						if (x < w - margin)
							x += (characters.at(*i).Advance >> 6) * scale;
						else
						{
							GLint t = characters.at('A').Size[1];
							t += characters.at('g').Bearing[1];
							x = xi;
							y -= t * scale;
						}
						break;
					}
					default:
					{
						Character ch = characters.at(*i);

						float xpos = x + ch.Bearing[0] * scale;
						float xl = xpos + ch.Size[0] * scale;
						float yh = y + ch.Bearing[1] * scale;
						float ypos = yh - ch.Size[1] * scale;

						float vertices[24] =
						{
							xpos, yh,   0.f, 0.f,
							xpos, ypos, 0.f, 1.f,
							xl,   ypos, 1.f, 1.f,
							xpos, yh,   0.f, 0.f,
							xl,   ypos, 1.f, 1.f,
							xl,   yh,   1.f, 0.f,
						};

						glBindTexture(GL_TEXTURE_2D, ch.ID);
						glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices);

						glDrawArrays(GL_TRIANGLES, 0, 6);

						TEST("Font::Draw(const std::basic_string<T>&, float, float, float)");

						x += (ch.Advance >> 6) * scale;
					}
				}
		}
};

#endif // GUI_FONT_H




























