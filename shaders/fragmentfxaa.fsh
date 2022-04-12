//  Fragment shader for fast approximate anti-aliasing and dithering
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

#version 330 core

// Input (vertex shader)
in vec2 texCoord;

// Uniforms
uniform sampler2D screenTex;
uniform sampler2D blueNoiseTex;
uniform vec3 rgbMax = vec3(63);
uniform uint frame = 0u;

// FXAA parameters
uniform float relThreshold = 0.1;
uniform float absThreshold = 0.05;
uniform float mulReduce = 1./8; // 1/3 to 1/8 (best)
uniform float minReduce = 1./128; // 1/32 to 1/12 (best)
uniform float maxSpan = 8;

// Output (screen)
out vec4 color;

vec3 dither(vec3 col)
{
	const float goldenRatioReciproc = 0.61803398875;
	vec3 blueNoise = texture(blueNoiseTex, gl_FragCoord.xy / textureSize(blueNoiseTex, 0)).rgb;
	blueNoise = fract(blueNoise + ((frame & 511u) * goldenRatioReciproc));
	return floor(col*rgbMax + blueNoise)/rgbMax;
}

vec3 fxaa()
{
	vec3 rgbM = texture(screenTex, texCoord).rgb;

	vec2 texelStep = 1.f/textureSize(screenTex, 0); // textureSize() returns an ivec2!

	// Sampling neighbour texels. Offsets are adapted to OpenGL texture coordinates. 
	vec2 texCoordTrans = texCoord - 0.5f * texelStep;
	vec3 rgbSW = textureLodOffset(screenTex, texCoordTrans, 0, ivec2(0, 0)).rgb;
	vec3 rgbSE = textureLodOffset(screenTex, texCoordTrans, 0, ivec2(1, 0)).rgb;
	vec3 rgbNW = textureLodOffset(screenTex, texCoordTrans, 0, ivec2(0, 1)).rgb;
	vec3 rgbNE = textureLodOffset(screenTex, texCoordTrans, 0, ivec2(1, 1)).rgb;

	// Gather minimum and maximum.
	vec3 rgbMin = min(rgbM, min(min(rgbNW, rgbNE), min(rgbSW, rgbSE)));
	vec3 rgbMax = max(rgbM, max(max(rgbNW, rgbNE), max(rgbSW, rgbSE)));
	vec3 rgbThreshold =  max(vec3(absThreshold), rgbMax * relThreshold);

	// If contrast is lower than a maximum threshold ...
	if ( all(lessThanEqual(rgbMax - rgbMin, rgbThreshold)) )
		// ... do no AA and return.
		return rgbM;

	// see http://en.wikipedia.org/wiki/Grayscale
	const vec3 toLuma = vec3(0.299f, 0.587f, 0.114f);

	// Convert from RGB to luma.
	float lumaNW = dot(rgbNW, toLuma);
	float lumaNE = dot(rgbNE, toLuma);
	float lumaSW = dot(rgbSW, toLuma);
	float lumaSE = dot(rgbSE, toLuma);
	float lumaM = dot(rgbM, toLuma);

	float lumaMin = min(lumaM, min(min(lumaNW, lumaNE), min(lumaSW, lumaSE)));
	float lumaMax = max(lumaM, max(max(lumaNW, lumaNE), max(lumaSW, lumaSE)));

	// Sampling is done perpendicularly wrt the gradient (along the edge)
	mat2x3 rgbSamplingDir;
	rgbSamplingDir[0] = -((rgbNW + rgbNE) - (rgbSW + rgbSE));
	rgbSamplingDir[1] =  ((rgbNW + rgbSW) - (rgbNE + rgbSE));

	// Sampling step distance depends on the luma: The brighter the sampled texels, the smaller the final sampling step direction.
	// This results, that brighter areas are less blurred/more sharper than dark areas.  
	float samplingDirectionReduce = max((lumaNW + lumaNE + lumaSW + lumaSE) * 0.25f * mulReduce, minReduce);

	// Factor for norming the sampling direction plus adding the brightness influence. 
	vec3 minSamplingDirFactor = 1 / (min(abs(rgbSamplingDir[0]), abs(rgbSamplingDir[1])) + samplingDirectionReduce);

	// Calculate final sampling direction vector by reducing, clamping to a range and finally adapting to the texture size. 
	for (int i = 0; i < 2; ++i)
		rgbSamplingDir[i] = clamp(rgbSamplingDir[i] * minSamplingDirFactor, vec3(-maxSpan), vec3(maxSpan));

	vec2 samplingDir = vec2(rgbSamplingDir[0].r, rgbSamplingDir[1].r);
	vec2 samplingDir2 = vec2(rgbSamplingDir[0].g, rgbSamplingDir[1].g);
	float samplingLen = dot(samplingDir, samplingDir);
	float samplingLen2 = dot(samplingDir2, samplingDir2);
	if (samplingLen < samplingLen2)
	{
		samplingDir = samplingDir2;
		samplingLen = samplingLen2;
	}
	samplingDir2 = vec2(rgbSamplingDir[0].b, rgbSamplingDir[1].b);
	samplingLen2 = dot(samplingDir2, samplingDir2);
	if (samplingLen < samplingLen2)
		samplingDir = samplingDir2;

	samplingDir *= texelStep;

	// Inner samples on the tab.
	vec3 rgbSampleNeg = textureLod(screenTex, texCoord + samplingDir * (1/3.f - .5f), 0).rgb;
	vec3 rgbSamplePos = textureLod(screenTex, texCoord + samplingDir * (2/3.f - .5f), 0).rgb;

	vec3 rgbTwoTab = (rgbSamplePos + rgbSampleNeg) * 0.5f;  

	// Outer samples on the tab.
	rgbSampleNeg = textureLod(screenTex, texCoord + samplingDir * (0/3.f - .5f), 0).rgb;
	rgbSamplePos = textureLod(screenTex, texCoord + samplingDir * (3/3.f - .5f), 0).rgb;
	
	vec3 rgbFourTab = (rgbSamplePos + rgbSampleNeg) * 0.25f + rgbTwoTab * 0.5f;   
	
	// Calculate luma for checking against the minimum and maximum value.
	float lumaFourTab = dot(rgbFourTab, toLuma);
	
	// Are outer samples of the tab beyond the edge ... 
	if (lumaFourTab < lumaMin || lumaFourTab > lumaMax)
		// ... yes, so use only two samples.
		return rgbTwoTab; 
	else
		// ... no, so use four samples. 
		return rgbFourTab;
}

void main()
{
	color = vec4(dither(fxaa()), 1);
}

// see FXAA
// http://developer.download.nvidia.com/assets/gamedev/files/sdk/11/FXAA_WhitePaper.pdf
// http://iryoku.com/aacourse/downloads/09-FXAA-3.11-in-15-Slides.pdf
// http://horde3d.org/wiki/index.php5?title=Shading_Technique_-_FXAA
// FXAA tutorial:
// https://catlikecoding.com/unity/tutorials/advanced-rendering/fxaa/





















