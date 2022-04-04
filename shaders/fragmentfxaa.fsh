// FRAGMENT SHADER (FXAA + DITHER)
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

vec3 bluenoisedither(vec3 col)
{
	const float goldenRatioReciproc = 0.61803398875;
	vec3 blueNoise = texture(blueNoiseTex, gl_FragCoord.xy / textureSize(blueNoiseTex, 0)).rgb;
	blueNoise = fract(blueNoise + ((frame & 511u) * goldenRatioReciproc));
	return floor(col*rgbMax + blueNoise)/rgbMax;
}

vec3 dither(vec3 col)
{
	return bluenoisedither(col);
}

void main()
{
	vec3 rgbM = texture(screenTex, texCoord).rgb;

	vec2 texelStep = 1./textureSize(screenTex, 0);

	// Sampling neighbour texels. Offsets are adapted to OpenGL texture coordinates. 
	vec2 texCoordTrans = texCoord - 0.5f * texelStep;
	vec3 rgbSW = textureLodOffset(screenTex, texCoordTrans, 0, ivec2(0, 0)).rgb;
	vec3 rgbSE = textureLodOffset(screenTex, texCoordTrans, 0, ivec2(1, 0)).rgb;
	vec3 rgbNW = textureLodOffset(screenTex, texCoordTrans, 0, ivec2(0, 1)).rgb;
	vec3 rgbNE = textureLodOffset(screenTex, texCoordTrans, 0, ivec2(1, 1)).rgb;

	// see http://en.wikipedia.org/wiki/Grayscale
	const vec3 toLuma = vec3(0.299f, 0.587f, 0.114f);

	// Convert from RGB to luma.
	float lumaNW = dot(rgbNW, toLuma);
	float lumaNE = dot(rgbNE, toLuma);
	float lumaSW = dot(rgbSW, toLuma);
	float lumaSE = dot(rgbSE, toLuma);
	float lumaM = dot(rgbM, toLuma);

	// Gather minimum and maximum luma.
	float lumaMin = min(lumaM, min(min(lumaNW, lumaNE), min(lumaSW, lumaSE)));
	float lumaMax = max(lumaM, max(max(lumaNW, lumaNE), max(lumaSW, lumaSE)));
	float lumaThreshold = max(absThreshold, lumaMax * relThreshold);

	// If contrast is lower than a maximum threshold ...
	if (lumaMax - lumaMin <= lumaThreshold)
	{
		// ... do no AA and return.
		color = vec4(dither(rgbM), 1);
		
		return;
	}

	// Sampling is done perpendicularly wrt the gradient (along the edge)
	vec2 samplingDirection;	
	samplingDirection.x = -((lumaNW + lumaNE) - (lumaSW + lumaSE));
	samplingDirection.y =  ((lumaNW + lumaSW) - (lumaNE + lumaSE));

	// Sampling step distance depends on the luma: The brighter the sampled texels, the smaller the final sampling step direction.
	// This results, that brighter areas are less blurred/more sharper than dark areas.  
	float samplingDirectionReduce = max((lumaNW + lumaNE + lumaSW + lumaSE) * 0.25f * mulReduce, minReduce);

	// Factor for norming the sampling direction plus adding the brightness influence. 
	float minSamplingDirectionFactor = 1 / (min(abs(samplingDirection.x), abs(samplingDirection.y)) + samplingDirectionReduce);

	// Calculate final sampling direction vector by reducing, clamping to a range and finally adapting to the texture size. 
	samplingDirection = clamp(samplingDirection * minSamplingDirectionFactor, vec2(-maxSpan), vec2(maxSpan)) * texelStep;
	
	// Inner samples on the tab.
	vec3 rgbSampleNeg = textureLod(screenTex, texCoord + samplingDirection * (1/3.f - .5f), 0).rgb;
	vec3 rgbSamplePos = textureLod(screenTex, texCoord + samplingDirection * (2/3.f - .5f), 0).rgb;

	vec3 rgbTwoTab = (rgbSamplePos + rgbSampleNeg) * 0.5f;  

	// Outer samples on the tab.
	vec3 rgbSampleNegOuter = textureLod(screenTex, texCoord - samplingDirection * (0/3.f - .5f), 0).rgb;
	vec3 rgbSamplePosOuter = textureLod(screenTex, texCoord - samplingDirection * (3/3.f - .5f), 0).rgb;
	
	vec3 rgbFourTab = (rgbSamplePosOuter + rgbSampleNegOuter) * 0.25f + rgbTwoTab * 0.5f;   
	
	// Calculate luma for checking against the minimum and maximum value.
	float lumaFourTab = dot(rgbFourTab, toLuma);
	
	// Are outer samples of the tab beyond the edge ... 
	if (lumaFourTab < lumaMin || lumaFourTab > lumaMax)
		// ... yes, so use only two samples.
		color = vec4(dither(rgbTwoTab), 1.0f); 
	else
		// ... no, so use four samples. 
		color = vec4(dither(rgbFourTab), 1.0f);
}

// see FXAA
// http://developer.download.nvidia.com/assets/gamedev/files/sdk/11/FXAA_WhitePaper.pdf
// http://iryoku.com/aacourse/downloads/09-FXAA-3.11-in-15-Slides.pdf
// http://horde3d.org/wiki/index.php5?title=Shading_Technique_-_FXAA
// FXAA tutorial:
// https://catlikecoding.com/unity/tutorials/advanced-rendering/fxaa/





















