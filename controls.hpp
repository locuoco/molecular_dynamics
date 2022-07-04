//  Controls
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

#ifndef CONTROLS_H
#define CONTROLS_H

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "math/helper.hpp"

namespace controls
{
	physics::mat4d Proj;
	physics::mat4d View;
	physics::point3d pos(0, 0, 5);
	physics::point3d vel(0, 0, 0);
	double horVel = 0;
	double vertVel = 0;
	double horAngle = math::pi<double>();
	double vertAngle = 0;
	double FoV = 90;
	double gamma = 2.2+1./30;
	bool cursor_disabled = false;
}

inline void updateScroll(GLFWwindow*, double, double yoff)
{
	using namespace controls;

	FoV -= 5 * yoff;

	FoV = math::clamp(FoV, 30., 110.);
}

inline void updateControls(GLFWwindow* window, int w, int h, double dt)
{
	using namespace controls;
	using namespace math;

	constexpr double acc = 40*5, strengthSpeed = 0.5;
	constexpr double mouseAcc = 5.e-2;
	double friction = std::min(6., 1/dt), mouseFriction = std::min(17., 1/dt);

	double xpos, ypos;

	glfwGetCursorPos(window, &xpos, &ypos);

	bool rightBPressed = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS;

	if (rightBPressed && !cursor_disabled)
	{
		glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

		glfwSetCursorPos(window, w/2, h/2);

		cursor_disabled = true;
	}
	else if (rightBPressed && cursor_disabled)
	{
		horVel += mouseAcc * (w/2 - xpos);
		vertVel += mouseAcc * (h/2 - ypos);

		glfwSetCursorPos(window, w/2, h/2);
	}
	else if (!rightBPressed && cursor_disabled)
	{
		glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);

		cursor_disabled = false;
	}

	horAngle += horVel * dt;
	vertAngle += vertVel * dt;

	horVel -= horVel * mouseFriction * dt;
	vertVel -= vertVel * mouseFriction * dt;

	if (horAngle > pi<double>())
		horAngle -= two_pi<double>();
	if (horAngle < -pi<double>())
		horAngle += two_pi<double>();

	constexpr double maxVertAngle = half_pi<double>() * (1.-eps<float>());

	if (vertAngle < -maxVertAngle)
	{
		vertAngle = -maxVertAngle;
		vertVel = 0;
	}
	if (vertAngle > maxVertAngle)
	{
		vertAngle = maxVertAngle;
		vertVel = 0;
	}

	physics::point3d dir(
		cos(vertAngle) * sin(horAngle),
		sin(vertAngle),
		cos(vertAngle) * cos(horAngle)
	);

	physics::point3d horDir(
		sin(horAngle),
		0,
		cos(horAngle)
	);

	physics::point3d up(0,1,0);
	physics::point3d right = normalize(cross(dir, up));

	vel -= vel * (friction * dt);

	if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
		vel += horDir * (dt * acc);
	if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
		vel -= horDir * (dt * acc);
	if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
		vel += right * (dt * acc);
	if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
		vel -= right * (dt * acc);
	if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
		vel += up * (dt * acc);
	if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
		vel -= up * (dt * acc);

	pos += vel * dt;

	if (glfwGetKey(window, GLFW_KEY_G) == GLFW_PRESS)
	{
		if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS)
			gamma += dt * strengthSpeed;
		if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS)
			gamma -= dt * strengthSpeed;

		gamma = clamp(gamma, 1., 5.);
	}

	double aspect = (double)w / (double)h;

	double height = tan(deg2rad(FoV)/2);

	double f = 1 / height;

    Proj = physics::mat4d(
		f / aspect, 0, 0, 0,
		0, f, 0, 0,
		0, 0, 0, -1,
		0, 0, 1.e-4, 0
	); // perspective projection with infinite far plane for reverse depth buffer
	transpose(Proj);
	View = physics::look_at(
		pos,		// position
		pos + dir,	// target
		up			// up vector
	);
}

#endif
































