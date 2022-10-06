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

#ifndef GUI_CONTROLS_H
#define GUI_CONTROLS_H

#include <numbers> // pi
#include <limits> // numeric_limits<T>::epsilon
#include <cmath> // sin, cos, exp

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "../math/helper.hpp" // two_pi, half_pi, clamp, deg2rad

inline void scroll_callback(GLFWwindow*, double, double);
// see below definition for documentation

struct controls
// class for handling mouse and key input and updating the camera
// position and direction and other parameters (like gamma correction).
// Used in `graphics.hpp`.
{
	physics::mat4d proj; // (perspective) projection matrix
	physics::mat4d view; // view matrix
	physics::vec3d pos = physics::vec3d(0, 0, 5); // position of the camera
	physics::vec3d vel = 0; // velocity of the camera
	physics::vec3d dir; // direction of the camera
	physics::vec3d up = physics::vec3d(0, 1, 0); // up-direction vector of the camera
	physics::vec2d angles = physics::vec2d(std::numbers::pi, 0); // (spherical) angles of the camera: (azimuthal, polar)
	physics::vec2d ang_vel = 0; // angular velocity of the camera
	physics::vec2d last_cursor_pos = 0; // last coordinates of the cursor
	GLFWwindow* window;
	double fov = 90; // field of view angle (in degrees)
	double gamma = 2.2+1./30; // gamma correction exponent
	int w, h; // width and height of the window
	bool cursor_disabled = false;

	friend void scroll_callback(GLFWwindow*, double, double);
	// see below definition for documentation

	controls(GLFWwindow* window, int w, int h) : window(window), w(w), h(h)
	// constructor:
	//	`window` is a pointer to GLFWwindow returned by `glfwCreateWindow`.
	//	`w` and `h` are the width and height of the window respectively.
	{
		glfwPollEvents();

		// enable keyboard input
		glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
		// set user pointer to `this`, which will be used inside the
		// `scroll_callback` function
		glfwSetWindowUserPointer(window, this);
		// set `scroll_callback` as the callback managing the input from the scroll wheel
		glfwSetScrollCallback(window, scroll_callback);
	}

	void update(double dt)
	// update camera position and rotation
	// `dt` is the real-time step
	{
		// collect unhandled mouse and key events
		glfwPollEvents();

		update_mouse(dt);
		update_keys(dt);
		update_view_proj();
	}

	private:

		void update_fov(double yoff)
		// update the value of `fov` (field of view).
		// `yoff` measures how much the mouse wheel has rotated and in which direction (upward or downward).
		// This method is called from friend function `scroll_callback`.
		{
			fov -= 5 * yoff;
			fov = math::clamp(fov, 30., 110.);
		}

		void update_mouse(double dt)
		// update rotation of the camera from mouse input.
		// `dt` is the real-time step.
		{
			using std::sin;
			using std::cos;
			using std::exp;

			constexpr double angle_acc_coeff = 5.e-2; // scaling applied to camera angular acceleration when the cursor moved
			constexpr double angleFriction = 17.; // friction coefficient applied to camera angular components

			// angular acceleration of the camera
			physics::vec2d ang_acc = 0;
			// get cursor position
			physics::vec2d cursor_pos;
			glfwGetCursorPos(window, &cursor_pos[0], &cursor_pos[1]);

			bool mouse_right_pressed = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS;

			// if right button is being pressed and cursor has not yet been disabled...
			if (mouse_right_pressed && !cursor_disabled)
			{
				// ...disable (hide) cursor
				glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

				cursor_disabled = true;
			}
			// if right button is being pressed and cursor has already been disabled...
			else if (mouse_right_pressed && cursor_disabled)
				// ...calculate the camera angular acceleration from the displacement of the cursor
				ang_acc = angle_acc_coeff * (last_cursor_pos - cursor_pos) / dt;
			// if right button is not being pressed and cursor is disabled...
			else if (!mouse_right_pressed && cursor_disabled)
			{
				// ...enable (show) cursor and set its position to the center of the window
				glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
				glfwSetCursorPos(window, w/2, h/2);

				cursor_disabled = false;
			}
			last_cursor_pos = cursor_pos;

			// simple integrator with damping
			ang_vel += ang_acc * dt;
			angles += ang_vel * (dt/2);
			ang_vel *= exp(-angleFriction * dt);
			angles += ang_vel * (dt/2);

			// remapping horizontal (azimuthal) angle into the [-pi, pi] range
			while (angles[0] > std::numbers::pi)
				angles[0] -= math::two_pi<double>;
			while (angles[0] < -std::numbers::pi)
				angles[0] += math::two_pi<double>;

			constexpr double max_vert_angle = math::half_pi<double> * (1.-std::numeric_limits<float>::epsilon());

			// clamping vertical (polar) angle between [-pi/2 + eps, pi/2 - eps], with eps = small number
			// setting polar angular velocity to zero if the angle has been clamped
			if (angles[1] < -max_vert_angle)
			{
				angles[1] = -max_vert_angle;
				ang_vel[1] = 0;
			}
			if (angles[1] > max_vert_angle)
			{
				angles[1] = max_vert_angle;
				ang_vel[1] = 0;
			}

			// set direction of the camera
			dir = physics::vec3d(
				cos(angles[1]) * sin(angles[0]),
				sin(angles[1]),
				cos(angles[1]) * cos(angles[0])
			);
		}

		void update_keys(double dt)
		// update translation of the camera from key input.
		// `dt` is the real-time step.
		{
			using std::sin;
			using std::cos;
			using std::exp;

			constexpr double acc_coeff = 200; // scaling applied to camera linear acceleration when a key is pressed
			constexpr double friction = 6.; // friction coefficient applied to camera linear components
			constexpr double gammaSpeed = 0.5; // rate of change of gamma value when 'G'+up/down are pressed

			// linear acceleration of the camera
			physics::vec3d acc = 0;

			// forward-direction vector
			physics::vec3d forward(
				sin(angles[0]),
				0,
				cos(angles[0])
			);
			// right-direction vector
			physics::vec3d right = normalize(cross(dir, up));
			// y is the vertical axis, so `forward` and `right` directions
			// always lie along the xz plane

			// if 'W' is pressed, go forward
			if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
				acc += forward * acc_coeff;
			// if 'S' is pressed, go backward
			if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
				acc -= forward * acc_coeff;
			// if 'D' is pressed, go right
			if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
				acc += right * acc_coeff;
			// if 'A' is pressed, go left
			if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
				acc -= right * acc_coeff;
			// if 'space' is pressed, go up
			if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
				acc += up * acc_coeff;
			// if 'left shift' is pressed, go down
			if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
				acc -= up * acc_coeff;

			// simple integrator with damping
			vel += acc * dt;
			pos += vel * (dt/2);
			vel *= exp(-friction * dt);
			pos += vel * (dt/2);

			// if 'G' is pressed...
			if (glfwGetKey(window, GLFW_KEY_G) == GLFW_PRESS)
			{
				// ...and 'up' is pressed, increase gamma correction exponent
				if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS)
					gamma += dt * gammaSpeed;
				// ...and 'down' is pressed, decrease gamma correction exponent
				if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS)
					gamma -= dt * gammaSpeed;

				gamma = math::clamp(gamma, 1., 5.);
			}
		}

		void update_view_proj()
		// update view and projection matrices
		{
			double aspect = (double)w / h; // aspect ratio
			double height = tan(math::deg2rad(fov)/2);
			double yscal = 1 / height;

			// calculate perspective projection matrix with infinite far plane
			// for reverse depth buffer
			proj = physics::mat4d(
				yscal / aspect, 0, 0, 0,
				0, yscal, 0, 0,
				0, 0, 0, 1.e-6,
				0, 0, -1, 0
			);
			view = physics::look_at(
				pos,       // position
				pos + dir, // target
				up         // up vector
			);
		}
};

inline void scroll_callback(GLFWwindow* w, double, double yoff)
// a function wrapper that calls the method `update_fov` of a `controls` object,
// whose pointer was stored with `glfwSetWindowUserPointer`.
// `w` is a pointer to a valid GLFW window object.
// `yoff` argument is simply forwarded to `update_fov`.
// `glfwGetWindowUserPointer(w)` retrieves the pointer to the `controls` object (as a `void*`).
{
	static_cast<controls*>(glfwGetWindowUserPointer(w))->update_fov(yoff);
}

#endif // GUI_CONTROLS_H
































