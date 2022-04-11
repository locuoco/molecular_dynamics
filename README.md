# Molecular dynamics 

![Screenshot of the program](screenshot.png)

This is an interactive program for molecular dynamics using classical force fields. A force field is defined from the following potential (using CHARMM convention):
<!---
$V(\mathrm{r})&=\sum_{i \sim j}K_{ij}^r \left(r_{ij} - r_{ij}^0\right)^2
+ \sum_{i\sim j\sim k}K_{ijk}^{\theta} \left(\theta_{ijk} - \theta_{ijk}^0\right)^2
+ \sum_{i\sim\cdot\sim k}K_{ik}^{UB} \left(r_{ik} - r_{ik}^0\right)^2 \\
&+ \sum_{i\sim j\sim k\sim l, n}K_{ijkl}^{\chi} \left(1 + \cos\left(n\chi_{ijkl} - \delta_{ijkl}\right) \right)
+ \sum_{ijkl\ \text{impropers}}K_{ijkl}^{\psi} \left( \psi_{ijkl} - \psi_{ijkl}^0 \right)^2 \\
&+ \sum_{ij} \epsilon_{ij} \left( \left( \frac{R_{ij}}{r_{ij}} \right)^{12} - 2 \left(\frac{R_{ij}}{r_{ij}} \right)^{6} \right)
+ \sum_{ij} k_C \frac{q_i q_j}{r_{ij}^2}$
--->
<img src="https://render.githubusercontent.com/render/math?math=\displaystyle+V(\mathrm{r})=\sum_{i\sim+j}K_{ij}^r\left(r_{ij}-r_{ij}^0\right)^2%2B\sum_{i\sim+j\sim+k}K_{ijk}^{\theta} \left(\theta_{ijk}-\theta_{ijk}^0\right)^2%2B\sum_{i\sim\cdot\sim+k}K_{ik}^{UB}\left(r_{ik}-r_{ik}^0\right)^2\\">
<img src="https://render.githubusercontent.com/render/math?math=\displaystyle\qquad%2B\sum_{i\sim+j\sim+k\sim+l,n}K_{ijkl}^{\chi}+\left(1%2B\cos\left(n\chi_{ijkl}-\delta_{ijkl}\right)\right)%2B\sum_{ijkl\+\text{impropers}}K_{ijkl}^{\psi}\left(\psi_{ijkl}-\psi_{ijkl}^0\right)^2\\%2B\sum_{ij}\epsilon_{ij}\left(\left(\frac{R_{ij}}{r_{ij}}\right)^{12}-2\left(\frac{R_{ij}}{r_{ij}}\right)^{6}\right)%2B\sum_{ij}k_C\frac{q_i+q_j}{r_{ij}^2}">

The first term describes the bond potential, modeled as an elastic potential between bonded atoms (with elastic constant <img src="https://render.githubusercontent.com/render/math?math=k^r=2K^r">), the second term is the angle potential, the third term is the Urey-Bradley potential, acting between 1-3 atoms, the fourth and fifth terms describe the dihedral (proper and improper) angles potentials and the last two terms correspond to non-bonded potentials (Lennard-Jones and electrostatic respectively).

References:
* M. P. Allen, D. J. Tildesley, *Computer Simulation of Liquids*, Oxford University Press, 2017
* P. Spijker, B. Markvoort, P. Hilbers, *Parallel Utility for Modeling of Molecular Aggregation*, Biomodeling and bioinformatics, Eindhoven University of Technology, 2007

![Animation of a system of water molecules](animation.gif)

## The code
The code makes use of C++ templates and concepts (thus it requires C++20) and is organised in many header files, that can be included from a single compilation unit. It is organised in the following way:
* `physics`: directory which contains part of code relevant to the resolution of the physical/numerical problem.
  * `fmm.hpp`: fast multipole method for fast calculation of long-range forces (not implemented yet).
  * `integrator.hpp`: classes and concepts for integration of classical hamiltonian dynamical systems.
  * `molecule.hpp`: classes, methods and other utilites for managing molecular systems (for now, you have only water molecules, also the dihedral potential part must be added).
  * `physics.hpp`: mainly a header file that includes everything.
  * `point.hpp`: classes, aliases and data structures for use in generic dynamical systems.
* `shaders`: directory which contains shaders for the rendering of impostors (for fast rendering of spheres), post-processing filters (fast approximate anti-aliasing, blue noise dithering) and text.
* `controls.hpp`: includes a function to manage keyboard and mouse controls.
* `Font.hpp`: class for drawing text on the window.
* `graphics.hpp`: class that manages the graphics used in this specific program.
* `math.hpp`: some math helper functions.
* `shader.hpp`: some functions to load shaders from file.
* `main.cpp`: it just contains the main loop and some basic initialization of the system (it can be ignored or modified).

## Dependencies

The graphics part has some dependencies on public libraries that enable the usage of modern OpenGL (Open Graphics Library):
* [GLFW 3](https://www.glfw.org/) (Graphics Library Framework 3): an open source, multi-platform API for creating windows, contexts and managing input and events.
* [GLEW](http://glew.sourceforge.net/) (OpenGL Extension Wrangler): a cross-platform library that include core and extended OpenGL functionalities.
* [FreeType](https://freetype.org/): an OpenGL library to render fonts (used in `Font.hpp`).
* [GLM](https://glm.g-truc.net/0.9.9/) (OpenGL Mathematics): a simple math header-only library for convenient interface with OpenGL and GLSL (OpenGL Shading Language).
* [STB Image](https://github.com/nothings/stb/blob/master/stb_image.h): a single-header file, public domain library for loading images.

These dependencies are required only for these header files:
* `graphics.hpp`
* `Font.hpp` (included in `graphics.hpp`)
* `shader.hpp` (included in `graphics.hpp`)

## Compilation

To compile the program, simply do (with MinGW):

    g++ main.cpp -o mold -std=c++20 -I <includes> -L <libs> -lopengl32 -lglu32 -lglew32.dll -lglfw3dll -lfreetype -Wall -Wextra -pedantic -O3

On Linux, the library names could be different:

    g++ main.cpp -o mold -std=c++20 -I <includes> -L <libs> -lopengl -lglu -lglew -lglfw3 -lfreetype -Wall -Wextra -pedantic -O3

where `<includes>` and `<libs>` are the paths for installed libraries header files and static library files (if required). The executable will be called `mold`.
