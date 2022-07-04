# Molecular dynamics 

![Screenshot of the program](screenshot.png)

This is an interactive program for molecular dynamics using classical force fields. A force field is defined from the following potential (using CHARMM convention):
<!---
$V(\mathrm{r}_i)&=\sum_{i \sim j}K_{ij}^r \left(r_{ij} - r_{ij}^0\right)^2
+ \sum_{i\sim j\sim k}K_{ijk}^{\theta} \left(\theta_{ijk} - \theta_{ijk}^0\right)^2
+ \sum_{i\sim\cdot\sim k}K_{ik}^{UB} \left(r_{ik} - r_{ik}^0\right)^2 \\
&+ \sum_{i\sim j\sim k\sim l, n}K_{ijkl}^{\chi} \left(1 + \cos\left(n\chi_{ijkl} - \delta_{ijkl}\right) \right)
+ \sum_{ijkl\ \text{impropers}}K_{ijkl}^{\psi} \left( \psi_{ijkl} - \psi_{ijkl}^0 \right)^2 \\
&+ \sum_{ij} \epsilon_{ij} \left( \left( \frac{R_{ij}}{r_{ij}} \right)^{12} - 2 \left(\frac{R_{ij}}{r_{ij}} \right)^{6} \right)
+ \sum_{ij} k_C \frac{q_i q_j}{r_{ij}^2}$

Old URL:
https://render.githubusercontent.com/render/math?math=\displaystyle+

Plus (+) = %2B
Space = %20
Comma (,) = %2C
--->
<img src="https://latex.codecogs.com/gif.latex?V=\sum_{i\sim%20j}K_{ij}^r\left(r_{ij}-r_{ij}^0\right)^2%2B\sum_{i\sim%20j\sim%20k}K_{ijk}^{\theta}%20\left(\theta_{ijk}-\theta_{ijk}^0\right)^2%2B\sum_{i\sim\cdot\sim%20k}K_{ik}^{UB}\left(r_{ik}-r_{ik}^0\right)^2\\">
<img src="https://latex.codecogs.com/gif.latex?\qquad%2B\sum_{i\sim%20j\sim%20k\sim%20l%2Cn}K_{ijkl}^{\chi}%20\left(1%2B\cos\left(n\chi_{ijkl}-\delta_{ijkl}\right)\right)%2B\sum_{ijkl\%20\text{impropers}}K_{ijkl}^{\psi}\left(\psi_{ijkl}-\psi_{ijkl}^0\right)^2\\">
<img src="https://latex.codecogs.com/gif.latex?\qquad%2B\sum_{i<j}\epsilon_{ij}\left(\left(\frac{R_{ij}}{r_{ij}}\right)^{12}-2\left(\frac{R_{ij}}{r_{ij}}\right)^6\right)%2B\sum_{i<j}k_C\frac{q_i%20q_j}{r_{ij}^2}">

The first term describes the bond potential, modeled as an elastic potential between bonded atoms (with elastic constant <img src="https://render.githubusercontent.com/render/math?math=k^r=2K^r">), the second term is the angle potential, the third term is the Urey-Bradley potential, acting between 1-3 atoms, the fourth and fifth terms describe the dihedral (proper and improper) angles potentials and the last two terms correspond to non-bonded potentials (Lennard-Jones and electrostatic respectively).

References:
* M. P. Allen, D. J. Tildesley, *Computer Simulation of Liquids*, Oxford University Press, 2017
* P. Spijker, B. Markvoort, P. Hilbers, *Parallel Utility for Modeling of Molecular Aggregation*, Biomodeling and bioinformatics, Eindhoven University of Technology, 2007

![Animation of a system of water molecules](animation.gif)

### Ewald summation

In a cubic periodic system of side <img src="https://latex.codecogs.com/gif.latex?L">, the electrostatic potential is given by:

<img src="https://latex.codecogs.com/gif.latex?V=\frac{1}{2}\sum_{ij}^N\sum_{\mathrm{n}\in%20Z^3}^{%27}\frac{q_i%20q_j}{\left|\mathrm{r}_{ij}%2B\mathrm{n}L\right|}">

where the prime symbol means that the <img src="https://latex.codecogs.com/gif.latex?i=j"> term must be excluded for <img src="https://latex.codecogs.com/gif.latex?\mathrm{n}=\mathrm{0}">. As it is not practical to calculate all the contributions directly, it is more convenient to calculate the long-range part in Fourier space, leading to a formula which converges much faster than the previous equation. To do so, the Green's function must be separated in the following way:

<img src="https://latex.codecogs.com/gif.latex?\frac{1}{r}=\frac{f(r)}{r}%2B\frac{1-f(r)}{r}">

where <img src="https://latex.codecogs.com/gif.latex?f(r)"> is the known as the splitting function. Choosing <img src="https://latex.codecogs.com/gif.latex?f(r)"> to be the complementary error function, one obtains the classical Ewald summation, whose terms are given by:

<img src="https://latex.codecogs.com/gif.latex?V=V^{(r)}%2BV^{(k)}%2BV^{(s)}%2BV^{(d)}">
<img src="https://latex.codecogs.com/gif.latex?V^{(r)}=\frac{1}{2}\sum_{ij}^N\sum_{\mathrm{n}\in%20Z^3}^{%27}q_i%20q_j\frac{\text{erfc}\left(\kappa\left|\mathrm{r}_{ij}%2B\mathrm{n}L\right|\right)}{\left|\mathrm{r}_{ij}%2B\mathrm{n}L\right|}">
<img src="https://latex.codecogs.com/gif.latex?V^{(k)}=\frac{1}{2L^3}\sum_{\mathrm{n}\neq\mathrm{0}}\frac{4\pi}{k_n^2}e^{-k_n^2/4\kappa^2}\left|\tilde{\rho}(\mathrm{k}_n)\right|^2">

where <img src="https://latex.codecogs.com/gif.latex?V^{(r)}"> ... etc...

### PPPM method

TODO

References:
* M. Deserno, C. Holm, *How to mesh up Ewald sums (I): A theoretical and numerical comparison of various particle mesh routines*, Max-Planck-Institut fur Polymerforschung, Ackermannweg, Germany, 1998

![Animation of a system of water molecules](animation.gif)

## The code
The code makes use of C++ templates and concepts (thus it requires C++20) and is organised in many header files, that can be included from a single compilation unit. It is organised in the following way:
* `physics`: directory which contains part of code relevant to the resolution of the physical/numerical problem.
  * `fmm.hpp`: fast multipole method for fast calculation of long-range forces (not implemented yet).
  * `integrator.hpp`: classes and concepts for integration of classical (hamiltonian and not) dynamical systems.
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
* [STB Image](https://github.com/nothings/stb/blob/master/stb_image.h): a single-header file, public domain library for loading images.

These dependencies are required only for these header files:
* `graphics.hpp`
* `Font.hpp` (included in `graphics.hpp`)
* `shader.hpp` (included in `graphics.hpp`)
* `controls.hpp` (included in `graphics.hpp`)

## Compilation

To compile the program, simply do (with MinGW):

    g++ main.cpp -o mold -std=c++20 -I <includes> -L <libs> -lopengl32 -lglu32 -lglew32.dll -lglfw3dll -lfreetype -Wall -Wextra -pedantic -O3

On Linux, the library names could be different:

    g++ main.cpp -o mold -std=c++20 -I <includes> -L <libs> -lopengl -lglu -lglew -lglfw3 -lfreetype -Wall -Wextra -pedantic -O3

where `<includes>` and `<libs>` are the paths for installed libraries header files and static library files (if required). The executable will be called `mold`.

## Basic usage

To create a molecular system use the following:
```c++
#include "physics/physics.hpp"

int main()
{
    physics::molecular_system my_system;
}
```
To add a molecule to the system, use the method `add_molecule`:
```c++
    my_system.add_molecule(molecule);
```
Currently available molecules are:

* `physics::water_tip3p<>`: Water molecule, using a flexible TIP3P model with O–H and H–H Lennard-Jones parameters
* `physics::water_tip3p_lr<>`: Water molecule, using a flexible TIP3P model optimized for long-range interactions
* `physics::water_fba_eps<>`: Water molecule, using the FBA/&epsilon; model

To set the coordinates of the molecule:
```c++
    my_system.add_molecule(physics::water_tip3p<>, {1, 2, 3});
```
where the coordinates are given in angstrom.
To advance the system by one step, do:
```c++
    my_system.step();
```
It is possible to change the floating point type and the numerical integrator to be used for the simulation:
```c++
    physics::molecular_system<long double, physics::pefrl> my_system;
```
By default, the floating point type is `double` (64-bit floating point) and the integrator is `physics::leapfrog`. Some currently available numerical integrators are:

* `physics::symplectic_euler`: Symplectic Euler method (1st order, 1 stage)
* `physics::leapfrog`: Leapfrog method (2nd order, 1 stage)
* `physics::multi_timestep_leapfrog`: Multi-timestep leapfrog method (2nd order, 1 stage for long-range forces)
* `physics::stochastic_leapfrog`: Stochastic "leapfrog" method (1 stage)
* `physics::damped_leapfrog`: Damped "leapfrog" method (1 stage)
* `physics::isokinetic_leapfrog`: Isokinetic "leapfrog" method (2nd order, 1 stage)
* `physics::nose_hoover`: Nosé-Hoover thermostats chain integrator (2nd order, 1 stage): it approximates a canonical (NVT) ensemble
* `physics::pefrl`: Position-extended Forest-Ruth-like method (4th order, 4 stages)
* `physics::vefrl`: Velocity-extended Forest-Ruth-like method (4th order, 4 stages)
* Composition schemes (they are structure-preserving, and can be used to construct higher-order, also symplectic, methods starting from 2nd order ones):
  * `physics::forest_ruth`: Forest-Ruth method (4th order, 3 stages)
  * `physics::suzuki4`: Suzuki method (4th order, 5 stages)
  * `physics::yoshida6`: Yoshida method (6th order, 8 stages)
  * `physics::suzuki8`: Suzuki method (8th order, 15 stages)
  * `physics::kahan_li10a`: Kahan-Li method (10th order, 31 stages)
* Explicit Runge-Kutta schemes:
  * `physics::ralston2`: Ralston method (2th order, 2 stages)
  * `physics::ralston4`: Ralston method (4th order, 4 stages)
  * `physics::butcher6`: Butcher method (6th order, 7 stages)
  * `physics::verner8`: Verner method (8th order, 11 stages)

It is possible to set a custom time step (in picoseconds) by adding a parameter to the `step` method:
```c++
    my_system.step(5e-4);
```
The biggest value for the time step so that leapfrog integration is stable is `2e-3` (2 femtoseconds, 1 is the default). This value corresponds more or less to the vibration period of O–H bonds.

To create a window, simply do:
```c++
#include "graphics.hpp"

int main()
{
    graphics my_window;

    while (!my_window.should_close())
    {}
}
```
To draw a frame of our system inside it, call the method:
```c++
    my_window.draw(my_system);
```
If we want to simulate the system until the window is closed:
```c++
    while (!my_window.should_close())
    {
        my_system.step();
        my_window.draw(my_system);
    }
```
