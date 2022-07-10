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
<div align="center">
<img src="https://latex.codecogs.com/svg.image?V=\sum_{i\sim%20j}K_{ij}^r\left(r_{ij}-r_{ij}^0\right)^2%2B\sum_{i\sim%20j\sim%20k}K_{ijk}^{\theta}%20\left(\theta_{ijk}-\theta_{ijk}^0\right)^2%2B\sum_{i\sim\cdot\sim%20k}K_{ik}^{UB}\left(r_{ik}-r_{ik}^0\right)^2\\"/>
</div>
<div align="center">
<img src="https://latex.codecogs.com/svg.image?\qquad%2B\sum_{i\sim%20j\sim%20k\sim%20l%2Cn}K_{ijkl}^{\chi}%20\left(1%2B\cos\left(n\chi_{ijkl}-\delta_{ijkl}\right)\right)%2B\sum_{ijkl\%20\text{impropers}}K_{ijkl}^{\psi}\left(\psi_{ijkl}-\psi_{ijkl}^0\right)^2\\"/>
</div>
<div align="center">
<img src="https://latex.codecogs.com/svg.image?\qquad%2B\sum_{i<j}\epsilon_{ij}\left(\left(\frac{R_{ij}}{r_{ij}}\right)^{12}-2\left(\frac{R_{ij}}{r_{ij}}\right)^6\right)%2B\sum_{i<j}k_C\frac{q_i%20q_j}{r_{ij}^2}"/>
</div>

The first term describes the bond potential, modeled as an elastic potential between bonded atoms (with elastic constant <img src="https://latex.codecogs.com/svg.image?k^r=2K^r"/>), the second term is the angle potential, the third term is the Urey-Bradley potential, acting between 1-3 atoms, the fourth and fifth terms describe the dihedral (proper and improper) angles potentials and the last two terms correspond to non-bonded potentials (Lennard-Jones and electrostatic respectively).

References:
* M. P. Allen, D. J. Tildesley, *Computer Simulation of Liquids*, Oxford University Press, 2017
* P. Spijker, B. Markvoort, P. Hilbers, *Parallel Utility for Modeling of Molecular Aggregation*, Biomodeling and bioinformatics, Eindhoven University of Technology, 2007

#### Table of contents

* [Molecular dynamics](#molecular-dynamics)
  * [Ewald summation](#ewald-summation)
  * [PPPM method](#pppm-method)
  * [Nosé-Hoover thermostat](#nosé-hoover-thermostat)
  * [Integration schemes](#integration-schemes)
* [The code](#the-code)
  * [Dependencies](#dependencies)
  * [Compilation](#compilation)
  * [Basic usage](#basic-usage)

![Animation of a system of water molecules](animation.gif)

### Ewald summation

In a cubic periodic system of side <img src="https://latex.codecogs.com/svg.image?L"/>, the electrostatic potential is given by:

<div align="center">
<img src="https://latex.codecogs.com/svg.image?V_C=\frac{1}{2}\sum_{ij}^N\sum_{\mathrm{n}\in%20Z^3}^{%27}\frac{z_i%20z_j}{\left|\mathrm{r}_{ij}%2B\mathrm{n}L\right|}"/>
</div>

where the prime symbol means that the <img src="https://latex.codecogs.com/svg.image?i=j"/> term must be excluded for <img src="https://latex.codecogs.com/svg.image?\mathrm{n}=\mathrm{0}"/>, and <img src="https://latex.codecogs.com/svg.image?z_i=\sqrt{k_C}q_i"/>. As it is not practical to calculate all the contributions directly, it is more convenient to calculate the long-range part in Fourier space, leading to a formula which converges much faster than the previous equation. To do so, the Green's function must be separated in the following way:

<div align="center">
<img src="https://latex.codecogs.com/svg.image?\frac{1}{r}=\frac{f(r)}{r}%2B\frac{1-f(r)}{r}"/>
</div>

where <img src="https://latex.codecogs.com/svg.image?f(r)"/> is sometimes referred to as the splitting function. Choosing <img src="https://latex.codecogs.com/svg.image?f(r)"/> to be the complementary error function, one obtains the classical Ewald summation, whose terms are given by:

<div align="center">
<img src="https://latex.codecogs.com/svg.image?V_C=V^{(r)}%2BV^{(k)}%2BV^{(s)}%2BV^{(d)}"/>
</div>
<div align="center">
<img src="https://latex.codecogs.com/svg.image?V^{(r)}=\frac{1}{2}\sum_{ij}^N\sum_{\mathrm{n}\in%20Z^3}^{%27}z_i%20z_j\frac{\text{erfc}\left(\kappa\left|\mathrm{r}_{ij}%2B\mathrm{n}L\right|\right)}{\left|\mathrm{r}_{ij}%2B\mathrm{n}L\right|}"/>
</div>
<div align="center">
<img src="https://latex.codecogs.com/svg.image?V^{(k)}=\frac{1}{2L^3}\sum_{\mathrm{n}\neq\mathrm{0}}\frac{4\pi}{k_{\mathrm{n}}^2}e^{-k_{\mathrm{n}}^2/4\kappa^2}\left|\tilde{\rho}(\mathrm{k_n})\right|^2=\frac{1}{2L^3}\sum_{\mathrm{n}\neq\mathrm{0}}\tilde{g}(k_{\mathrm{n}})\tilde{\gamma}(k_{\mathrm{n}})\left|\tilde{\rho}(\mathrm{k_n})\right|^2"/>
</div>
<div align="center">
<img src="https://latex.codecogs.com/svg.image?V^{(s)}=-\frac{\kappa}{\sqrt{\pi}}\sum_i%20z_i^2"/>
</div>
<div align="center">
<img src="https://latex.codecogs.com/svg.image?V^{(d)}=\frac{2\pi}{(1+2\epsilon_r)L^3}\left(\sum_i%20z_i\mathrm{r}_i\right)^2"/>
</div>

where <img src="https://latex.codecogs.com/svg.image?V^{(r)}"/> is the contribution from real space (short-range), <img src="https://latex.codecogs.com/svg.image?V^{(k)}"/> is the contribution from the reciprocal space (long-range), <img src="https://latex.codecogs.com/svg.image?V^{(s)}"/> is the self-energy correction and <img src="https://latex.codecogs.com/svg.image?V^{(d)}"/> is the dipole correction. The Fourier transform of the charge density <img src="https://latex.codecogs.com/svg.image?\tilde{\rho}"/> is defined as:

<div align="center">
<img src="https://latex.codecogs.com/svg.image?\tilde{\rho}(\mathrm{k})=\int_V%20\rho(\mathrm{r})e^{-i\mathrm{k}\cdot\mathrm{r}}d^3r=\sum_{j=1}^Nz_je^{-i\mathrm{k}\cdot\mathrm{r}_j"/>
</div>

The <img src="https://latex.codecogs.com/svg.image?\mathrm{k_n}"/>-vectors are given by <img src="https://latex.codecogs.com/svg.image?\mathrm{k_n}=2\pi\mathrm{n}/L"/>, while <img src="https://latex.codecogs.com/svg.image?\epsilon_r"/> is the relative dielectric constant (equal to 1 in vacuum) and <img src="https://latex.codecogs.com/svg.image?\kappa"/> is a free parameter, known as the Ewald parameter.

References:
* P. Ewald, *Die Berechnung optischer und elektrostatischer Gitterpotentiale*, Annalen der Physik, 369, pp. 253-287, 1921

### PPPM method

The particle-particle, particle-mesh (PPPM, or P<sup>3</sup>M) method can be applied to speed up the calculation of the reciprocal space term of Ewald summation thanks to fast Fourier transform (FFT) algorithms. Since FFT is based on discrete Fourier transforms (DFT), it requires sample points to be equally spaced, so a necessary preparatory step is to interpolate the charges on a 3-dimensional lattice (called mesh) with spacing <img src="https://latex.codecogs.com/svg.image?h"/>. The charge of a single mesh point <img src="https://latex.codecogs.com/svg.image?\mathrm{r_p}"/> is given by:

<div align="center">
<img src="https://latex.codecogs.com/svg.image?z_M(\mathrm{r_p})=\int_VW(\mathrm{r_p}-\mathrm{r})\rho(\mathrm{r})d^3r=\sum_{i=1}^Nz_iW(\mathrm{r_p}-\mathrm{r}_i)"/>
</div>

where <img src="https://latex.codecogs.com/svg.image?\mathrm{p}=h\mathrm{n}"/> and <img src="https://latex.codecogs.com/svg.image?W"/> is the charge assignment function, which is chosen so that the sum behaves as a convolution with a small window (so that the cost of the computation of <img src="https://latex.codecogs.com/svg.image?z_M"/> is <img src="https://latex.codecogs.com/svg.image?O(N)"/>). Its Fourier transform can be written as a DFT:

<div align="center">
<img src="https://latex.codecogs.com/svg.image?\tilde{z}_M(\mathrm{k_n})=\sum_{\mathrm{r_p}\in%20M}z_M(\mathrm{r_p})e^{-i\mathrm{k_n}\cdot\mathrm{r_p}"/>
</div>

The reciprical space term of the potential is then given by:

<div align="center">
<img src="https://latex.codecogs.com/svg.image?V^{(k)}=\frac{1}{2L^3}\sum_{\mathrm{n}\neq\mathrm{0}}\tilde{G}_{opt}(k_{\mathrm{n}})\left|\tilde{z}_M(\mathrm{k_n})\right|^2"/>
</div>

where <img src="https://latex.codecogs.com/svg.image?\tilde{G}_{opt}"/> is the optimal influence function:

<div align="center">
<img src="https://latex.codecogs.com/svg.image?\tilde{G}_{opt}(\mathrm{k})=\frac{\tilde{\mathrm{D}}(\mathrm{k})\cdot\sum_{\mathrm{m}\in%20Z^3}\tilde{U}^2\left(\mathrm{k}%2B\frac{2\pi}{h}\mathrm{m}\right)\tilde{\mathrm{R}}\left(\mathrm{k}%2B\frac{2\pi}{h}\mathrm{m}\right)}{\left|\tilde{\mathrm{D}}(\mathrm{k})\right|^2\left[\sum_{\mathrm{m}\in%20Z^3}\tilde{U}^2\left(\mathrm{k}%2B\frac{2\pi}{h}\mathrm{m}\right)\right]^2}"/>
</div>

where <img src="https://latex.codecogs.com/svg.image?\tilde{U}(\mathrm{k})=\tilde{W}(\mathrm{k})/h^3"/>, while <img src="https://latex.codecogs.com/svg.image?\tilde{\mathrm{D}}(\mathrm{k})"/> is the Fourier transform of the differential operator and <img src="https://latex.codecogs.com/svg.image?\tilde{\mathrm{R}}(\mathrm{k})"/> is the Fourier transform of the true reference force:

<div align="center">
<img src="https://latex.codecogs.com/svg.image?\tilde{\mathrm{R}}(\mathrm{k})=-i\mathrm{k}\tilde{g}(k)\tilde{\gamma}(k)"/>
</div>

Note that the influence function does not depend on the particles positions and thus can be calculated just once at the beginning of the simulation, as long as the volume of the simulation box does not vary. The series are highly convergent and can usually be truncated at <img src="https://latex.codecogs.com/svg.image?|\mathrm{m}|=2"/>. Note also that the same steps can be performed for the calculation of the dispersion forces (i.e. the ones arising from
<img src="https://latex.codecogs.com/svg.image?1/r^6"/> part of the Lennard-Jones potential or the Buckingham potential).

To obtain the forces, two schemes can be used depending on the application: in the *ik*-differentiation scheme the gradient is computed analytically in Fourier space, while in the *ad*-differentiation scheme the gradient is calculated in real space. The first method conserves the total momentum but not the energy and it is preferred in canonical ensembles because it is more accurate, while the second method conserves the energy but not the momentum, which is preferred in microcanonical ensembles (Ballenegger et al., 2012).

References:
* R. W. Hockney, J. W. Eastwood, *Computer Simulation Using Particles*, Bristol: Adam Hilger, 1988
* M. Deserno, C. Holm, *How to mesh up Ewald sums (I): A theoretical and numerical comparison of various particle mesh routines*, Journal of Chemical Physics, 109, 7678-7693, 1998
* V. Ballenegger, J. J. Cerdà, C. Holm, *How to Convert SPME to P3M: Influence Functions and Error Estimates*, Journal of Chemical Theory and Computation, 2012
* R. E. Isele-Holder, W. Mitchell, A. E. Ismail, *Development and application of a particle-particle particle-mesh Ewald method for dispersion interactions*, Journal of Chemical Physics, 137, 174107, 2012

### Nosé-Hoover thermostat

Simulating a system by integrating the Hamilton's equations of motion will naturally result in a microcanonical (NVE) ensemble, in which the total energy is conserved. In many applications, the internal energy of the system is not known a priori and it is more useful to control the temperature, and, in particular, one is interested in simulating a canonical (NVT) ensemble. One way to control the temperature is to employ isokinetic equations of motions, which are derived to constrain the kinetic energy to be conserved through a friction coefficient <img src="https://latex.codecogs.com/svg.image?\xi"/>:

<div align="center">
<img src="https://latex.codecogs.com/svg.image?\dot{\mathrm{r}}_i=\frac{\mathrm{p}_i}{m_i}%2C\qquad\dot{\mathrm{p}}_i=\mathrm{f}_i-\xi\mathrm{p}_i\\"/>
</div>
<div align="center">
<img src="https://latex.codecogs.com/svg.image?\xi=\frac{\sum_i\mathrm{p}_i\cdot\mathrm{f}_i/m_i}{\sum_ip_i^2/m_i}"/>
</div>

In this way, the temperature can be controlled by rescaling the kinetic energy accordingly at the beggining of the simulation. Using this method, configurations sample the canonical ensemble, but not the momenta (not following the Maxwell-Boltzmann distribution).

The Nosé-Hoover thermostat equations, instead, are (without scaling, given by Hoover, 1985):

<div align="center">
<img src="https://latex.codecogs.com/svg.image?\dot{\mathrm{r}}_i=\frac{\mathrm{p}_i}{m_i}%2C\qquad\dot{\mathrm{p}}_i=\mathrm{f}_i-\frac{p_{\eta}}{Q}\mathrm{p}_i"/>
</div>
<div align="center">
<img src="https://latex.codecogs.com/svg.image?\dot{\eta}=\frac{p_{\eta}}{Q}%2C\qquad\dot{p}_{\eta}=\sum_i\frac{p_i^2}{m_i}-gk_BT"/>
</div>

where <img src="https://latex.codecogs.com/svg.image?g"/> is the number of degrees of freedom while <img src="https://latex.codecogs.com/svg.image?Q"/> is a thermal intertia. These equations allow small oscillations in instantaneous temperature (whose magnitude depends on the number of the degrees of freedom), but, at equilibrium, they are able to sample the canonical ensemble for both configurations and momenta. They conserve an energy-like quantity, corresponding to the total Hamiltonian of the system containing the particles and the thermostat:

<div align="center">
<img src="https://latex.codecogs.com/svg.image?\mathcal{H}_{tot}=\mathcal{T}%2B\mathcal{V}%2B\frac{p_{\eta}^2}{2Q}%2Bgk_BT\eta"/>
</div>

For systems with few degrees of freedom, the system shows lack of ergodicity and the Nosé-Hoover thermostat shows problems in sampling the equilibrium distribution. It is possible, in this case, to use a thermostats chain:

<div align="center">
<img src="https://latex.codecogs.com/svg.image?\dot{\mathrm{r}}_i=\frac{\mathrm{p}_i}{m_i}%2C\qquad\dot{\mathrm{p}}_i=\mathrm{f}_i-\frac{p_{\eta_1}}{Q_1}\mathrm{p}_i"/>
</div>
<div align="center">
<img src="https://latex.codecogs.com/svg.image?\dot{\eta}_j=\frac{p_{\eta_j}}{Q_j}%2C\qquad\dot{p}_{\eta_j}=G_j-\frac{p_{\eta_{j+1}}}{Q_{j+1}}p_{\eta_j}%2C\qquad\dot{p}_{\eta_M}=G_M"/>
</div>

where the driving forces are defined as:

<div align="center">
<img src="https://latex.codecogs.com/svg.image?G_1=\sum_i\frac{p_i^2}{m_i}-gk_BT"/>
</div>
<div align="center">
<img src="https://latex.codecogs.com/svg.image?G_j=\frac{p_{\eta_{j-1}}^2}{Q_{j-1}}-k_BT"/>
</div>

and the conserved quantity is then:

<div align="center">
<img src="https://latex.codecogs.com/svg.image?\mathcal{H}_{tot}=\mathcal{T}%2B\mathcal{V}%2B\sum_{j=1}^M\frac{p_{\eta_j}^2}{2Q_j}%2Bgk_BT\eta_1%2B\sum_{j=2}^Mk_BT\eta_j"/>
</div>

while <img src="https://latex.codecogs.com/svg.image?M"/> is number of thermostats. This effectively adds <img src="https://latex.codecogs.com/svg.image?2M"/> degrees of freedom to the system, which may be useful to achieve ergodicity.

References:
* S. Nosé, *A unified formulation of the constant temperature molecular-dynamics methods*, Journal of Chemical Physics, 81 (1): pp. 511-519, 1984
* W. G. Hoover, *Canonical dynamics: equilibrium phase-space distributions*, Physical Review A, 31, pp. 1695-1697, 1985
* G. J. Martyna, M. L. Klein, *Nosé-Hoover chains: The canonical ensemble via continuous dynamics*, Journal of Chemical Physics, 97, 2635, 1992
* I. Fukuda, K. Moritsugu, *Coupled Nosé-Hoover equations of motions without time scaling*, Journal of Physics A, 50, 015002, 2016

### Integration schemes

Differential equations can be integrated using various integration schemes. In particular, for a dynamical Hamiltonian system, symplectic integrators can be used. Leapfrog method is the most widely used integrator for molecular dynamics, which has the symplectic and time-reversibility properties. An accurate symplectic fourth-order method with 4 stages can be constructed (Omelyan et al., 2002). However, to simulate a canonical ensemble, Nosé-Hoover equations need to be integrated. Although the same numerical schemes can be used for this purpose, specialized second-order integrators might be more appropriate and efficient (Itoh et al., 2013). Higher-order methods can be constructed starting from a second-order one through structure-preserving composition schemes (Kahan and Li, 1997). Unfortunately, composition methods cannot be used to accelerate most molecular dynamics simulation since the time-step must be smaller than a certain threshold to maintain stability (however they can be employed to get more reliable simulations).

Since the main sources of instability in a molecular dynamics simulation are the bonded terms of the potential (due to having much higher frequency than non-bonded terms), a way to accelerate simulations is to use different time-steps for bonded terms and non-bonded ones, since the computation of long-range forces require a considerable amount of time in the evaluation of the force field.

References:
* I. P. Omelyan, I. M. Mryglod, R. Folk, *Optimized Forest-Ruth- and Suzuki-like algorithms for integration of motion in many-body systems*, Computer Physics Communications, 2002
* W. Kahan, R.-C. Li, *Composition constants for raising the orders of unconventional schemes for ordinary differential equations*, Mathematics of Computation 66(219), pp. 1089-1099, 1997
* E. Hairer, G. Wanner, C. Lubich, *Geometric Numerical Integration. Structure-Preserving Algorithms for Ordinary Differential Equations*, Springer Series in Computational Mathematics, 2006
* [Mathematics Source Library C & ASM, Runge-Kutta Methods](http://www.mymathlib.com/diffeq/runge-kutta/)
* S. G. Itoh, T. Morishita, H. Okumura, *Decomposition-order effects of time integrator on ensemble averages for the Nosé-Hoover thermostat*, The Journal of Chemical Physics, 139, 2013

## The code
The code makes use of C++ templates and concepts (thus it requires C++20) and is organised in many header files, that can be included from one or more compilation units. It is organised in the following way:
* `gui`: directory which contains the management of the GUI (Graphical User Interface).
  * `shaders`: directory which contains shaders for the rendering of impostors (for fast rendering of spheres), post-processing filters (fast approximate anti-aliasing, blue noise dithering) and text.
  * `controls.hpp`: includes a function to manage keyboard and mouse controls.
  * `Font.hpp`: class for drawing text on the window.
  * `graphics.hpp`: class that manages the graphics used in this specific program.
  * `shader.hpp`: some functions to load shaders from file.
  * `stb_image.h`: a single-header file, public domain library for loading images (taken from [here](https://github.com/nothings/stb/blob/master/stb_image.h)).
* `math`: directory which contains helper functions and FFT implementation.
  * `dft.hpp`: contains the `dft` class which implements a simple radix-2 FFT algorithm for both real and complex inputs, and also multidimensional variants.
  * `helper.hpp`: some math helper functions.
* `physics`: directory which contains part of code relevant to the resolution of the physical/numerical problem.
  * `ewald.hpp`: Ewald summation and PPPM method for fast calculation of long-range forces.
  * `integrator.hpp`: classes and concepts for integration of (hamiltonian and not) dynamical systems.
  * `molecule.hpp`: classes, methods and other utilites for managing molecular systems (for now, you have only water molecules, also the dihedral potential part must be added).
  * `physics.hpp`: mainly a header file that includes everything.
  * `point.hpp`: classes, aliases and data structures for vectors, matrices and tensors with some helper function.
* `utils`: directory which contains some utilities used in this project.
  * `parallel_sort.hpp`: multi-threaded parallel sorting algorithm based on `std::sort`, `std::inplace_merge` and `utils::thread_pool`.
  * `thread_pool.hpp`: a simple multi-queue thread pool based on `std::jthread`. It reuses threads for many tasks instead of creating and destroying them continuously and thus avoids some overhead.
* `main.cpp`: it just contains the main loop and some basic initialization of the system (it can be ignored or modified).

### Dependencies

The graphics part has some dependencies on public libraries that enable the usage of modern OpenGL (Open Graphics Library):
* [GLFW 3](https://www.glfw.org/) (Graphics Library Framework 3): an open source, multi-platform API for creating windows, contexts and managing input and events.
* [GLEW](http://glew.sourceforge.net/) (OpenGL Extension Wrangler): a cross-platform library that include core and extended OpenGL functionalities.
* [FreeType](https://freetype.org/): an OpenGL library to render fonts (used in `Font.hpp`).

These dependencies are required only for these header files:
* `graphics.hpp`
* `Font.hpp` (included in `graphics.hpp`)
* `shader.hpp` (included in `graphics.hpp`)
* `controls.hpp` (included in `graphics.hpp`)

To install all dependencies at once on Ubuntu, you can run the shell script `dependencies_ubuntu.sh`.

### Compilation

To compile the program, simply do (with MinGW):

    g++ main.cpp -o mold -std=c++20 -I <includes> -L <libs> -lopengl32 -lglu32 -lglew32.dll -lglfw3dll -lfreetype -Wall -Wextra -pedantic -Ofast

On Linux, the library names could be different:

    g++ main.cpp -o mold -std=c++20 -I <includes> -L <libs> -lGL -lGLU -lGLEW -lglfw -lfreetype -Wall -Wextra -pedantic -Ofast

where `<includes>` and `<libs>` are the paths for installed libraries header files and static library files (if required). The executable will be called `mold`. If GCC is used for compilation, version 10+ is required for full C++20 support. Running the graphical part of the program requires OpenGL 3.3+.

To compile on Ubuntu, you can also run the shell script `compile.sh`.

### Basic usage

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
