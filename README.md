# Molecular dynamics 

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
