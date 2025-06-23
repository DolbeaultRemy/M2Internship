# M2Internship
Julia simulations of N 2-level atoms on a lattice interacting via the electric and magnetic interaction.

The evolution of the density matrix is computed thanks to the Lindblad master equation:

$$\partial_t \rho = -\frac{i}{\hbar} [H, \rho] +\mathcal{L}[\rho]$$

where the Lindbladian is:

$$\mathcal{L}[\rho] = \sum_{i, j} \Gamma_{i, j}(J_i \rho J_j^\dagger - \frac{1}{2} \{J_i^\dagger J_j, \rho\})$$


## Electric dipole-dipole interaction:

The Hamiltonian of this interaction can be written as:

$$H_{\text{ElecDD}} = \sum_{i, j, i \neq j} \Omega_{i, j}\sigma_i^+ \sigma_j^-$$

It realizes a spin-flip exchange between two atoms.


The matrices $\Omega_{i, j}$ and $\Gamma_{i, j}$ (which is the collective decay matrix of the Lindbladian) are computed thanks to the \textit{CollectiveSpins.jl} library.


## Magnetic dipole-dipole interaction

In our regime

$$H_{\text{MgtDD}} = \sum_{<i, j>} \Omega_{i, j}^{mgt} \sigma_z^i \sigma_z^j$$
