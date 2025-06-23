# M2Internship
Julia simulations of N 2-level atoms on a lattice interacting via the electric and magnetic interaction.

We simulate N interacting atoms of Erbium on an optical lattice thanks to the _QuantumOptics.jl_, _CollectiveSpins.jl_, _QuantumCumulants.jl_ libraries. The atoms of Erbium are prepared in the -7/-6 excited/ground state.

The evolution of the density matrix is computed thanks to the Lindblad master equation:

$$\partial_t \rho = -\frac{i}{\hbar} [H, \rho] +\mathcal{L}[\rho]$$

where the Lindbladian is:

$$\mathcal{L}[\rho] = \sum_{i, j} \Gamma_{i, j}(J_i \rho J_j^\dagger - \frac{1}{2} \{J_i^\dagger J_j, \rho\})$$


## Electric dipole-dipole interaction:

The Hamiltonian of this interaction can be written as:

$$H_{\text{ElecDD}} = \sum_{i, j, i \neq j} \Omega_{i, j}\sigma_i^+ \sigma_j^-$$

It realizes a spin-flip exchange between two atoms.

The matrices $\Omega_{i, j}$ and $\Gamma_{i, j}$ (which is the collective decay matrix of the Lindbladian) are computed thanks to the _CollectiveSpins.jl_ library.


## Magnetic dipole-dipole interaction

The effect of the magnetic dipole-dipole interaction is taken into account by adding the following Hamiltonian to the Lindblald master equation:

$$H_{\text{MgtDD}} = \sum_{<i, j>} \Omega_{i, j}^{mgt} \sigma_z^i \sigma_z^j$$

In the configuration of our experiment, this Hamiltonian can be rewritten as:

$$ H_{\text{MgtDD}} = \sum_{<i, j>} (a n^{\uparrow \uparrow}\_{ij} + b(n^{\uparrow \downarrow}\_{ij} + n^{\downarrow \downarrow}\_{ij}) + c n^{\downarrow \downarrow}\_{ij})$$

where  [a, b, c] = [53, 42, 33] Hz

## Running the simulations

In order to compute the evolution of a large numbre of atoms
