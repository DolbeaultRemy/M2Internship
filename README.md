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

where  [a, b, c] = [53, 42, 33] Hz.

## Running the simulations

In order to compute the time evolution of a large number of atoms N, meanfield plus correlations (MPC) simulations can be performed. Indeed, the full resolution of the master equation is quickly impossible, as all the operators scale as $2^N$.

The derivation of the (MPC) is done by the _QuantumCumulants.jl_ library, which derives the symbolic equations for the desired operators thanks to a cumulant expansion method, and automatically closes the set of differential equations.

In order to work with a symbolic number of atoms, we need to use the IndexedOperators objects, which are only implemented for transition and creation/annihilation operators. We thus first rewrite the Hamiltonians with these operators, then derive our set of differential equations, and finally expand the sums with our real number of atoms and convert the parameters into their numerical value. We then save the functions in C files, where each file correspond to one differential equation.

The equations are converted in C for 2 reasons:

- C functions can be compiled in parallel thanks to the Make file
- When a new Julia session is created, there is no trivial way of loading compiled functions, so the compilation has to be done again. With the C functions, we do not have this issue.

Once compiled, the C functions are linked to the dispatcher C function, which will call all the subfunctions (or differential equations) in a single call from Julia. Indeed, calling each subfunction can slow down the simulations and critically fill the memory.


