# M2Internship
Julia simulations of N 2-level atoms on a lattice interacting via the electric and magnetic interaction.

A lot of different system can be simulated with this method by changing the Hamiltonians and parameters of your system.

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

### Installation
To run the simulations, specific version of _QuantumCumulants.jl_ and _CollectiveSpins.jl_ are needed. Please run the Setup.jl file to install the correct dependencies. It will create a virtual environment (change the working directory in the Setup file if you are not already inside). The _CollectiveSpins.jl_ is deprecated, so the modified _CollectiveSpins.jl_ package on this Github has to be used. Download and unzip it, then go in you julia/dev directory and change the package by the one of this github.

In order to compute the time evolution of a large number of atoms N, meanfield plus correlations (MPC) simulations can be performed. Indeed, the full resolution of the master equation is quickly impossible, as all the operators scale as $2^N$.

The derivation of the (MPC) is done by the _QuantumCumulants.jl_ library, which derives the symbolic equations for the desired operators thanks to a cumulant expansion method, and automatically closes the set of differential equations.

In order to work with a symbolic number of atoms, we need to use the IndexedOperators objects, which are only implemented for transition and creation/annihilation operators. We thus first rewrite the Hamiltonians with these operators, then derive our set of differential equations, and finally expand the sums with our real number of atoms and convert the parameters into their numerical value. We then save the functions in C files, where each file correspond to one differential equation.

The equations are converted in C for 2 reasons:

- C functions can be compiled in parallel thanks to the Make file
- When a new Julia session is created, there is no trivial way of loading compiled functions, so the compilation has to be done again. With the C functions, we do not have this issue.

Once compiled, the C functions are linked to the dispatcher C function, which will call all the subfunctions (or differential equations) in a single call from Julia. Indeed, calling each subfunction can slow down the simulations and critically fill the memory.

### Run

After running the Setup.jl file to load your virtual environment, run the ElecMgtDD_QC_CFunctions_Op.ipynb file. It will create a Cfunctions directory, where all the functions will be stored, a dispatcher that will call all the subfunction in a single call from Julia, and an objs.txt file that helps the compilation of the Cfunction. Then, run in a command shell the command _make -f MakefileMac -jNbrCores_, where Nbrcores is the number of cores you want to use to compile your functions. If your on windows, run _make -f MakefileWindows -jNbrCores_. Finally, run the QC_solve.ipynb to solve the differential equations. If you are using windows, change the name of the library in the ccall with "liballfuncs.dll".


