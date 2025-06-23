using Pkg
Pkg.activate(".") # Change with working directory if you are not already inside
Pkg.add(PackageSpec(name="QuantumCumulants", version="0.3.8"))
Pkg.add(PackageSpec(name="JLD2", version="0.5.13"))
Pkg.add(PackageSpec(name="OrdinaryDiffEq", version="6.98.0"))
Pkg.add(PackageSpec(name="PyPlot", version="2.11.6"))
Pkg.add(PackageSpec(name="QuantumOptics", version="1.2.3"))
Pkg.add(PackageSpec(name="Statistics", version="1.11.1"))
Pkg.add(PackageSpec(name="Symbolics", version="6.31.1"))

# CollectiveSpins is outdated. Import it in dev mode, then go in you julia/dev directory and change the package by the one of this github (after unzipping it)
Pkg.develop("CollectiveSpins")