using QuantumCumulants
using CollectiveSpins
using LinearAlgebra
using Symbolics
using JLD2
using PhysicalConstants.CODATA2018: μ_B, μ_0, h
using Unitful

""" Mgt DDI for general vectors """
function DDI_general(r, m1, m2)
    rnorm = norm(r)
    prefac = ustrip(μ_0) * ustrip(μ_B)^2 / (4 * π * rnorm^3 * (λ*1e-9)^3)
    return prefac * (dot(m1, m2) - 3 * dot(m1, r) * dot(m2, r) / rnorm^2) / ustrip(h)
end

""" Return the A, B, and C arrays for the computation of the HMgtDD Hamiltonian """
function MgtCoeffsArray(system, μ_es, μ_gs)
    A, B, C = zeros(N, N), zeros(N, N), zeros(N, N)
    for i = 1:N
        ri = system.spins[i].position
        for j = 1:N 
            if j ≠ i
                d =  system.spins[j].position .- system.spins[i].position
                A[i, j] = DDI_general(d, μ_es, μ_es)
                B[i, j] = DDI_general(d, μ_es, μ_gs)
                C[i, j] = DDI_general(d, μ_gs, μ_gs)
            end
        end
    end
    return A, B, C
end

""" Create the big function that will call all the subfunctions to avoid doing a lot of ccalls """
function generate_full_dispatcher(filename::String, n::Int)
    open(filename, "w") do io
        println(io, "#include <complex.h>\n#include <math.h>\n")
        # External subfunctions
        for i in 1:(n)
            println(io, "extern void diffeqf_$i(double complex* du, const double complex* RHS1);")
        end

        println(io, "\nvoid diffeqf(double complex* du, const double complex* RHS1) {")

        # Call the subfunctions
        for i in 1:(n)
            println(io, "    diffeqf_$i(du, RHS1);")
        end

        println(io, "}")
    end
end

""" Create the object file (with the name of all the functions and their corresponding files) to avoid compilation issue """
function objs_file()
    open("objs.txt", "w") do io
         println(io, "dispatcher.o")
        for i in 1:length(eqs_eval_sum)
            println(io, "Cfunctions/diffeqf_$i.o")
        end
    end
end

""" Function creating the makefile for the corresponding HT frequency """
function change_mkfile(thetaBext, phiBext)
    write("Makefile", replace(read("MakefileTemplate", String), "liballfuncs.dll"=>"libs/liballfuncs_thetaBext_$(round(thetaBext, digits=2))_phiBext_$(round(phiBext, digits=2)).dll"))
end

""" Function that create the eqs for a certain angle of the quantization axis """
function CreateEquations(thetaBext, phiBext, eqs)
    # Define the CollectiveSpins system
    e = [sin(thetaBext)*cos(phiBext), sin(thetaBext)*sin(phiBext), cos(thetaBext)]   # Quantization axis
    a_dim,b_dim,c_dim = [d_xy,d_xy,d_z]/λ
    geo = geometry.box(a_dim,b_dim,c_dim;Nx=Nx,Ny=Ny,Nz=Nz)
    system = CollectiveSpins.SpinCollection(geo, e, gammas=1.)
    Ω_CS = interaction.OmegaMatrix(system)
    Γ_CS = interaction.GammaMatrix(system)

    # Magnetic interaction
    gJ_gs = 1.163801 # Lande g-factors
    gJ_ex = 1.2599
    μ_gs = -6*gJ_gs .* e # Magnetic dipole moment of the GS
    μ_es = -7*gJ_ex .* e # Magnetic dipole moment of the ES


    # Recompute parameters
    Γij_ = [Γ_CS[i, j] for i = 1:N for j=1:N]
    Ωij_ = [Ω_CS[i, j] for i = 1:N for j=1:N if i≠j]
    Amgt, Bmgt, Cmgt = MgtCoeffsArray(system, μ_es, μ_gs)
    A_ = [Amgt[i, j] for i = 1:N for j=1:N if i≠j]
    B_ = [Bmgt[i, j] for i = 1:N for j=1:N if i≠j]
    C_ = [Cmgt[i, j] for i = 1:N for j=1:N if i≠j]
    Bs_, Cs_ = sum(Bmgt, dims=1)', sum(Cmgt, dims=1)'
    p0 = [Γij_; Ωij_; A_; B_; C_; Bs_; Cs_]
    p0 = ps .=> p0
    eqs_eval = substitute(eqs, Dict(p0))

     # Build dir
    isdir("Cfunctions") || mkdir("Cfunctions")
    isdir("Operators") || mkdir("Operators")
    isdir("libs") || mkdir("libs")

    # Build and save variables
    op_list = []
    var_array = []
    for i in 1:length(eqs_eval)
        var = eqs_eval[i].lhs
        push!(var_array, var)
        
        v_str = string(var)
        em = eachmatch(r"σ(\d+)", v_str)
        ind = [m.captures[1] for m in em]
        push!(op_list, [parse(Int, i) for i in ind])
    end
    @save "Operators/op_list_thetaBext_$(round(thetaBext, digits=2))_phiBext_$(round(phiBext, digits=2)).jdl2" op_list

    # Build and save C functions
    for i in 1:length(eqs_eval)
        # Save the C function
        code = Symbolics.build_function([eqs_eval[i].rhs], var_array, target=Symbolics.CTarget())
        c_code = replace(code, 
            "im" => "*I", "double* du" => "double complex* du", "const double* RHS1" => "const double complex* RHS1", "du[0]" => "du[$(i-1)]", "diffeqf" => "diffeqf_$i")
        open("Cfunctions/diffeqf_$i.c", "w") do io
            print(io, "#include <complex.h>\n") # Include complex package
            print(io, c_code)
        end
        # Free RAM memomry
        code = nothing
        c_code = nothing
    end

    # Generate the dispatcher
    generate_full_dispatcher("dispatcher.c", length(eqs_eval))

     # Generate the objs.txt
    objs_file()

end

# Define geometry of system
Nx,Ny,Nz = [4,4,1]
N = Nx*Ny*Nz

println("N = "*string(N))
d_xy, d_z = 266., 532. # Optical lattice spacing in nm
λ = 1299.
Γ0 = 1. # In Hz

theta_Bext_array = range(50*pi/180, 70*pi/180, 10)
phi_Bext_array = range(35*pi/180, 55*pi/180, 10)


# Symbolic eqations

@cnumbers Nsymb bn cn

hspace = NLevelSpace(Symbol(:atom),2)

Γ(i,j) = IndexedVariable(:Γ,i,j)
Ω(i,j) = IndexedVariable(:Ω,i,j;identical=false)
A(i,j) = IndexedVariable(:A,i,j;identical=false)
B(i,j) = IndexedVariable(:B,i,j;identical=false)
C(i,j) = IndexedVariable(:C,i,j;identical=false)
Bsum(i) = IndexedVariable(:Bs,i)
Csum(i) = IndexedVariable(:Cs,i)

i = Index(hspace,:i,Nsymb,hspace)
j = Index(hspace,:j,Nsymb,hspace)
k = Index(hspace,:k,Nsymb,hspace)
l = Index(hspace,:l,Nsymb,hspace)

σ(x,y,z) = IndexedOperator(Transition(hspace,:σ,x,y),z)

# Define the symbolic Hamiltonians
H_elec = Σ(Ω(i,j)*σ(2,1,i)*σ(1,2,j), j, i)
H_mgt1 = Σ(A(i,j)*σ(2,2,i)*σ(2,2,j)
        - 2*B(i,j)*(σ(2,2,i)*σ(2,2,j))
        + C(i,j)*(σ(2,2,i)*σ(2,2,j)), j, i)
        
H_mgt2 = 2*Σ((Bsum(i)*σ(2,2,i)-Csum(i)*σ(2,2,i)), i)

H = Symbolics.simplify(H_elec+H_mgt1+H_mgt2)

J = [σ(1,2,i)] # σ-
rates = [Γ(i,j)]

ops = [σ(2, 2, k), σ(2, 1, k)] # n_up/σ+

println("Calculate equations")
eqs = meanfield(ops,H,J;rates=rates,order=2)
println("Complete equations")
complete!(eqs)
println("Evaluate sums")
eqs_eval_sum = evaluate(eqs; limits=(Nsymb=>N))

# List of symbolic parameters
Γij_symb = [Γ(i,j) for i = 1:N for j=1:N]
Ωij_symb = [Ω(i,j) for i = 1:N for j=1:N if i≠j]
A_symb = [A(i,j) for i = 1:N for j=1:N if i≠j]
B_symb = [B(i,j) for i = 1:N for j=1:N if i≠j]
C_symb = [C(i,j) for i = 1:N for j=1:N if i≠j]
Bs_symb = [Bsum(i) for i = 1:N]
Cs_symb = [Csum(i) for i = 1:N]
ps = [Γij_symb; Ωij_symb; A_symb; B_symb; C_symb; Bs_symb; Cs_symb]


for i = 1:length(theta_Bext_array)
    println(string(i)*"/"*string(length(theta_Bext_array)))
    for j = 1:length(phi_Bext_array)
        println(string(j)*"/"*string(length(phi_Bext_array)))
        println("Create equations")
        t = @elapsed begin
            CreateEquations(theta_Bext_array[i], phi_Bext_array[j], eqs_eval_sum)
        end
        println("Time taken $i : $t secondes")
        println("Compile Cfunctions")
        t = @elapsed begin
            change_mkfile(theta_Bext_array[i], phi_Bext_array[j])
            run(`make -j11`)
        end
        println("Time taken $i : $t secondes")
        println()
    end
end