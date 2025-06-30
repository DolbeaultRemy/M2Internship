using QuantumCumulants
using CollectiveSpins
using Symbolics
using JLD2
using ProgressMeter
using Random

""" Test if the atom i and j are NN, returns 1 if yes, 0 else """
function NNTest(i, j)
    Matidx = reshape([1:N;], Nx, Ny, Nz)
    idxi = findall(x->x==i, Matidx)[1]
    x, y, z = idxi[1], idxi[2], idxi[3]

    idxj = findall(x->x==j, Matidx)[1]
    xj, yj, zj = idxj[1], idxj[2], idxj[3]

    if (abs(x - xj) == 1 && y == yj && z == zj) ||
       (abs(y - yj) == 1 && x == xj && z == zj) ||
       (abs(z - zj) == 1 && x == xj && y == yj)
        return 1
    else
        return 0
    end
end

function nbrNN(i::Int)
    Matidx = reshape([1:N;], Nx, Ny, Nz)
    idxi = findall(x->x==i, Matidx)[1]
    x, y, z = idxi[1], idxi[2], idxi[3]
    nbr = 0
    for j = 1:N
        idxj = findall(x->x==j, Matidx)[1]
        xj, yj, zj = idxj[1], idxj[2], idxj[3]
    
        if (abs(x - xj) == 1 && y == yj && z == zj) ||
           (abs(y - yj) == 1 && x == xj && z == zj) ||
           (abs(z - zj) == 1 && x == xj && y == yj)
            nbr +=1
        end
    end
    return nbr
end

""" Choose (1-filling_fraction)*N elements to remove in a list of N elements """
function idx_rm(N, filling_fraction)
    # Choose randomly the (1-filling_fraction*N) spins to remove
    idx_pos_rm = randperm(N)[1:round(Int, (1-filling_fraction)*N)]
    return idx_pos_rm
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
        for i in 1:length(eqs_eval)
            println(io, "Cfunctions/diffeqf_$i.o")
        end
    end
end

""" Function creating the makefile for the corresponding HT frequency """
function change_mkfile(var_filling_fraction, r)
    write("Makefile", replace(read("MakefileTemplate", String), "liballfuncs.dll"=>"libs/liballfuncs_FF_$(round(var_filling_fraction, digits=2))_rep_$r.dll"))
end

""" Function that create the eqs for a certain FF """
function CreateEquations(var_filling_fraction, eqs, r)
    # Choose atoms that will be removed
    irm = idx_rm(N, var_filling_fraction)

    # Set to 0 all the parameters with spin removed from geometry:
    Γij_ = [(i ∈ irm || j ∈ irm) ? 0 : Γ_CS[i, j] for i = 1:N, j = 1:N]
    Ωij_ = [(i ∈ irm || j ∈ irm) ? 0 : Ω_CS[i, j] for i = 1:N, j = 1:N if i≠j]
    
    a_ = [(i ∈ irm || j ∈ irm) ? 0 : Ω_ZZ_coeff[1]*NNTest(i, j) for i = 1:N for j=1:N if i≠j]
    b_ = [(i ∈ irm || j ∈ irm) ? 0 : Ω_ZZ_coeff[2]*NNTest(i, j) for i = 1:N for j=1:N if i≠j]
    c_ = [(i ∈ irm || j ∈ irm) ? 0 : Ω_ZZ_coeff[3]*NNTest(i, j) for i = 1:N for j=1:N if i≠j]
    nn_ = [(i ∈ irm) ? 0 : nbrNN(i) for i = 1:N]

    # Evaluate equations
    p0 = [Γij_...; Ωij_...; a_; b_; c_; nn_; Ω_ZZ_coeff[2]; Ω_ZZ_coeff[3]]
    p0 = ps .=> p0;
    eqs_eval = substitute(eqs, Dict(p0));

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
    @save "Operators/op_list_FF_$(round(var_filling_fraction, digits=2))_rep_$r.jdl2" op_list

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
Nx,Ny,Nz = [5,5,1]
N = Nx*Ny*Nz

println("N = "*string(N))
d_xy, d_z = 266., 532. # Optical lattice spacing in nm
λ = 1299.
Γ0 = 1. # In Hz
e = [0,0,1]   # Quantization axis

# Magnetic interaction coeffs
Ω_ZZ_coeff = [53, 42, 33]./Γ0

a_dim,b_dim,c_dim = [d_xy,d_xy,d_z]/λ
geo = geometry.box(a_dim,b_dim,c_dim;Nx=Nx,Ny=Ny,Nz=Nz)
system = CollectiveSpins.SpinCollection(geo, e, gammas=1.)

Ω_CS = interaction.OmegaMatrix(system)
Γ_CS = interaction.GammaMatrix(system)

repetition = 1; # Nbr of repetition per FF

@cnumbers Nsymb bn cn

h = NLevelSpace(Symbol(:atom),2)

Γ(i,j) = IndexedVariable(:Γ,i,j)
Ω(i,j) = IndexedVariable(:Ω,i,j;identical=false)
a(i,j) = IndexedVariable(:a,i,j;identical=false)
b(i,j) = IndexedVariable(:b,i,j;identical=false)
c(i,j) = IndexedVariable(:c,i,j;identical=false)
nn(i) = IndexedVariable(:n, i) # Nbr of NN of atom i
Ωht(i) = IndexedVariable(:Ω_ht,i)

i = Index(h,:i,Nsymb,h)
j = Index(h,:j,Nsymb,h)
k = Index(h,:k,Nsymb,h)
l = Index(h,:l,Nsymb,h)

σ(x,y,z) = IndexedOperator(Transition(h,:σ,x,y),z)

H_elec = Σ(Ω(i,j)*σ(2,1,i)*σ(1,2,j), j, i)
H_mgt1 = Σ(a(i,j)*σ(2,2,i)*σ(2,2,j)
        - 2*b(i,j)*σ(2,2,i)*σ(2,2,j)
        + c(i,j)*σ(2,2,i)*σ(2,2,j), j, i)
        
H_mgt2 = 2*Σ(nn(i)*(bn*σ(2,2,i)-cn*σ(2,2,i)), i)

H = Symbolics.simplify(H_elec+H_mgt1+H_mgt2);

J = [σ(1,2,i)] # σ-
rates = [Γ(i,j)]

ops = [σ(2, 2, k), σ(2, 1, k)] # n_up/σ+

eqs = meanfield(ops,H,J;rates=rates,order=2)
complete!(eqs);

eqs_eval = evaluate(eqs; limits=(Nsymb=>N));

Γij_symb = [Γ(i,j) for i = 1:N for j=1:N]
Ωij_symb = [Ω(i,j) for i = 1:N for j=1:N if i≠j]

a_symb = [a(i,j) for i = 1:N for j=1:N if i≠j]
b_symb = [b(i,j) for i = 1:N for j=1:N if i≠j]
c_symb = [c(i,j) for i = 1:N for j=1:N if i≠j]

nn_symb = [nn(i) for i = 1:N]

ps = [Γij_symb; Ωij_symb; a_symb; b_symb; c_symb; nn_symb; bn; cn];

for i = 1:1:1
    println(string(i)*"/"*string(N))
    var_filling_fraction = i/N
    for r in 1:repetition
        println(string(r)*"/"*string(repetition))
        println("Create equations")
        t = @elapsed begin
            CreateEquations(var_filling_fraction, eqs_eval, r)
        end
        println("Time taken $i : $t secondes")
        println("Compile Cfunctions")
        t = @elapsed begin
            change_mkfile(var_filling_fraction, r)
            run(`make -j11`)
        end
        println("Time taken $i : $t secondes")
        println()
    end
end