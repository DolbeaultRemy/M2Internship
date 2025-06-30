using Printf

# Input and output setup
input_file = "dispatcher.c"
output_dir = "subdispatchers"

# Ensure output directory exists
isdir(output_dir) || mkdir(output_dir)

# Read the dispatcher file
lines = readlines(input_file)

# Extract only diffeqf_i declarations
pattern = r"diffeqf_(\d+)"
matches = [match(pattern, line) for line in lines if occursin("diffeqf_", line)]
function_ids = [parse(Int, m.captures[1]) for m in matches if m !== nothing]
sort!(function_ids)

# Constants
chunk_size = 2000
num_chunks = ceil(Int, length(function_ids) / chunk_size)

# Generate subdispatch files
for i in 1:num_chunks
    start_idx = (i - 1) * chunk_size + 1
    end_idx = min(i * chunk_size, length(function_ids))
    group = function_ids[start_idx:end_idx]

    file_path = joinpath(output_dir, @sprintf("subdispatcher_%03d.c", i))
    open(file_path, "w") do io
        println(io, "#include <complex.h>\n#include <math.h>\n")
        for fid in group
            println(io, @sprintf("extern void diffeqf_%d(double complex* du, const double complex* RHS1);", fid))
        end
        println(io, "")
        println(io, @sprintf("void call_diffeq_group_%03d(double complex* du, const double complex* RHS1) {", i))
        for fid in group
            println(io, @sprintf("    diffeqf_%d(du, RHS1);", fid))
        end
        println(io, "}")
    end
end

println("✅ Generated $(num_chunks) subdispatch files in '$output_dir'.")