using Printf

# Entrée / sortie
input_file = "objs.txt"
output_dir = "subdispatchers_objs"
isdir(output_dir) || mkdir(output_dir)

# Lecture
lines = readlines(input_file)

# Dispatcher
dispatcher_line = filter(l -> occursin("dispatcher.o", l), lines)[1]

# Extraire, parser, trier
diffeq_lines = filter(l -> occursin("diffeqf_", l), lines)
diffeq_tuples = [(parse(Int, match(r"diffeqf_(\d+)", l).captures[1]), l) for l in diffeq_lines]
sorted_diffeq_lines = [x[2] for x in sort(diffeq_tuples, by = x -> x[1])]

# Chunking
chunk_size = 1000
num_chunks = ceil(Int, length(sorted_diffeq_lines) / chunk_size)

# Génération
for i in 1:num_chunks
    start_idx = (i - 1) * chunk_size + 1
    end_idx = min(i * chunk_size, length(sorted_diffeq_lines))
    chunk = sorted_diffeq_lines[start_idx:end_idx]

    file_path = joinpath(output_dir, @sprintf("subdispatcher_%03d.txt", i))

    # Remplacer dispatcher.o par subdispatcher/subdispatcher_XXX.o
    new_dispatcher_line = replace(dispatcher_line, "dispatcher.o" => @sprintf("subdispatchers/subdispatcher_%03d.o", i))

    open(file_path, "w") do io
        println(io, new_dispatcher_line)
        for line in chunk
            println(io, line)
        end
    end
end

println("✅ Fichiers créés dans le dossier '$output_dir' ✅")
