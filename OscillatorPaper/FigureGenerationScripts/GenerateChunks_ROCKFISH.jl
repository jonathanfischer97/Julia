using Combinatorics  # For the `combinations` function

function generate_triplet_chunks(num_chunks::Int)
    names = [:ka1, :kb1, :kcat1, :ka2, :kb2, :ka3, :kb3, :ka4, :kb4, :ka7, :kb7, :kcat7, :L, :K, :P, :A]

    triplets = collect(combinations(names, 3))
    chunk_size = div(length(triplets), num_chunks)
    remainder = length(triplets) % num_chunks
    
    chunks = []
    start_idx = 1
    for i in 1:num_chunks
        end_idx = start_idx + chunk_size - 1
        if remainder > 0
            end_idx += 1
            remainder -= 1
        end
        push!(chunks, (start_idx, end_idx))
        println("$start_idx:$end_idx")
        start_idx = end_idx + 1
    end
    return chunks
end

generate_triplet_chunks(parse(Int,ARGS[1]))