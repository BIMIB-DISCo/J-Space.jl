#function that compute num of mutation for all edges
function evolution_seq(Tree::AbstractMetaGraph, Neutral_mut_rate::AbstractFloat,
                                                                length_ROI::Int)
    mutations = []
    #println("evolution_seq funzione")
    for e in edges(Tree)
        len = get_prop(Tree, dst(e), :Time)
    #    println("len: ",len)
        λ = Neutral_mut_rate * length_ROI * len
    #    println("valore λ: ", λ)
        n_mut = rand(Poisson(λ), 1)[1]
    #    println("num mut dalla distribuzione di poisson: ",n_mut)
        push!(mutations,n_mut)
        set_prop!(Tree, src(e), dst(e), :N_Mutations, n_mut)
    end
    return mutations
end

#change genome in position on pos
function transform_genome(genome::LongDNASeq, pos::Vector{Any})
    transition_matrix = DataFrame(A = [0.0, 0.33333333, 0.33333333, 0.33333333],
                                  C = [0.33333333, 0.0, 0.33333333, 0.33333333],
                                  G = [0.33333333, 0.33333333, 0.0, 0.33333333],
                                  T = [0.33333333, 0.33333333, 0.33333333, 0.0])
    for p in pos
        #println("pos: ",p)
        nucleotide = genome[p]
        #println("subsequence: ",genome[p-2:p+2])
        #println("nuclotide: ",nucleotide)
        #come nella simulatione, ho una certa prob che avvenga una simulatione
        prob_cum = cumsum(transition_matrix[!, string(nucleotide)])
        #println("prob_cum: ",prob_cum)
        k = round( rand(), digits=7) # with 7 decimals
        target_nucl = collect(k .<= prob_cum)
        #println("target_nucl: ",target_nucl)
        min = findfirst(target_nucl)
        #println("min: ",min)
        new_nucleotide = collect(names(transition_matrix)[min])[1]#type char
        #println("new_nucleotide: ",new_nucleotide)
        genome[p] = DNA(new_nucleotide)
    end
    return genome
end

#check that num of mutation for each path (root -> leaf) is less than len_ROI
function check_num_path_mutation(Tree::AbstractMetaGraph, len_ROI::Int)
    leafs = get_leafs(Tree)
    check = false
    i = 1
    while check == false && i <= length(leafs)
        yen_k = yen_k_shortest_paths(Tree, 1, leafs[i])
        path = yen_k.paths[1]
        sum = 0
        for v in 1:length(path)-1
            sum = sum + get_prop(Tree, path[v], path[v+1], :N_Mutations)
        end
        if sum > len_ROI
            check = true
        end
        i += 1
    end
end


#create a input for tool ART -> FASTA file and tree on format newick
function SC_experiment(Tree::AbstractMetaGraph, neural_mut_rate::Float64;
                                         path::String = "", len_ROI::Int = 6000,
                                         single_cell::Bool=true)

    #crete/load reference genome
    if path != ""
        open(FASTA.Reader, path) do reader
            for record in reader
                g_seq = FASTX.sequence(record)
            end
        end
        len_ROI = length(g_seq)
    else
        g_seq = randdnaseq(len_ROI)
        rec = FASTA.Record("Reference", g_seq)
        w = FASTA.Writer(open("my-out.fasta", "w"))
        write(w, rec)
        close(w)
    end

    #compute mutations ∀ node
    n_mutations = evolution_seq(Tree, neural_mut_rate, len_ROI)
    check = check_num_path_mutation(Tree, len_ROI)
    if check == true #qui in realtà è il path che pesa di più, con più mutazioni
        return "ERROR: mutations are more than lenght ROI -> (file fasta)"
    end
    possible_position = Set(1:len_ROI)
    #g_seq_e = copy(g_seq) #reference g_seq not change
    position_used = []
    #create fasta ∀ nodes
    for e in edges(Tree)
        g_seq_e = LongDNASeq()
        if has_prop(Tree, src(e), :Fasta)
            #println("il source ha un file Fasta")
            g_seq_e = copy(get_prop(Tree, src(e), :Fasta))
        else
            g_seq_e = copy(g_seq)
        end
        n_mut = get_prop(Tree, e, :N_Mutations)
        pos_edge = []
        for i=1:n_mut
            pos = rand(possible_position)
            delete!(possible_position, pos)
            push!(pos_edge, pos)
            push!(position_used, pos)
        end
        g_seq_e = transform_genome(g_seq_e, pos_edge)
        set_prop!(Tree, dst(e), :Fasta, g_seq_e)
    end

    #return fasta of leaf nodes
    leafs = get_leafs(Tree)
    fasta_samples = []
    for l in leafs
        f = get_prop(Tree, l, :Fasta)
        push!(fasta_samples, f)
    end

    #write fasta on files if single_cell is true
    if single_cell
        mkpath("Fasta output")#create folder
        for i in 1:length(leafs)
            w = FASTA.Writer(open("Fasta output\\sample" * string(leafs[i])
                                                               * ".fasta", "w"))
            rec = FASTA.Record("Sample" * string(leafs[i]), fasta_samples[i])
            write(w, rec)
            close(w)
        end
    end
    return g_seq, fasta_samples, position_used
end
