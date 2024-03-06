### -*- Mode: Julia -*-

### SingleCellExperiment.jl

using BioSequences
using FASTX
using LinearAlgebra


"""
Save fasta.
"""
function save_Fasta_W(Ref::LongSequence,
                      fasta_samples::Vector{Any},
                      tree::AbstractMetaGraph,
                      path_save_file::String)

      leafs = get_leafs(tree)
      mkpath(path_save_file * "\\Fasta output") # Create folder
      path_for_fasta = path_save_file * "\\Fasta output"
      for i in 1:length(leafs)
          w = FASTA.Writer(open(path_for_fasta *
                                "\\sample" *
                                string(leafs[i]) *
                                ".fasta",
                                "w"))
          rec = FASTA.Record("Sample" * string(leafs[i]), fasta_samples[i])
          write(w, rec)
          close(w)
      end
      w = FASTA.Writer(open(path_for_fasta *
                            "\\reference" *
                            ".fasta",
                            "w"))
      rec = FASTA.Record("Reference", Ref)
      write(w, rec)
      close(w)
end

function save_Fasta_L(Ref::LongSequence,
                      fasta_samples::Vector{Any},
                      tree::AbstractMetaGraph,
                      path_save_file::String)

      leafs = get_leafs(tree)
      mkpath(path_save_file * "/Fasta output") # Create folder
      path_for_fasta = path_save_file * "/Fasta output"
      for i in 1:length(leafs)
          w = FASTA.Writer(open(path_for_fasta *
                                "/sample" *
                                string(leafs[i]) *
                                ".fasta",
                                "w"))
          rec = FASTA.Record("Sample" * string(leafs[i]), fasta_samples[i])
          write(w, rec)
          close(w)
      end
      w = FASTA.Writer(open(path_for_fasta *
                            "/reference" *
                            ".fasta",
                            "w"))
      rec = FASTA.Record("Reference", Ref)
      write(w, rec)
      close(w)
end

"""
Computes number of mutations for all edges.
"""
function evolution_seq(Tree::AbstractMetaGraph,
                       Neutral_mut_rate::AbstractFloat,
                       length_ROI::Int,
                       max_time::AbstractFloat,
                       seed::MersenneTwister)

    mutations = []
    λ = Neutral_mut_rate * length_ROI
    for e in edges(Tree)
        Tfinal = get_prop(Tree, dst(e), :Time) - get_prop(Tree, src(e), :Time)
        t_curr = 0
        n_mut = 0
        #Tfinal = Tfinal/max_time
        while t_curr <= Tfinal
            t_curr = t_curr + rand(seed, Exponential(1 / λ), 1)[1]
            n_mut += 1
        end
        push!(mutations, n_mut)
        set_prop!(Tree, src(e), dst(e), :N_Mutations, n_mut)
    end
    return mutations
end


"""
Change genome in position pos.
"""
function transform_genome(genome::LongSequence,
                          pos,
                          seed::MersenneTwister,
                          mutations_tot::DataFrame,
                          sample::Int;
                          funct::Int = 0,
                          position_used::Vector{Any} = [],
                          muts_driver = [])

    transition_matrix =
        DataFrame(A = [0.0, 0.33333333, 0.33333333, 0.33333333],
                  C = [0.33333333, 0.0, 0.33333333, 0.33333333],
                  G = [0.33333333, 0.33333333, 0.0, 0.33333333],
                  T = [0.33333333, 0.33333333, 0.33333333, 0.0])
    substitution = []
    index = 1
    for p in pos
        nucleotide = genome[p]

        prob_cum = cumsum(transition_matrix[!, string(nucleotide)])

        k = round( rand(seed), digits = 7)
        target_nucl = collect(k .<= prob_cum)
        min = findfirst(target_nucl)

        new_nucleotide = collect(names(transition_matrix)[min])[1] # type char

        genome[p] = DNA(new_nucleotide)
        sub = string(p)*"_"*string(nucleotide)*">"*string(new_nucleotide)
        if funct == 1
            push!(mutations_tot, [p,
                                  string(nucleotide),
                                  string(new_nucleotide),
                                  sub,
                                  muts_driver[index],
                                  sample])
            push!(substitution, sub)
            push!(position_used, p)
            index += 1
        else
            push!(mutations_tot, [p,
                                  string(nucleotide),
                                  string(new_nucleotide),
                                  sub,
                                  false,
                                  sample])
        end
    end

    if funct == 1
        return genome, substitution, position_used, mutations_tot
    else
        return genome, mutations_tot
    end
end


"""
Checks that num of mutation for each path (root -> leaf) is less than
len_ROI
"""
function check_num_path_mutation(Tree::AbstractMetaGraph, len_ROI::Int)
    leafs = get_leafs(Tree)
    check = false
    i = 1
    while check == false && i <= length(leafs)
        yen_k = yen_k_shortest_paths(Tree, 1, leafs[i])
        path = yen_k.paths[1]
        sum = 0
        for v in 1:length(path)-1
            sum = sum + get_prop(Tree, path[v], path[v + 1], :N_Mutations)
        end
        if sum > len_ROI
            check = true
        end
        i += 1
    end
    return check
end


"""
Creates a input for tool ART -> FASTA file
"""
function Molecular_evolution_ISA(Tree::AbstractMetaGraph,
                             neural_mut_rate::Float64,
                             seed::MersenneTwister,
                             g_seq::LongSequence,
                             len_ROI::Int,
                             set_mut::Vector{Any})

    mutations_tot = DataFrame(Pos = Any[],
                              Reference = String[],
                              Alternative = String[],
                              Mut_id = Any[],
                              Driver = Any[],
                              Sample = Any[])

    max_time = max_time_nodes(Tree, get_leafs(Tree))
    ## Compute mutations ∀ node
    for i in [1,2,3]
        check = false
        n_mutations = evolution_seq(Tree,
                                    neural_mut_rate,
                                    len_ROI,
                                    max_time,
                                    seed)
        check = check_num_path_mutation(Tree, len_ROI)

        if check == true && i >= 3
            return "ERROR: mutations are more than length ROI -> (file fasta)."
        else
            break
        end
    end

        possible_position = Set(1:len_ROI)

        position_used = []

        # create fasta ∀ nodes
        for e in edges(Tree)

            g_seq_e = LongDNA{4}()
            #g_seq_e = LongDNASeq()

            if has_prop(Tree, src(e), :Fasta)
                g_seq_e = copy(get_prop(Tree, src(e), :Fasta))
            else
                g_seq_e = copy(g_seq)
            end

            n_mut = get_prop(Tree, e, :N_Mutations)
            pos_edge = []

            for i = 1:n_mut
                pos = rand(seed, possible_position)
                delete!(possible_position, pos)
                push!(pos_edge, pos)
                push!(position_used, pos)
            end

            subpop_father = get_prop(Tree, src(e), :Subpop_Child)
            subpop_child = get_prop(Tree, dst(e), :Subpop_Child)

            mut_f = set_mut[subpop_father]
            mut_child = set_mut[subpop_child]

            g_seq_e, mutations_tot = transform_genome(g_seq_e,
                                                      pos_edge,
                                                      seed,
                                                      mutations_tot,
                                                      dst(e))

            if length(mut_f) != length(mut_child)
                if typeof(mut_f) == String
                    num_mut_driver = length(mut_child) - 1
                else
                    num_mut_driver = length(mut_child) - length(mut_f)
                end
                new_muts = copy(mut_child)
                if typeof(mut_f) == String
                    setdiff!(new_muts, [mut_f])
                else
                    setdiff!(new_muts, mut_f)
                end

                #check not same mutation driver
                muts_already_write = []
                for nm in new_muts
                    if nm ∈ mutations_tot.Driver
                        push!(muts_already_write,nm)
                    end
                end
                setdiff!(new_muts, muts_already_write)
                num_mut_driver = length(new_muts)

                pos_edge_d = []
                for i = 1:length(new_muts)
                    pos = rand(seed, possible_position)
                    delete!(possible_position, pos)
                    push!(pos_edge_d, pos)
                end
                g_seq_e, sub, position_used, mutations_tot =
                               transform_genome(g_seq_e,
                                                pos_edge_d,
                                                seed,
                                                mutations_tot,
                                                dst(e),
                                                funct=1,
                                                position_used=position_used,
                                                muts_driver = new_muts)
            end

            set_prop!(Tree, dst(e), :Fasta, g_seq_e)
        end

        # Return fasta of leaf nodes
        leafs = get_leafs(Tree)
        fasta_samples = []

        for l in leafs
            f = get_prop(Tree, l, :Fasta)
            push!(fasta_samples, f)
        end

        return g_seq, fasta_samples, position_used, mutations_tot
end


"""
Call function ISA
"""
function experiment_ISA(tree::AbstractMetaGraph,
                        neural_mut_rate::Float64,
                        seed::MersenneTwister,
                        len_ROI::Int,
                        set_mut::Vector{Any};
                        frequency_dna::Vector{Float64} = [0.3,0.2,0.3])
    ## Create reference genome
    sw = SamplerWeighted(dna"ACGT", frequency_dna)
    g_seq = randseq(seed, DNAAlphabet{4}(), sw, len_ROI)
    # g_seq = randdnaseq(seed, len_ROI)
    rec = FASTA.Record("Reference", g_seq)
    w = FASTA.Writer(open("Reference.fasta", "w"))
    write(w, rec)
    close(w)
    
    g_seq, fasta_samples, position_used, mutations_tot =
        Molecular_evolution_ISA(tree,
                                neural_mut_rate,
                                seed,
                                g_seq,
                                len_ROI,
                                set_mut)

    mutations_tot_2 = copy(mutations_tot)
    leafs = get_leafs(tree)
    paths_tot = []
    for l in leafs
        yen_k = yen_k_shortest_paths(tree, 1, l)
        path = yen_k.paths[1]
        push!(paths_tot, path)
    end
    paths_tot
    
    for i in 1:length(mutations_tot_2.Sample)
        sample = mutations_tot_2.Sample[i]
        if sample ∉ leafs
            ls = leafs[sample .∈ paths_tot]
            mutations_tot_2.Sample[i] = ls
        end
    end
    return g_seq, fasta_samples, position_used, mutations_tot_2
end


function experiment_ISA(tree::AbstractMetaGraph,
                        neural_mut_rate::Float64,
                        seed::MersenneTwister,
                        path::String,
                        set_mut::Vector{Any};
                        frequency_dna::Vector{Any} = [])

    ## load reference genome
    g_seq = LongDNA{4}()
    ## g_seq_e = LongDNASeq()
    open(FASTA.Reader, path) do reader
        for record in reader
            g_seq = FASTX.sequence(record)
        end
    end

    len_ROI = length(g_seq)

    g_seq, fasta_samples, position_used, mutations_tot =
        Molecular_evolution_ISA(tree,
                                neural_mut_rate,
                                seed,
                                g_seq,
                                len_ROI,
                                set_mut)
    mutations_tot_2 = copy(mutations_tot)
    leafs = get_leafs(tree)
    paths_tot = []
    for l in leafs
        yen_k = yen_k_shortest_paths(tree, 1, l)
        path = yen_k.paths[1]
        push!(paths_tot, path)
    end
    paths_tot

    for i in 1:length(mutations_tot_2.Sample)
        sample = mutations_tot_2.Sample[i]
        if sample ∉ leafs
            ls = leafs[sample .∈ paths_tot]
            mutations_tot_2.Sample[i] = ls
        end
    end
    return g_seq, fasta_samples, position_used, mutations_tot_2
end


"""
    it check that the drivers mutation are not overwritten
"""
function not_pos_driver(pos_used::Vector{Any}, indexs::Vector{Int})
    for p in indexs[1]:indexs[2]
        if p ∈ pos_used
            return false
        end
    end
    return true
end


function not_pos_driver(pos_used::Vector{Any}, indexs::Int)
    if indexs ∈ pos_used
        return false
    else
        return true
    end
end

"""
    comment gemonic evolution INDEL + SVN
"""
function genomic_evolution(Seq_f::LongSequence,
                           rate_Indel::AbstractFloat,
                           size_indel_arr::Vector{Float64},
                           branch_length::AbstractFloat,
                           Model_Selector_matrix::DataFrame,
                           prob_A::Vector{Float64},
                           prob_C::Vector{Float64},
                           prob_G::Vector{Float64},
                           prob_T::Vector{Float64},
                           seed::MersenneTwister,
                           mutations_tot::DataFrame,
                           sample::Int;
                           position_used::Vector{Any}= [],
                           num_mut_driver::Int = -1,
                           muts = [0])

    sequence_father = copy(Seq_f)
    len_father = length(sequence_father)
    len_num_mut_driver = num_mut_driver

    As = findall(x -> x == 'A', string(sequence_father))
    n_A = length(As)
    Cs = findall(x -> x == 'C', string(sequence_father))
    n_C = length(Cs)
    Gs = findall(x -> x == 'G', string(sequence_father))
    n_G = length(Gs)
    Ts = findall(x -> x == 'T', string(sequence_father))
    n_T = length(Ts)
    curr_time = 0

    while curr_time <= branch_length

        λ_A = n_A * sum(Model_Selector_matrix[:,1])
        λ_C = n_C * sum(Model_Selector_matrix[:,2])
        λ_G = n_G * sum(Model_Selector_matrix[:,3])
        λ_T = n_T * sum(Model_Selector_matrix[:,4])
        λ_indel = (len_father+1) * rate_Indel
        λ_tot = sum([λ_A, λ_C, λ_G, λ_T, λ_indel])

        curr_time = curr_time + rand(seed, Exponential(1 / λ_tot), 1)[1]

        prob_vet = vcat(λ_A, λ_C, λ_G, λ_T, λ_indel) ./ λ_tot
        prob_cum = cumsum(prob_vet)
        k = rand(seed)

        mutation = collect(k .<= prob_cum)
        min_mutation = findfirst(mutation)

        if num_mut_driver >= 0
            if num_mut_driver == 0
                num_mut_driver = -1
            else
                mut_sub = muts[len_num_mut_driver - num_mut_driver + 1]

                num_mut_driver -= 1

            end
        end

        if min_mutation == 5 #indel
            e = rand(seed, ["insertion", "deletion"])

            init_pos = rand(seed, 1:len_father)

            if e == "insertion"
                init_pos = rand(seed, 1:len_father)
                prob_cum_size = cumsum(size_indel_arr)
                k = rand(seed)
                id_size_indel = collect(k .<= prob_cum_size)
                length_ins = findfirst(id_size_indel)
                if len_father - init_pos < length_ins
                    length_ins = len_father - init_pos
                end

                check = not_pos_driver(position_used, [init_pos, length_ins])
                trials = 1
                while check == false && trials < 100
                    init_pos = rand(seed, 1:len_father)
                    k = rand(seed)
                    id_size_indel = collect(k .<= prob_cum_size)
                    length_ins = findfirst(id_size_indel)
                    check = not_pos_driver(position_used,[init_pos, length_ins])
                    trials += 1
                end

                 init_pos_ins = rand(seed, 1:len_father)
                 insertion_sequence =
                                   sequence_father[init_pos:init_pos+length_ins]
                 new_sequence = sequence_father[1:init_pos_ins]

                 append!(new_sequence, insertion_sequence)
                 append!(new_sequence, sequence_father[init_pos_ins + 1 : end])
                 sequence_father = new_sequence
                 len_father = len_father + length_ins
                 fine = init_pos_ins + length_ins
                 if num_mut_driver >= 0
                     push!(mutations_tot, [[init_pos_ins, fine],
                                           "-",
                                           string(insertion_sequence),
                                           "insertion",
                                           mut_sub,
                                           sample])
                     [push!(position_used, p) for p in init_pos:fine]
                 else
                     push!(mutations_tot, [[init_pos_ins, fine],
                                           "-",
                                           string(insertion_sequence),
                                           "insertion",
                                           false,
                                           sample])
                 end


             else #deletion
                 prob_cum_size = cumsum(size_indel_arr)
                 k = rand(seed)
                 id_size_indel = collect(k .<= prob_cum_size)
                 length_del = findfirst(id_size_indel)
                 len_father = length(sequence_father)

                 check = not_pos_driver(position_used, [init_pos, length_del])
                 trials = 1
                 while check == false && trials < 100
                     init_pos = rand(seed, 1:len_father)
                     k = rand(seed)
                     id_size_indel = collect(k .<= prob_cum_size)
                     length_del = findfirst(id_size_indel)
                     check = not_pos_driver(position_used,
                                            [init_pos, length_del])
                     trials += 1
                 end

                 init_pos_del = rand(seed, 1:len_father)
                 end_pos_del = min(len_father, init_pos_del + length_del)
                 new_sequence = sequence_father[1:init_pos_del]

                 if end_pos_del != len_father
                     append!(new_sequence, sequence_father[end_pos_del + 1:end])
                 end

                 deletion_sequence = sequence_father[init_pos_del:end_pos_del]
                 sequence_father = new_sequence
                 len_father = len_father - length_del
                 if num_mut_driver >= 0
                     fine = init_pos_ins + length_ins
                     push!(mutations_tot, [[init_pos_del, end_pos_del],
                                           string(deletion_sequence),
                                           "-",
                                           "deletion",
                                           mut_sub,
                                           sample])
                 else
                     push!(mutations_tot, [[init_pos_del, end_pos_del],
                                           string(deletion_sequence),
                                           "-",
                                           "deletion",
                                           false,
                                           sample])
                 end
            end


        elseif min_mutation == 4 #T
            pos_mutation = rand(seed, Ts)

            #not overwrite driver mut
            check = not_pos_driver(position_used, pos_mutation)
            trials = 1
            while check == false && trials < length(Ts)
                pos_mutation = rand(seed, Ts)
                check = not_pos_driver(position_used,  pos_mutation)
                trials += 1
            end

            prob_cum_T = cumsum(prob_T)
            k = rand(seed)
            mutation = collect(k .<= prob_cum_T)
            ff = findfirst(mutation)
            new_nucleotide = collect(names(Model_Selector_matrix)[ff])[1]

            sequence_father[pos_mutation] = DNA(new_nucleotide)

            mut_id = string(pos_mutation)*"_T>"*string(new_nucleotide)
            if num_mut_driver >= 0
                push!(mutations_tot, [pos_mutation,
                                      "T",
                                      string(new_nucleotide),
                                      mut_id,
                                      mut_sub,
                                      sample])
                push!(position_used, pos_mutation)
            else
                push!(mutations_tot, [pos_mutation,
                                     "T",
                                     string(new_nucleotide),
                                     mut_id,
                                     false,
                                     sample])
            end


        elseif min_mutation == 3 #G
            pos_mutation = rand(seed, Gs)

            #not overwrite driver mut
            check = not_pos_driver(position_used, pos_mutation)
            trials = 1
            while check == false && trials < length(Gs)
                pos_mutation = rand(seed, Gs)
                check = not_pos_driver(position_used,  pos_mutation)
                trials += 1
            end

            prob_cum_G = cumsum(prob_G)
            k = rand(seed)
            mutation = collect(k .<= prob_cum_G)
            ff = findfirst(mutation)
            new_nucleotide = collect(names(Model_Selector_matrix)[ff])[1]
            sequence_father[pos_mutation] = DNA(new_nucleotide)

            mut_id = string(pos_mutation)*"_G>"*string(new_nucleotide)
            if num_mut_driver >= 0
                push!(mutations_tot, [pos_mutation,
                                      "G",
                                      string(new_nucleotide),
                                      mut_id,
                                      mut_sub,
                                      sample])
                push!(position_used, pos_mutation)
            else
                push!(mutations_tot, [pos_mutation,
                                      "G",
                                      string(new_nucleotide),
                                      mut_id,
                                      false,
                                      sample])
            end

        elseif min_mutation == 2 #C
            pos_mutation = rand(seed, Cs)

            #not overwrite driver mut
            check = not_pos_driver(position_used, pos_mutation)
            trials = 1
            while check == false && trials < length(Cs)
                pos_mutation = rand(seed, Cs)
                check = not_pos_driver(position_used,  pos_mutation)
                trials += 1
            end

            prob_cum_C = cumsum(prob_C)
            k = rand(seed)
            mutation = collect(k .<= prob_cum_C)
            ff = findfirst(mutation)
            new_nucleotide = collect(names(Model_Selector_matrix)[ff])[1]
            sequence_father[pos_mutation] = DNA(new_nucleotide)

            mut_id = string(pos_mutation)*"_C>"*string(new_nucleotide)
            if num_mut_driver >= 0
                push!(mutations_tot, [pos_mutation,
                                      "C",
                                      string(new_nucleotide),
                                      mut_id,
                                      mut_sub,
                                      sample])
                push!(position_used, pos_mutation)
            else
                push!(mutations_tot, [pos_mutation,
                                      "C",
                                      string(new_nucleotide),
                                      mut_id,
                                      false,
                                      sample])
            end


        elseif min_mutation == 1 #A
            pos_mutation = rand(seed, As)

            #not overwrite driver mut
            check = not_pos_driver(position_used, pos_mutation)
            trials = 1
            while check == false && trials < length(As)
                pos_mutation = rand(seed, As)
                check = not_pos_driver(position_used,  pos_mutation)
                trials += 1
            end

            prob_cum_A = cumsum(prob_A)
            k = rand(seed)
            mutation = collect(k .<= prob_cum_A)
            ff = findfirst(mutation)
            new_nucleotide = collect(names(Model_Selector_matrix)[ff])[1]
            sequence_father[pos_mutation] = DNA(new_nucleotide)
            mut_id = string(pos_mutation)*"_A>"*string(new_nucleotide)
            if num_mut_driver >= 0
                push!(mutations_tot, [pos_mutation,
                                      "A",
                                      string(new_nucleotide),
                                      mut_id,
                                      mut_sub,
                                      sample])
                push!(position_used, pos_mutation)
            else
                push!(mutations_tot, [pos_mutation,
                                      "A",
                                      string(new_nucleotide),
                                      mut_id,
                                      false,
                                      sample])
            end
        end

        #update values
        As = findall(x -> x == 'A', string(sequence_father))
        n_A = length(As)
        Cs = findall(x -> x == 'C', string(sequence_father))
        n_C = length(Cs)
        Gs = findall(x -> x == 'G', string(sequence_father))
        n_G = length(Gs)
        Ts = findall(x -> x == 'T', string(sequence_father))
        n_T = length(Ts)
    end

    if len_num_mut_driver != -1
         return sequence_father, num_mut_driver, position_used, mutations_tot
     else
         return sequence_father, mutations_tot
     end
end

"""
    compute distribution for size indel
"""
function size_indel_dist(len_g::Int,
                         size_indel::Int,
                         lavalette_par::AbstractFloat)

    if size_indel > len_g
        size_indel = len_g - 1
    end

    size_indel_arr = [( (u * size_indel) / (size_indel-u +1)) ^ (-lavalette_par)
                                                          for u in 1:size_indel]
    size_indel_arr = size_indel_arr ./ sum(size_indel_arr)
    return size_indel_arr
end

"""
    Molecular evolution with several substitution models
"""
function Molecular_evolution_NoISA(Tree::AbstractMetaGraph,
                                   Ref::LongSequence,
                                   Len::Int,
                                   Selector::String,
                                   params::IdDict,
                                   rate_Indel::AbstractFloat,
                                   size_indel::Int,
                                   seed::MersenneTwister,
                                   set_mut::Vector{Any},
                                   lavalette_par::AbstractFloat,
                                   approx_snv_indel::Int)

    Tree_SC = copy(Tree)
    set_prop!(Tree_SC, 1, :Subpop_Child, 1)

    mutations_tot = DataFrame(Pos = Any[],
                              Reference = String[],
                              Alternative = String[],
                              Mut_id = Any[],
                              Driver = Any[],
                              Sample = Any[])


    if rate_Indel != 0.0
        size_indel_arr = size_indel_dist(Len, size_indel, lavalette_par)
    else
        size_indel_arr = [0.0]
    end
    #Model_Selector
    Model_Selector_matrix = Q(Selector, params)
    if typeof(Model_Selector_matrix) == String
        return Ref, [], Tree_SC
    end

    prob_A = Model_Selector_matrix[:,1] ./ sum(Model_Selector_matrix[:,1])
    prob_C = Model_Selector_matrix[:,2] ./ sum(Model_Selector_matrix[:,2])
    prob_G = Model_Selector_matrix[:,3] ./ sum(Model_Selector_matrix[:,3])
    prob_T = Model_Selector_matrix[:,4] ./ sum(Model_Selector_matrix[:,4])

    position_used = []
    #max_time = max_time_nodes(Tree, get_leafs(Tree))
    for e in edges(Tree_SC)

        g_seq_e = LongDNA{4}()
        #g_seq_e = LongSequence()
        if has_prop(Tree_SC, src(e), :Fasta)
            g_seq_e = copy(get_prop(Tree_SC, src(e), :Fasta))
        else
            g_seq_e = copy(Ref)
        end

        branch_length = get_prop(Tree_SC, dst(e), :Time) -
                        get_prop(Tree_SC, src(e), :Time)
        #branch_length = branch_length/max_time
        subpop_father = get_prop(Tree_SC, src(e), :Subpop_Child)
        subpop_child = get_prop(Tree_SC, dst(e), :Subpop_Child)
        mut_f = set_mut[subpop_father]
        mut_child = set_mut[subpop_child]

        if length(mut_f) == length(mut_child)
            if approx_snv_indel == 0
                sequence, mutations_tot = genomic_evolution(g_seq_e,
                                             rate_Indel,
                                             size_indel_arr,
                                             branch_length,
                                             Model_Selector_matrix,
                                             prob_A,
                                             prob_C,
                                             prob_G,
                                             prob_T,
                                             seed,
                                             mutations_tot,
                                             dst(e))
            else
                sequence, mutations_tot = genomic_evolution_SNV(g_seq_e,
                                                   branch_length,
                                                   Model_Selector_matrix,
                                                   prob_A,
                                                   prob_C,
                                                   prob_G,
                                                   prob_T,
                                                   seed,
                                                   mutations_tot,
                                                   dst(e))
                if rate_Indel != 0.0
                    sequence, mutations_tot = genomic_evolution_INDEL(sequence,
                                                   rate_Indel,
                                                   size_indel_arr,
                                                   branch_length,
                                                   seed,
                                                   mutations_tot,
                                                   dst(e))
                end
            end

        else
            if typeof(mut_f) == String
                num_mut_driver = length(mut_child) - 1
            else
                num_mut_driver = length(mut_child) - length(mut_f)
            end
            new_muts = copy(mut_child)
            if typeof(mut_f) == String
                setdiff!(new_muts, [mut_f])
            else
                setdiff!(new_muts, mut_f)
            end

            muts_already_write = []
            #check not same mutation driver
            for nm in new_muts
                if nm ∈ mutations_tot.Driver
                    push!(muts_already_write,nm)
                end
            end
            setdiff!(new_muts, muts_already_write)
            num_mut_driver = length(new_muts)

            if approx_snv_indel == 0
                sequence, num_mut_driver, position_used, mutations_tot =
                         genomic_evolution(g_seq_e,
                                           rate_Indel,
                                           size_indel_arr,
                                           branch_length,
                                           Model_Selector_matrix,
                                           prob_A,
                                           prob_C,
                                           prob_G,
                                           prob_T,
                                           seed,
                                           mutations_tot,
                                           dst(e),
                                           position_used=position_used,
                                           num_mut_driver=num_mut_driver,
                                           muts=new_muts)
            else
                sequence, num_mut_driver, position_used, mutations_tot =
                          genomic_evolution_SNV(g_seq_e,
                                                branch_length,
                                                Model_Selector_matrix,
                                                prob_A,
                                                prob_C,
                                                prob_G,
                                                prob_T,
                                                seed,
                                                mutations_tot,
                                                dst(e),
                                                position_used=position_used,
                                                num_mut_driver=num_mut_driver,
                                                muts=new_muts)
                if num_mut_driver == -1
                    num_mut_driver = 0
                end
                if rate_Indel != 0.0
                    sequence, num_mut_driver, position_used, mutations_tot =
                        genomic_evolution_INDEL(sequence,
                                                rate_Indel,
                                                size_indel_arr,
                                                branch_length,
                                                seed,
                                                mutations_tot,
                                                dst(e),
                                                position_used=position_used,
                                                num_mut_driver=num_mut_driver,
                                                muts=new_muts)
                end

            end
            if num_mut_driver > 0
                pos_rand = sample(seed,
                                  1:length(sequence),
                                  num_mut_driver,
                                  replace = false)
                new_muts = new_muts[end-num_mut_driver:end]
                sequence, substitution, position_used, mutations_tot =
                               transform_genome(sequence,
                                                pos_rand,
                                                seed,
                                                mutations_tot,
                                                dst(e),
                                                funct = 1,
                                                position_used=position_used,
                                                muts_driver = new_muts)

                println("WARNING -> the mutational rates are low")

                println("Low mutational rate or very short branch length.
                Few mutations were generated in this branch, and driver
                mutations were added by force.")

                for sub in substitution
                    println(sub)
                end
            end
        end
        set_prop!(Tree_SC, dst(e), :Fasta, sequence)
    end

    # Return fasta of leaf nodes
    leafs = get_leafs(Tree_SC)
    fasta_samples = []
    for l in leafs
        f = get_prop(Tree_SC, l, :Fasta)
        push!(fasta_samples, f)
    end


    return Ref, fasta_samples, Tree_SC, mutations_tot
end

"""
Call function Molecular_evolution_NoISA
"""
function experiment_noISA(Tree::AbstractMetaGraph,
                          Len::Int,
                          Selector::String,
                          params::IdDict,
                          rate_Indel::AbstractFloat,
                          size_indel::Int,
                          seed::MersenneTwister,
                          set_mut::Vector{Any},
                          lavalette_par::AbstractFloat,
                          approx_snv_indel::Int;
                          frequency_dna::Vector{Float64} = [0.3,0.2,0.2])

        ## Create reference genome
        sw = SamplerWeighted(dna"ACGT", frequency_dna)
        Ref = randseq(seed, DNAAlphabet{4}(), sw, Len)
        #Ref = randdnaseq(seed, Len)
        rec = FASTA.Record("Reference", Ref)
        w = FASTA.Writer(open("Reference.fasta", "w"))
        write(w, rec)
        close(w)

        g_seq, fasta_samples, position_used, mutations_tot =
                                     Molecular_evolution_NoISA(Tree,
                                                               Ref,
                                                               Len,
                                                               Selector,
                                                               params,
                                                               rate_Indel,
                                                               size_indel,
                                                               seed,
                                                               set_mut,
                                                               lavalette_par,
                                                               approx_snv_indel)
        mutations_tot_2 = copy(mutations_tot)
        leafs = get_leafs(Tree)
        paths_tot = []
        for l in leafs
            yen_k = yen_k_shortest_paths(Tree, 1, l)
            path = yen_k.paths[1]
            push!(paths_tot, path)
        end
        paths_tot

        for i in 1:length(mutations_tot_2.Sample)
            sample = mutations_tot_2.Sample[i]
            if sample ∉ leafs
                ls = leafs[sample .∈ paths_tot]
                mutations_tot_2.Sample[i] = ls
            end
        end
        return g_seq, fasta_samples, position_used, mutations_tot_2
end

function experiment_noISA(Tree::AbstractMetaGraph,
                          path::String,
                          Selector::String,
                          params::IdDict,
                          rate_Indel::AbstractFloat,
                          size_indel::Int,
                          seed::MersenneTwister,
                          set_mut::Vector{Any},
                          lavalette_par::AbstractFloat,
                          approx_snv_indel::Int,
                          frequency_dna::Vector{Any} = [])

        ## load reference genome
        Ref = LongDNA{4}()
        #Ref = LongDNASeq()
        ## load reference genome
        open(FASTA.Reader, path) do reader
            for record in reader
                Ref = FASTX.sequence(record)
            end
        end

        Len = length(Ref)

        Ref, fasta_samples, position_used, mutations_tot =
                                    Molecular_evolution_NoISA(Tree,
                                                              Ref,
                                                              Len,
                                                              Selector,
                                                              params,
                                                              rate_Indel,
                                                              size_indel,
                                                              seed,
                                                              set_mut,
                                                              lavalette_par,
                                                              approx_snv_indel)

        mutations_tot_2 = copy(mutations_tot)
        leafs = get_leafs(Tree)
        paths_tot = []
        for l in leafs
            yen_k = yen_k_shortest_paths(Tree, 1, l)
            path = yen_k.paths[1]
            push!(paths_tot, path)
        end
        paths_tot

        for i in 1:length(mutations_tot_2.Sample)
            sample = mutations_tot_2.Sample[i]
            if sample ∉ leafs
                ls = leafs[sample .∈ paths_tot]
                mutations_tot_2.Sample[i] = ls
            end
        end
        return Ref, fasta_samples, position_used, mutations_tot_2
end

"""
    Compute probability for signature
"""
function linear_com_matrix(vector_change_points::Vector{Float64},
                           vector_activities::Matrix{Float64},
                           Matrix_sign::Matrix{String},
                           t_curr::Float64,
                           ratio_background_signature::Float64,
                           background_matrix::Vector{Float64})


        #change values in float
        Sign_matr = map((x) -> parse(Float64, x), Matrix_sign[:,2:end])

        #matrix 96x1, dot αᵢ for signature associated e
        sub_prob_matrix = []
        for time in 1:length(vector_change_points)
            sub_prob = Sign_matr * vector_activities[time, :]

            sub_prob = 0.5*(sub_prob * ratio_background_signature) +
                            ((1-ratio_background_signature) * background_matrix)
            push!(sub_prob_matrix, sub_prob)
        end
         #signature
         variation = Matrix_sign[:,1]
         labels = Array{String,2}(undef, 96, 6)
         i = 1

         for s in variation
             v = []
             s_1 = split(s, ['['])
             s_2 = split(s_1[2], [']'])
             s_3 = split(s_2[1], '>')
             start_tri = s_1[1]*s_3[1]*s_2[2]
             push!(v, start_tri)
             push!(v, s_1[1])#nucl sx
             push!(v, s_3[1])#nucl central
             push!(v, s_2[2])#nucl dx
             push!(v, s_3[2])#nucl final
             finish_tri = s_1[1]*s_3[2]*s_2[2]
             push!(v, finish_tri)
             labels[i, :] = v
             i = i+1
         end

         #add labels
         for l in 1:length(sub_prob_matrix)
             sub_prob_matrix[l] = hcat(labels, sub_prob_matrix[l])
         end

         return sub_prob_matrix, labels
end

"""
    Matrix that sum prob of pair trinucleotides
"""
function collapsed_matrix(sub_prob_matrix::Vector{Any})
    coll_sub_prob_tot = []
    labels = sub_prob_matrix[1][:, 1:6]

    for i in 1:length(sub_prob_matrix)
        value = map((x) -> convert(Float64, x), sub_prob_matrix[i][:,7])

        coll_sub_prob = Array{Any,2}(undef,0,2)
        for first in ["A","C","G","T"], second in ["A","C","G","T"]
            #take all element first[sub]second
            check1 = findall(x-> x == first, labels[:, 2])
            check1 = labels[check1, :]
            check2 = findall(x-> x == second, check1[:, 4])

            #add same values
            labels[check2, :]
            #slide and sum
            sum1 = []
            sum2 = []
            for i in 1:size(labels[check2, :])[1]
                row = labels[check2, :][i,:]
                if row[3] == "C"
                    push!(sum1, value[check2[i]])
                else
                    push!(sum2, value[check2[i]])
                end
            end
            label1 = first*"C"*second*"/"*first*"G"*second
            label2 = first*"A"*second*"/"*first*"T"*second
            coll_sub_prob = [coll_sub_prob;[label1 sum(sum1)]]
            coll_sub_prob = [coll_sub_prob;[label2 sum(sum2)]]
        end
        push!(coll_sub_prob_tot, coll_sub_prob)
    end

    return coll_sub_prob_tot #
end

"""
    Compute count trinucleotide
"""
function vector_count_genome(Ref::LongSequence,
                             labels::Vector{Any})
    Count_trinuclotide = []
    Count_trinuclotide = DataFrame(Labels = Any[],
                         Count = Int64[])
    Count_tri_split = DataFrame(Labels = Any[],
                         Count = Int64[])
    for l in labels
        l1 = split(l, '/')[1]
        l2 = split(l, '/')[2]
        n_l1 = length(collect(eachmatch(Regex(l1), convert(String, Ref))))
        n_l2 = length(collect(eachmatch(Regex(l2), convert(String, Ref))))
        push!(Count_trinuclotide, [l, n_l1 + n_l2])
        push!(Count_tri_split, [l1,n_l1])
        push!(Count_tri_split, [l2,n_l2])
    end
    return Count_trinuclotide, Count_tri_split
end

"""
Call function noISA for signature
"""
function experiment_noISA_sign(Tree::AbstractMetaGraph,
                               Len::Int, #length genome
                               Selector::String,#96-SBS
                               mut_rate_average::AbstractFloat,
                               rate_Indel::AbstractFloat,
                               size_indel::Int,
                               seed::MersenneTwister,
                               set_mut::Vector{Any},
                               lavalette_par::AbstractFloat,
                               used_sign::Vector{String},
                               vector_change_points::Vector{Float64},
                               vector_activities::Matrix{Float64},
                               ratio_background_signature::Float64;
                               frequency_dna::Vector{Float64} = [0.3,0.2,0.2])

        ## Create reference genome
        sw = SamplerWeighted(dna"ACGT", frequency_dna)
        Ref = randseq(seed, DNAAlphabet{4}(), sw, Len)
        #Ref = randdnaseq(seed, Len)
        rec = FASTA.Record("Reference", Ref)
        w = FASTA.Writer(open("Reference.fasta", "w"))
        write(w, rec)
        close(w)

        g_seq, fasta_samples, Tree, mutations_tot =
                    Molecular_evolution_NoISA_sign(Tree,
                                                   Ref,
                                                   Len, #length genome,
                                                   Selector,#96-SBS-37 or -38
                                                   mut_rate_average,
                                                   rate_Indel,
                                                   size_indel,
                                                   seed,
                                                   set_mut,
                                                   lavalette_par,
                                                   used_sign,
                                                   vector_change_points,
                                                   vector_activities,
                                                   ratio_background_signature)
        mutations_tot_2 = copy(mutations_tot)
        leafs = get_leafs(Tree)
        paths_tot = []
        for l in leafs
            yen_k = yen_k_shortest_paths(Tree, 1, l)
            path = yen_k.paths[1]
            push!(paths_tot, path)
        end
        paths_tot

        for i in 1:length(mutations_tot_2.Sample)
            sample = mutations_tot_2.Sample[i]
            if sample ∉ leafs
                ls = leafs[sample .∈ paths_tot]
                mutations_tot_2.Sample[i] = ls
            end
        end
        return g_seq, fasta_samples, Tree, mutations_tot_2
end

function experiment_noISA_sign(Tree::AbstractMetaGraph,
                               path::String, #length genome
                               Selector::String,#96-SBS-37 or -38
                               mut_rate_average::AbstractFloat,
                               rate_Indel::AbstractFloat,
                               size_indel::Int,
                               seed::MersenneTwister,
                               set_mut::Vector{Any},
                               lavalette_par::AbstractFloat,
                               used_sign::Vector{String},
                               vector_change_points::Vector{Float64},
                               vector_activities::Matrix{Float64},
                               ratio_background_signature::Float64;
                               frequency_dna::Vector{Any} = [])

        ## load reference genome
        Ref = LongDNA{4}()
        #Ref = LongDNASeq()
        open(FASTA.Reader, path) do reader
            for record in reader
                Ref = FASTX.sequence(record)
            end
        end

        Len = length(Ref)

        Ref, fasta_samples, Tree, mutations_tot =
                    Molecular_evolution_NoISA_sign(Tree,
                                                   Ref,
                                                   Len, #length genome,
                                                   Selector,#96-SBS-37 or -38
                                                   mut_rate_average,
                                                   rate_Indel,
                                                   size_indel,
                                                   seed,
                                                   set_mut,
                                                   lavalette_par,
                                                   used_sign,
                                                   vector_change_points,
                                                   vector_activities,
                                                   ratio_background_signature)

        mutations_tot_2 = copy(mutations_tot)
        leafs = get_leafs(Tree)
        paths_tot = []
        for l in leafs
            yen_k = yen_k_shortest_paths(Tree, 1, l)
            path = yen_k.paths[1]
            push!(paths_tot, path)
        end
        paths_tot

        for i in 1:length(mutations_tot_2.Sample)
            sample = mutations_tot_2.Sample[i]
            if sample ∉ leafs
                ls = leafs[sample .∈ paths_tot]
                mutations_tot_2.Sample[i] = ls
            end
        end
        return Ref, fasta_samples, Tree, mutations_tot_2
end


"""
    TO DO ["SBS-37","SBS-38"]
"""
function Molecular_evolution_NoISA_sign(Tree::AbstractMetaGraph,
                                        Ref::LongSequence,
                                        Len::Int, #length genome
                                        Selector::String,#96-SBS-37 or -38
                                        mut_rate_average::AbstractFloat,
                                        rate_Indel::AbstractFloat,
                                        size_indel::Int,
                                        seed::MersenneTwister,
                                        set_mut::Vector{Any},
                                        lavalette_par::AbstractFloat,
                                        used_sign::Vector{String},
                                        vector_change_points::Vector{Float64},
                                        vector_activities::Matrix{Float64},
                                        ratio_background_signature::Float64)

    Tree_SC = copy(Tree)
    set_prop!(Tree_SC, 1, :Subpop_Child, 1)

    mutations_tot = DataFrame(Pos = Any[],
                              Reference = String[],
                              Alternative = String[],
                              Mut_id = Any[],
                              Driver = Any[],
                              Sample = Any[])

    t_global = 0.0
    if rate_Indel != 0.0
        size_indel_arr = size_indel_dist(Len, size_indel, lavalette_par)
    else
        size_indel_arr = [0.0]
    end
    #SIGNATURE
    if Sys.iswindows()
        if Selector == "SBS-37"
            SBS = readdlm(".\\utility\\SBS_GRCh37.txt", String)
        else
            SBS = readdlm(".\\utility\\SBS_GRCh38.txt", String)
        end
    else
        if Selector == "SBS-37"
            SBS = readdlm("./utility/SBS_GRCh37.txt", String)
        else
            SBS = readdlm("./utility/SBS_GRCh38.txt", String)
        end
    end
    #labels
    names_sing = SBS[1, 2:end]
    #signature selected
    SBS = SBS[2:end ,:]
    idx = findall(x -> x ∈ used_sign, names_sing)
    idx = idx .+ 1
    pushfirst!(idx,1)
    Matrix_signature = SBS[:, idx]

    homologous_sub = ["G>T" "C>A";
                      "G>C" "C>G";
                      "G>A" "C>T";
                      "A>T" "T>A";
                      "A>G" "T>C";
                      "A>C" "T>G"]

    background_matrix = fill(1/96, 96)
    sub_prob_matrix, labels = J_Space.linear_com_matrix(vector_change_points,
                              vector_activities,
                              Matrix_signature,
                              t_global,
                              ratio_background_signature,
                              background_matrix)


    coll_sub_prob = J_Space.collapsed_matrix(sub_prob_matrix)

    position_used = []
    #max_time = max_time_nodes(Tree, get_leafs(Tree))
    for e in edges(Tree_SC)

        g_seq_e = LongDNA{4}()
        #g_seq_e = LongDNASeq()
        if has_prop(Tree_SC, src(e), :Fasta)
            g_seq_e = copy(get_prop(Tree_SC, src(e), :Fasta))
        else
            g_seq_e = copy(Ref)
        end
        t_global = get_prop(Tree_SC, src(e), :Time)
        time_possible = collect(t_global .>= vector_change_points)
        id_time = findlast(time_possible)
        if id_time === nothing
                id_time = length(vector_change_points)
        end
        branch_length = get_prop(Tree_SC, dst(e), :Time) -
                        get_prop(Tree_SC, src(e), :Time)
        #branch_length = branch_length/max_time
        subpop_father = get_prop(Tree_SC, src(e), :Subpop_Child)

        subpop_child = get_prop(Tree_SC, dst(e), :Subpop_Child)

        mut_f = set_mut[subpop_father]
        mut_child = set_mut[subpop_child]
        if length(mut_f) == length(mut_child)
                sequence, mutations_tot, id_time =
                            genomic_evolution_sign(g_seq_e,
                                                   branch_length,
                                                   coll_sub_prob,
                                                   seed,
                                                   position_used,
                                                   mutations_tot,
                                                   dst(e),
                                                   t_global,
                                                   vector_change_points,
                                                   vector_activities,
                                                   ratio_background_signature,
                                                   mut_rate_average,
                                                   sub_prob_matrix,
                                                   id_time,
                                                   homologous_sub)
                if rate_Indel != 0.0
                    sequence, mutations_tot = genomic_evolution_INDEL(sequence,
                                                       rate_Indel,
                                                       size_indel_arr,
                                                       branch_length,
                                                       seed,
                                                       mutations_tot,
                                                       dst(e))
                end

        else
            if typeof(mut_f) == String
                num_mut_driver = length(mut_child) - 1
            else
                num_mut_driver = length(mut_child) - length(mut_f)
            end
            new_muts = copy(mut_child)
            if typeof(mut_f) == String
                setdiff!(new_muts, [mut_f])
            else
                setdiff!(new_muts, mut_f)
            end
            muts_already_write = []
            #check not same mutation driver
            for nm in new_muts
                if nm ∈ mutations_tot.Driver
                    push!(muts_already_write,nm)
                end
            end
            setdiff!(new_muts, muts_already_write)
            num_mut_driver = length(new_muts)

            sequence, mutations_tot, id_time ,num_mut_driver, position_used =
                          genomic_evolution_sign(g_seq_e,
                                                 branch_length,
                                                 coll_sub_prob,
                                                 seed,
                                                 position_used,
                                                 mutations_tot,
                                                 dst(e),
                                                 t_global,
                                                 vector_change_points,
                                                 vector_activities,
                                                 ratio_background_signature,
                                                 mut_rate_average,
                                                 sub_prob_matrix,
                                                 id_time,
                                                 homologous_sub,
                                                 num_mut_driver=num_mut_driver,
                                                 muts=new_muts)
                if num_mut_driver == -1
                    num_mut_driver = 0
                end
                if rate_Indel != 0.0
                    sequence, num_mut_driver, position_used, mutations_tot =
                        genomic_evolution_INDEL(sequence,
                                                rate_Indel,
                                                size_indel_arr,
                                                branch_length,
                                                seed,
                                                mutations_tot,
                                                dst(e),
                                                position_used=position_used,
                                                num_mut_driver=num_mut_driver,
                                                muts=new_muts)
                end

            if num_mut_driver > 0
                pos_rand = sample(seed,
                                  1:length(sequence),
                                  num_mut_driver,
                                  replace = false)
                new_muts = new_muts[end-num_mut_driver+1:end]
                sequence, substitution, position_used, mutations_tot =
                               transform_genome(sequence,
                                                pos_rand,
                                                seed,
                                                mutations_tot,
                                                dst(e),
                                                funct = 1,
                                                position_used=position_used,
                                                muts_driver = new_muts)

                println("WARNING -> the mutational rates are low")
                println("Low mutational rate or very short branch length.
                Few mutations were generated in this branch, and driver
                mutations were added by force.")

                for sub in substitution
                    println(sub)
                end
            end
        end
        set_prop!(Tree_SC, dst(e), :Fasta, sequence)
    end

    # Return fasta of leaf nodes
    leafs = get_leafs(Tree_SC)
    fasta_samples = []
    for l in leafs
        f = get_prop(Tree_SC, l, :Fasta)
        push!(fasta_samples, f)
    end

    return Ref, fasta_samples, Tree_SC, mutations_tot
end

"""
    simulate signature on genome
"""
function genomic_evolution_sign(Seq_f::LongSequence,
                               branch_length::AbstractFloat,
                               coll_sub_prob::Vector{Any},
                               seed::MersenneTwister,
                               position_used::Vector{Any},
                               mutations_tot::DataFrame,
                               sample::Int,
                               t_global::Float64,
                               vector_change_points::Vector{Float64},
                               vector_activities::Matrix{Float64},
                               ratio_background_signature::Float64,
                               mut_rate_avg::Float64,
                               sub_prob_matrix::Vector{Any},
                               id_time::Int,
                               homologous_sub::Matrix{String};
                               num_mut_driver::Int = -1,
                               muts = [0])

      sequence_father = copy(Seq_f)
      len_num_mut_driver = num_mut_driver
      labels = sub_prob_matrix[1][:, 1:6]
      Count_trinuclotide, Count_tri_split =
                                   vector_count_genome(sequence_father[3:end-2],
                                                    coll_sub_prob[id_time][:,1])
      t_curr = 0
      while t_curr <= branch_length
          t_global_old = t_global
          coll_sub_prob_value = coll_sub_prob[id_time][:,2]
          coll_sub_prob_value = map((x) -> convert(Float64, x),
                                    coll_sub_prob_value)
          λ_tot = dot(coll_sub_prob_value,Count_trinuclotide.Count)*mut_rate_avg
          new_time = rand(seed, Exponential(1 / λ_tot), 1)[1]

          t_curr += new_time
          t_global += new_time

          #check driver mutation
          if num_mut_driver >= 0
              if num_mut_driver == 0
                  num_mut_driver = -1
              else
                  mut_sub = muts[len_num_mut_driver - num_mut_driver + 1]
                  num_mut_driver -= 1
              end
          end

          #selecting the substituted trinucleodite and norm
          middle_matrix = coll_sub_prob_value .* Count_trinuclotide.Count
          middle_matrix = middle_matrix .* mut_rate_avg
          Prob_matrix = middle_matrix ./ sum(middle_matrix)

          prob_cum = cumsum(Prob_matrix)
          k = rand(seed)
          trinucleotide = collect(k .<= prob_cum)
          id_trinucleotide = findfirst(trinucleotide)

          coll_sub_prob[id_time][id_trinucleotide, :][1]

          if isodd(id_trinucleotide)
              trinucl_single = [Count_tri_split[id_trinucleotide, :Count],
                           Count_tri_split[id_trinucleotide+1, :Count]]
          else
              trinucl_single = [Count_tri_split[id_trinucleotide-1,:Count],
                             Count_tri_split[id_trinucleotide, :Count]]
          end
          trinucl_single = trinucl_single ./ sum(trinucl_single)

          #choose trinucleotide
          prob_cum = cumsum(trinucl_single)
          k = rand(seed)
          trinucleotide = collect(k .<= prob_cum)
          single_trinucl = findfirst(trinucleotide)

          trinucleotide = split(
                Count_trinuclotide[id_trinucleotide, :][1], '/')[single_trinucl]
          range(m::RegexMatch) = m.offset .+ (0:length(m.match)-1)
          pos_range = [range(e) for e ∈ eachmatch(Regex(trinucleotide),
                                              string(sequence_father[3:end-2]))]

          occ = rand(seed, 1:length(pos_range))
          pos = pos_range[occ] .+ 2

          #not overwrite driver mut
          check = not_pos_driver(position_used, pos[2])
          trials = 1
          while check == false && trials < length(pos_range)
              mv = Count_tri_split[Count_tri_split.Labels .== trinucleotide,:]
              occ = rand(seed,1:mv.Count[1])
              pos = pos_range[occ]
              check = not_pos_driver(position_used, pos[2])
              trials += 1
          end

          check1 = findall(x-> x == string(trinucleotide[1]), labels[:, 2])
          middle_labels = labels[check1, :]
          check2 = findall(x-> x == string(trinucleotide[3]), labels[:, 4])

          id_ps = findall(x->x in check2, check1)
          possible_sub = sub_prob_matrix[id_time][check1[id_ps], :]

          if trinucleotide[2] in ['C','G']
              possible_sub = possible_sub[findall(x -> x .== "C",
                                                possible_sub[:,3]), :]
          else
              possible_sub = possible_sub[findall(x -> x .== "T",
                                                possible_sub[:,3]), :]
          end

          #norm
          prob_sub = possible_sub[:,7] ./sum(possible_sub[:,7])

          prob_cum = cumsum(prob_sub)
          k = rand(seed)
          bool_sub = collect(k .<= prob_cum)
          sub_id = findfirst(bool_sub)

          if trinucleotide[2] == 'G' || trinucleotide[2] == 'A'
                  sub = possible_sub[sub_id,:][3]*">"*possible_sub[sub_id,:][5]

                  id = findall(x -> x .== sub, homologous_sub[:,2])

                  sub_trinucl =  homologous_sub[id,1]

                  final_trinucl = trinucleotide[1]*
                                  sub_trinucl[1][3]*
                                  trinucleotide[3]
          else
                  final_trinucl = possible_sub[sub_id,:][6]
          end

          #updated trinucleodite_counts
          sub_seq = sequence_father[pos[2]-2:pos[2]+2]
          sub_seq_split = [sub_seq[i:i+2] for i in 1:length(sub_seq)-2]
          for tri in sub_seq_split
              pos_range = [range(e) for e ∈ eachmatch(Regex(trinucleotide),
                                              string(sequence_father[3:end-2]))]
              id = findall(x-> x == string(tri),
                                            Count_tri_split.Labels)[1]
              Count_tri_split[id,:].Count = length(pos_range)
          end

          new_sub_seq = sub_seq
          if typeof(final_trinucl) == String
              final_trinucl = LongSequence{DNAAlphabet{4}}(final_trinucl)
          end
          new_sub_seq[3] = final_trinucl[2]
          sub_seq_split = [new_sub_seq[i:i+2] for i in 1:length(new_sub_seq)-2]

          for tri in sub_seq_split
              id = findall(x-> x == string(tri),
                                             Count_tri_split.Labels)[1]
              pos_range = [range(e) for e ∈ eachmatch(Regex(trinucleotide),
                                              string(sequence_father[3:end-2]))]
              Count_tri_split[id,:].Count = length(pos_range)

          end

          mut = string(pos[2])*
                "_"*
                string(trinucleotide)*
                ">"*
                string(final_trinucl)

          if num_mut_driver >= 0
              push!(mutations_tot, [pos[2],
                                    string(trinucleotide),
                                    string(final_trinucl),
                                    mut,
                                    mut_sub,
                                    sample])
              push!(position_used, pos[2])
          else
              push!(mutations_tot, [pos[2],
                                    string(trinucleotide),
                                    string(final_trinucl),
                                    mut,
                                    false,
                                    sample])
          end
          sequence_father[pos] = final_trinucl

          if (vector_change_points[id_time] > t_global_old) &&
             (vector_change_points[id_time] <= t_global)
             id_time += 1
         end
      end
      if len_num_mut_driver != -1
          return sequence_father,
                 mutations_tot,
                 id_time,
                 num_mut_driver,
                 position_used
      else
          return sequence_father, mutations_tot, id_time
      end
end

"""
Simulate SNV on genome
"""
function genomic_evolution_SNV(Seq_f::LongSequence,
                               branch_length::AbstractFloat,
                               Model_Selector_matrix::DataFrame,
                               prob_A::Vector{Float64},
                               prob_C::Vector{Float64},
                               prob_G::Vector{Float64},
                               prob_T::Vector{Float64},
                               seed::MersenneTwister,
                               mutations_tot::DataFrame,
                               sample::Int;
                               position_used::Vector{Any} = [],
                               num_mut_driver::Int = -1,
                               muts = [0])

    sequence_father = copy(Seq_f)
    len_num_mut_driver = num_mut_driver

    As = findall(x -> x == 'A', string(sequence_father))
    n_A = length(As)
    Cs = findall(x -> x == 'C', string(sequence_father))
    n_C = length(Cs)
    Gs = findall(x -> x == 'G', string(sequence_father))
    n_G = length(Gs)
    Ts = findall(x -> x == 'T', string(sequence_father))
    n_T = length(Ts)

    curr_time = 0
    while curr_time <= branch_length
        #evaluate rate
        λ_A = n_A * sum(Model_Selector_matrix[:,1])
        λ_C = n_C * sum(Model_Selector_matrix[:,2])
        λ_G = n_G * sum(Model_Selector_matrix[:,3])
        λ_T = n_T * sum(Model_Selector_matrix[:,4])

        λ_tot = sum([λ_A, λ_C, λ_G, λ_T])

        curr_time = curr_time + rand(seed, Exponential(1 / λ_tot), 1)[1]
        prob_vet = vcat(λ_A, λ_C, λ_G, λ_T) ./ λ_tot
        prob_cum = cumsum(prob_vet)

        k = rand(seed)
        mutation = collect(k .<= prob_cum)
        min_mutation = findfirst(mutation)

        if num_mut_driver >= 0
            if num_mut_driver == 0
                num_mut_driver = -1
            else
                mut_sub = muts[len_num_mut_driver - num_mut_driver + 1]
                num_mut_driver -= 1
            end
        end

        if min_mutation == 4 #T
            pos_mutation = rand(seed, Ts)

            #not overwrite driver mut
            check = not_pos_driver(position_used, pos_mutation)
            trials = 1
            while check == false && trials < length(Ts)
                pos_mutation = rand(seed, Ts)
                check = not_pos_driver(position_used,  pos_mutation)
                trials += 1
            end

            prob_cum_T = cumsum(prob_T)
            k = rand(seed)
            mutation = collect(k .<= prob_cum_T)
            ff = findfirst(mutation)
            new_nucleotide = collect(names(Model_Selector_matrix)[ff])[1]
            sequence_father[pos_mutation] = DNA(new_nucleotide)

            mut_id = string(pos_mutation)*"_T>"*string(new_nucleotide)
            if num_mut_driver >= 0
                push!(mutations_tot, [pos_mutation,
                                      "T",
                                      string(new_nucleotide),
                                      mut_id,
                                      mut_sub,
                                      sample])
                push!(position_used, pos_mutation)
            else
                push!(mutations_tot, [pos_mutation,
                                      "T",
                                      string(new_nucleotide),
                                      mut_id,
                                      false,
                                      sample])
            end

            if new_nucleotide == "A"
             n_A=n_A+1
             push!(As, pos_mutation)
            end

            if new_nucleotide == "C"
             n_C=n_C+1
             push!(Cs, pos_mutation)
            end

            if new_nucleotide == "G"
             n_G=n_G+1
             push!(Gs, pos_mutation)
            end

            n_T=n_T-1
            filter!(p -> p != pos_mutation, Ts)

        elseif min_mutation == 3 #G
            pos_mutation = rand(seed, Gs)

            #not overwrite driver mut
            check = not_pos_driver(position_used, pos_mutation)
            trials = 1
            while check == false && trials < length(Gs)
                pos_mutation = rand(seed, Gs)
                check = not_pos_driver(position_used,  pos_mutation)
                trials += 1
            end

            prob_cum_G = cumsum(prob_G)
            k = rand(seed)
            mutation = collect(k .<= prob_cum_G)
            ff = findfirst(mutation)
            new_nucleotide = collect(names(Model_Selector_matrix)[ff])[1]
            sequence_father[pos_mutation] = DNA(new_nucleotide)

            mut_id = string(pos_mutation)*"_G>"*string(new_nucleotide)
            if num_mut_driver >= 0
                push!(mutations_tot, [pos_mutation,
                                      "G",
                                      string(new_nucleotide),
                                      mut_id,
                                      mut_sub,
                                      sample])
                push!(position_used, pos_mutation)
            else
                push!(mutations_tot, [pos_mutation,
                                      "G",
                                      string(new_nucleotide),
                                      mut_id,
                                      false,
                                      sample])
            end

            if new_nucleotide == "A"
             n_A=n_A+1
             push!(As, pos_mutation)
            end

            if new_nucleotide == "C"
             n_C=n_C+1
             push!(Cs, pos_mutation)
            end

            if new_nucleotide == "T"
             n_T=n_T+1
             push!(Ts, pos_mutation)
            end

            n_G=n_G-1
            filter!(p -> p != pos_mutation, Gs)

        elseif min_mutation == 2 #C
            pos_mutation = rand(seed, Cs)

            #not overwrite driver mut
            check = not_pos_driver(position_used, pos_mutation)
            trials = 1
            while check == false && trials < length(Cs)
                pos_mutation = rand(seed, Cs)
                check = not_pos_driver(position_used,  pos_mutation)
                trials += 1
            end

            prob_cum_C = cumsum(prob_C)
            k = rand(seed)
            mutation = collect(k .<= prob_cum_C)
            ff = findfirst(mutation)
            new_nucleotide = collect(names(Model_Selector_matrix)[ff])[1]
            sequence_father[pos_mutation] = DNA(new_nucleotide)

            mut_id = string(pos_mutation)*"_C>"*string(new_nucleotide)
            if num_mut_driver >= 0
                push!(mutations_tot, [pos_mutation,
                                      "C",
                                      string(new_nucleotide),
                                      mut_id,
                                      mut_sub,
                                      sample])
                push!(position_used, pos_mutation)
            else
                push!(mutations_tot, [pos_mutation,
                                      "C",
                                      string(new_nucleotide),
                                      mut_id,
                                      false,
                                      sample])
            end

            if new_nucleotide == "A"
                n_A=n_A+1
                push!(As, pos_mutation)
            end

            if new_nucleotide == "G"
                n_G=n_G+1
                push!(Gs, pos_mutation)
            end

            if new_nucleotide == "T"
                n_T=n_T+1
                push!(Ts, pos_mutation)
            end

            n_C=n_C-1
            filter!(p -> p != pos_mutation,  Cs)

        elseif min_mutation == 1 #A
            pos_mutation = rand(seed, As)

            #not overwrite driver mut
            check = not_pos_driver(position_used, pos_mutation)
            trials = 1
            while check == false && trials < length(As)
                pos_mutation = rand(seed, As)
                check = not_pos_driver(position_used,  pos_mutation)
                trials += 1
            end

            prob_cum_A = cumsum(prob_A)
            k = rand(seed)
            mutation = collect(k .<= prob_cum_A)
            ff = findfirst(mutation)
            new_nucleotide = collect(names(Model_Selector_matrix)[ff])[1]
            sequence_father[pos_mutation] = DNA(new_nucleotide)

            mut_id = string(pos_mutation)*"_A>"*string(new_nucleotide)
            if num_mut_driver >= 0
                push!(mutations_tot, [pos_mutation,
                                      "A",
                                      string(new_nucleotide),
                                      mut_id,
                                      mut_sub,
                                      sample])
                push!(position_used, pos_mutation)
            else
                push!(mutations_tot, [pos_mutation,
                                      "A",
                                      string(new_nucleotide),
                                      mut_id,
                                      false,
                                      sample])
            end

            if new_nucleotide == "C"
                n_C=n_C+1
                push!(Cs, pos_mutation)
            end

            if new_nucleotide == "G"
                n_G=n_G+1
                push!(Gs, pos_mutation)
            end

            if new_nucleotide == "T"
                n_T=n_T+1
                push!(Ts, pos_mutation)
            end

            n_A=n_A-1
            filter!(p -> p != pos_mutation, As)
        end

    end
    if len_num_mut_driver != -1
         return sequence_father, num_mut_driver, position_used, mutations_tot
     else
         return sequence_father, mutations_tot
     end
end

"""
Simulate INDEL on genome
"""
function genomic_evolution_INDEL(Seq_f::LongSequence,
                                 rate_Indel::AbstractFloat,
                                 size_indel_arr::Vector{Float64},
                                 branch_length::AbstractFloat,
                                 seed::MersenneTwister,
                                 mutations_tot::DataFrame,
                                 sample::Int;
                                 position_used::Vector{Any} = [],
                                 num_mut_driver::Int = -1,
                                 muts = [0])


    sequence_father = copy(Seq_f)
    len_father = length(sequence_father)
    len_num_mut_driver = num_mut_driver
    curr_time = 0.0

    while curr_time <= branch_length
        λ_tot = (len_father+1) * rate_Indel
        curr_time = curr_time + rand(seed, Exponential(1 / λ_tot), 1)[1]
        e = rand(seed, ["insertion", "deletion"])

        #choose initial position
        init_pos = rand(seed, 1:len_father)
        if num_mut_driver >= 0
            if num_mut_driver == 0
                num_mut_driver = -1
            else
                mut_sub = muts[len_num_mut_driver - num_mut_driver + 1]
                num_mut_driver -= 1
            end
        end
        if e == "insertion"
            prob_cum_size = cumsum(size_indel_arr)
            k = rand(seed)
            id_size_indel = collect(k .<= prob_cum_size)
            length_ins = findfirst(id_size_indel)

            if len_father - init_pos < length_ins
                length_ins = len_father - init_pos
            end

            #not overwrite driver mut
            check = not_pos_driver(position_used, [init_pos, length_ins])
            trials = 1
            while check == false && trials < 100
                init_pos = rand(seed, 1:len_father)
                k = rand(seed)
                id_size_indel = collect(k .<= prob_cum_size)
                length_ins = findfirst(id_size_indel)
                check = not_pos_driver(position_used,[init_pos, length_ins])
                trials += 1
            end

            init_pos_ins = rand(seed, 1:len_father)
            insertion_sequence = sequence_father[init_pos:init_pos+length_ins]
            new_sequence = sequence_father[1:init_pos_ins]

            append!(new_sequence, insertion_sequence)
            append!(new_sequence, sequence_father[init_pos_ins + 1 : end])
            sequence_father = new_sequence
            len_father = len_father + length_ins
            fine = init_pos_ins + length_ins
            if num_mut_driver >= 0
                push!(mutations_tot, [[init_pos_ins, fine],
                                      "-",
                                      string(insertion_sequence),
                                      "insertion",
                                      mut_sub,
                                      sample])
                [push!(position_used, p) for p in init_pos:fine]
            else
                push!(mutations_tot, [[init_pos_ins, fine],
                                      "-",
                                      string(insertion_sequence),
                                      "insertion",
                                      false,
                                      sample])
            end


        else #deletion
            prob_cum_size = cumsum(size_indel_arr)
            k = rand(seed)
            id_size_indel = collect(k .<= prob_cum_size)
            length_del = findfirst(id_size_indel)

            #not overwrite driver mut
            check = not_pos_driver(position_used, [init_pos, length_del])
            trials = 1
            while check == false && trials < 100
                init_pos = rand(seed, 1:len_father)
                k = rand(seed)
                id_size_indel = collect(k .<= prob_cum_size)
                length_del = findfirst(id_size_indel)
                check = not_pos_driver(position_used,  [init_pos, length_del])
                trials += 1
            end

            init_pos_del = rand(seed, 1:length(sequence_father))
            end_pos_del = min(len_father, init_pos_del + length_del)
            new_sequence = sequence_father[1:init_pos_del]
            if end_pos_del != len_father
                append!(new_sequence, sequence_father[end_pos_del + 1:end])
            end

            deletion_sequence = sequence_father[init_pos_del:end_pos_del]
            sequence_father = new_sequence
            len_father = len_father - length_del

            if num_mut_driver >= 0
                fine = init_pos_ins + length_ins
                push!(mutations_tot, [[init_pos_del, end_pos_del],
                                       string(deletion_sequence),
                                       "-",
                                       "deletion",
                                       mut_sub,
                                       sample])
            else
                push!(mutations_tot, [[init_pos_del, end_pos_del],
                                      string(deletion_sequence),
                                      "-",
                                      "deletion",
                                      false,
                                      sample])
            end
        end
    end
    if len_num_mut_driver != -1
         return sequence_father, num_mut_driver, position_used, mutations_tot
     else
         return sequence_father, mutations_tot
     end
end
    ### end of file -- SingleCellExperiment.jl
