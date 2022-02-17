### -*- Mode: Julia -*-

### SingleCellExperiment.jl
import BioSequences
using FASTX

"""
Save fasta.
"""
function save_Fasta_W(Ref::LongDNASeq,
                      fasta_samples::Vector{Any},
                      tree::AbstractMetaGraph,
                      path_save_file::String)

      leafs = get_leafs(tree)
      mkpath(path_save_file * "\\Fasta output") # Create folder
      path_for_fasta = path_save_file * "\\Fasta output"
      for i in 1:length(leafs)
          w = FASTA.Writer(open(path_for_fasta
                                * "\\sample"
                                * string(leafs[i])
                                * ".fasta",
                                "w"))
          rec = FASTA.Record("Sample"
                             * string(leafs[i]),
                             fasta_samples[i])
          write(w, rec)
          close(w)
      end
      w = FASTA.Writer(open(path_for_fasta
                            * "\\reference"
                            * ".fasta",
                            "w"))
      rec = FASTA.Record("Reference", Ref)
      write(w, rec)
      close(w)
end

"""
Save fasta.
"""
function save_Fasta_L(Ref::LongDNASeq,
                      fasta_samples::Vector{Any},
                      tree::AbstractMetaGraph,
                      path_save_file::String)

      leafs = get_leafs(tree)
      mkpath(path_save_file * "/Fasta output") # Create folder
      path_for_fasta = path_save_file * "/Fasta output"
      for i in 1:length(leafs)
          w = FASTA.Writer(open(path_for_fasta
                                * "/sample"
                                * string(leafs[i])
                                * ".fasta",
                                "w"))
          rec = FASTA.Record("Sample"
                             * string(leafs[i]),
                             fasta_samples[i])
          write(w, rec)
          close(w)
      end
      w = FASTA.Writer(open(path_for_fasta
                            * "/reference"
                            * ".fasta",
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
                       seed::MersenneTwister)

    mutations = []
    λ = Neutral_mut_rate * length_ROI

    for e in edges(Tree)
        Tfinal = get_prop(Tree, dst(e), :Time) - get_prop(Tree, src(e), :Time)
        t_curr = 0
        n_mut = 0

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
function transform_genome(genome::LongDNASeq,
                          pos::Vector{Any},
                          seed::MersenneTwister;
                          funct::Int = 0,
                          mutation_driver::Dict = Dict{}(),
                          position_used::Vector{Any} = [])

    transition_matrix =
        DataFrame(A = [0.0, 0.33333333, 0.33333333, 0.33333333],
                  C = [0.33333333, 0.0, 0.33333333, 0.33333333],
                  G = [0.33333333, 0.33333333, 0.0, 0.33333333],
                  T = [0.33333333, 0.33333333, 0.33333333, 0.0])
    substitution = []

    for p in pos
        nucleotide = genome[p]

        prob_cum = cumsum(transition_matrix[!, string(nucleotide)])

        k = round( rand(seed), digits = 7)
        target_nucl = collect(k .<= prob_cum)
        min = findfirst(target_nucl)

        new_nucleotide = collect(names(transition_matrix)[min])[1] # type char

        genome[p] = DNA(new_nucleotide)
        if funct == 1
            sub = string(p)*"_"*string(nucleotide)*">"*string(new_nucleotide)
            mutation_driver[p] = sub
            push!(substitution, sub)
            push!(position_used, p)
        end
    end

    if funct == 1
        return genome, substitution, position_used, mutation_driver
    else
        return genome
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
Creates a input for tool ART -> FASTA file (WITHOUT FASTA).
"""
function Molecular_evolution_ISA(Tree::AbstractMetaGraph,
                             neural_mut_rate::Float64,
                             seed::MersenneTwister,
                             len_ROI::Int,
                             set_mut::Vector{Any})
    ## Create reference genome
    g_seq = randdnaseq(seed, len_ROI)
    rec = FASTA.Record("Reference", g_seq)
    w = FASTA.Writer(open("Reference.fasta", "w"))
    write(w, rec)
    close(w)
    mutation_driver = Dict()

    ## Compute mutations ∀ node
    for i in [1,2,3]
        check = false
        n_mutations = evolution_seq(Tree, neural_mut_rate, len_ROI, seed)
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

            g_seq_e = LongDNASeq()

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

            g_seq_e = transform_genome(g_seq_e, pos_edge, seed)

            if length(mut_f) != length(mut_child)
                num_mut_driver = length(mut_child) - length(mut_f)
                new_muts = copy(mut_child)
                filter!(m -> m ∉ mut_f, new_muts)
                pos_edge_d = []
                for i = 1:length(new_muts)
                    pos = rand(seed, possible_position)
                    delete!(possible_position, pos)
                    push!(pos_edge_d, pos)
                end
                g_seq_e, sub, position_used, mutation_driver =
                               transform_genome(g_seq_e,
                                                pos_edge_d,
                                                seed,
                                                funct=1,
                                                mutation_driver=mutation_driver,
                                                position_used=position_used)
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

        return g_seq, fasta_samples, position_used, mutation_driver
end

"""
Creates a input for tool ART -> FASTA file (WITH REF FASTA).
"""
function Molecular_evolution_ISA(Tree::AbstractMetaGraph,
                                 neural_mut_rate::Float64,
                                 seed::MersenneTwister,
                                 path::String,
                                 set_mut::Vector{Any})
    ## load reference genome
    g_seq = LongDNASeq()
    open(FASTA.Reader, path) do reader
        for record in reader
            g_seq = FASTX.sequence(record)
        end
    end

    len_ROI = length(g_seq)
    mutation_driver = Dict()

    ## Compute mutations ∀ node
    for i in [1,2,3]
        check = false
        n_mutations = evolution_seq(Tree, neural_mut_rate, len_ROI, seed)
        check = check_num_path_mutation(Tree, len_ROI)

        if check == true && i >= 3
            return "ERROR: mutations are more than length ROI -> (file fasta)."
        end
    end

        possible_position = Set(1:len_ROI)

        position_used = []

        # create fasta ∀ nodes
        for e in edges(Tree)

            g_seq_e = LongDNASeq()

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

            g_seq_e = transform_genome(g_seq_e, pos_edge, seed)

            if length(mut_f) != length(mut_child)
                num_mut_driver = length(mut_child) - length(mut_f)
                new_muts = copy(mut_child)
                filter!(m -> m ∉ mut_f, new_muts)
                pos_edge_d = []
                for i = 1:length(new_muts)
                    pos = rand(seed, possible_position)
                    delete!(possible_position, pos)
                    push!(pos_edge_d, pos)
                    push!(position_used, pos)
                end
                g_seq_e,sub, position_used, mutation_driver =
                               transform_genome(g_seq_e,
                                                pos_edge_d,
                                                seed,
                                                funct=1,
                                                mutation_driver=mutation_driver,
                                                position_used=position_used)
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

        return g_seq, fasta_samples, position_used, mutation_driver
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
    comment gemonic evolution
"""
function genomic_evolution(Seq_f::LongDNASeq,
                           rate_Indel::AbstractFloat,
                           size_indel_arr::Vector{Float64},
                           branch_length::AbstractFloat,
                           Model_Selector_matrix::DataFrame,
                           prob_A::Vector{Float64},
                           prob_C::Vector{Float64},
                           prob_G::Vector{Float64},
                           prob_T::Vector{Float64},
                           seed::MersenneTwister,
                           position_used::Vector{Any};
                           mutation_driver::Dict = Dict{}(),
                           num_mut_driver::Int = -1,
                           muts::Vector{Int} = [0])

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

                 if num_mut_driver >= 0
                     fine = init_pos_ins + length_ins
                     mutation_driver[mut_sub] = "ins $init_pos_ins-$fine"
                     [push!(position_used, p) for p in init_pos:fine]
                 end


             else #deletion
                 prob_cum_size = cumsum(size_indel_arr)
                 k = rand(seed)
                 id_size_indel = collect(k .<= prob_cum_size)
                 length_del = findfirst(id_size_indel)


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
                     mutation_driver[mut_sub] = "del $init_pos_del-$end_pos_del"
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

            if num_mut_driver >= 0
                mutation_driver[mut_sub] = string(pos_mutation)
                                           * "_T>"
                                           * string(new_nucleotide)
                push!(position_used, pos_mutation)
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

            if num_mut_driver >= 0
                mutation_driver[mut_sub] = string(pos_mutation)
                                           * "_G>"
                                           * string(new_nucleotide)
                push!(position_used, pos_mutation)
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

            if num_mut_driver >= 0
                mutation_driver[mut_sub] = string(pos_mutation)
                                           *"_C>"
                                           *string(new_nucleotide)
                push!(position_used, pos_mutation)
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
            if num_mut_driver >= 0
                mutation_driver[mut_sub] = string(pos_mutation)
                                           *"_A>"
                                           *string(new_nucleotide)
                push!(position_used, pos_mutation)
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
         return sequence_father, num_mut_driver, position_used, mutation_driver
     else
         return sequence_father
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
    Molecular evolution with several substitution models (With ref)
"""
function Molecular_evolution_NoISA(Tree::AbstractMetaGraph,
                                   path::String,
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
    Ref = LongDNASeq()
    mutation_driver = Dict()

    ## load reference genome
    open(FASTA.Reader, path) do reader
        for record in reader
            Ref = FASTX.sequence(record)
        end
    end

    size_indel_arr = size_indel_dist(length(Ref), size_indel, lavalette_par)

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
    for e in edges(Tree_SC)

        g_seq_e = LongDNASeq()

        if has_prop(Tree_SC, src(e), :Fasta)
            g_seq_e = copy(get_prop(Tree_SC, src(e), :Fasta))
        else
            g_seq_e = copy(Ref)
        end

        branch_length = get_prop(Tree_SC, dst(e), :Time) -
                        get_prop(Tree_SC, src(e), :Time)
        subpop_father = get_prop(Tree_SC, src(e), :Subpop_Child)
        subpop_child = get_prop(Tree_SC, dst(e), :Subpop_Child)
        mut_f = set_mut[subpop_father]
        mut_child = set_mut[subpop_child]

        if length(mut_f) == length(mut_child)
            if approx_snv_indel == 0
                sequence = genomic_evolution(g_seq_e,
                                             rate_Indel,
                                             size_indel_arr,
                                             branch_length,
                                             Model_Selector_matrix,
                                             prob_A,
                                             prob_C,
                                             prob_G,
                                             prob_T,
                                             seed,
                                             position_used)
            else
                sequence_1 = genomic_evolution_SNV(g_seq_e,
                                                   branch_length,
                                                   Model_Selector_matrix,
                                                   prob_A,
                                                   prob_C,
                                                   prob_G,
                                                   prob_T,
                                                   seed,
                                                   position_used)
                sequence = genomic_evolution_INDEL(sequence_1,
                                                   rate_Indel,
                                                   size_indel_arr,
                                                   branch_length,
                                                   seed,
                                                   position_used)
            end

        else
            num_mut_driver = length(mut_child) - length(mut_f)
            new_muts = copy(mut_child)
            filter!(m -> m ∉ mut_f, new_muts)

            if approx_snv_indel == 0
                sequence, num_mut_driver, position_used, mutation_driver =
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
                                            position_used,
                                            mutation_driver=mutation_driver,
                                            num_mut_driver=num_mut_driver,
                                            muts=new_muts)
            else
                sequence_1, num_mut_driver, position_used, mutation_driver =
                          genomic_evolution_SNV(g_seq_e,
                                                branch_length,
                                                Model_Selector_matrix,
                                                prob_A,
                                                prob_C,
                                                prob_G,
                                                prob_T,
                                                seed,
                                                position_used,
                                                mutation_driver=mutation_driver,
                                                num_mut_driver=num_mut_driver,
                                                muts=new_muts)
                if num_mut_driver == -1
                    num_mut_driver = 0
                end
                sequence, num_mut_driver, position_used, mutation_driver =
                        genomic_evolution_INDEL(sequence_1,
                                                rate_Indel,
                                                size_indel_arr,
                                                branch_length,
                                                seed,
                                                position_used,
                                                mutation_driver=mutation_driver,
                                                num_mut_driver=num_mut_driver,
                                                muts=new_muts)

            end
            if num_mut_driver > 0
                pos_rand = sample(seed,
                                  1:length(sequence),
                                  num_mut_driver,
                                  replace = false)
                sequence, substitution, position_used, mutation_driver =
                               transform_genome(sequence,
                                                pos_rand,
                                                seed,
                                                funct = 1,
                                                mutation_driver=mutation_driver,
                                                position_used=position_used)

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


    return Ref, fasta_samples, Tree_SC, mutation_driver
end

"""
    Molecular evolution with several substitution models (Without ref)
"""
function Molecular_evolution_NoISA(Tree::AbstractMetaGraph,
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

    ## Create reference genome
    Ref = randdnaseq(seed, Len)
    rec = FASTA.Record("Reference", Ref)
    w = FASTA.Writer(open("Reference.fasta", "w"))
    write(w, rec)
    close(w)

    mutation_driver = Dict()

    size_indel_arr = size_indel_dist(Len, size_indel, lavalette_par)

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
    for e in edges(Tree_SC)
        g_seq_e = LongDNASeq()

        if has_prop(Tree_SC, src(e), :Fasta)
            g_seq_e = copy(get_prop(Tree_SC, src(e), :Fasta))
        else
            g_seq_e = copy(Ref)
        end

        branch_length = get_prop(Tree_SC, dst(e), :Time) -
                        get_prop(Tree_SC, src(e), :Time)
        subpop_father = get_prop(Tree_SC, src(e), :Subpop_Child)

        subpop_child = get_prop(Tree_SC, dst(e), :Subpop_Child)

        mut_f = set_mut[subpop_father]
        mut_child = set_mut[subpop_child]

        if length(mut_f) == length(mut_child)
            if approx_snv_indel == 0
                sequence = genomic_evolution(g_seq_e,
                                             rate_Indel,
                                             size_indel_arr,
                                             branch_length,
                                             Model_Selector_matrix,
                                             prob_A,
                                             prob_C,
                                             prob_G,
                                             prob_T,
                                             seed,
                                             position_used)
            else
                sequence_1 = genomic_evolution_SNV(g_seq_e,
                                                   branch_length,
                                                   Model_Selector_matrix,
                                                   prob_A,
                                                   prob_C,
                                                   prob_G,
                                                   prob_T,
                                                   seed,
                                                   position_used)
                sequence = genomic_evolution_INDEL(sequence_1,
                                                   rate_Indel,
                                                   size_indel_arr,
                                                   branch_length,
                                                   seed,
                                                   position_used)
            end

        else
            num_mut_driver = length(mut_child) - length(mut_f)
            new_muts = copy(mut_child)
            filter!(m -> m ∉ mut_f, new_muts)

            if approx_snv_indel == 0
                sequence, num_mut_driver, position_used, mutation_driver =
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
                                                position_used,
                                                mutation_driver=mutation_driver,
                                                num_mut_driver=num_mut_driver,
                                                muts=new_muts)
            else
                sequence_1, num_mut_driver, position_used, mutation_driver =
                          genomic_evolution_SNV(g_seq_e,
                                                branch_length,
                                                Model_Selector_matrix,
                                                prob_A,
                                                prob_C,
                                                prob_G,
                                                prob_T,
                                                seed,
                                                position_used,
                                                mutation_driver=mutation_driver,
                                                num_mut_driver=num_mut_driver,
                                                muts=new_muts)
                if num_mut_driver == -1
                    num_mut_driver = 0
                end
                sequence, num_mut_driver, position_used, mutation_driver =
                        genomic_evolution_INDEL(sequence_1,
                                                rate_Indel,
                                                size_indel_arr,
                                                branch_length,
                                                seed,
                                                position_used,
                                                mutation_driver=mutation_driver,
                                                num_mut_driver=num_mut_driver,
                                                muts=new_muts)
            end

            if num_mut_driver > 0
                pos_rand = sample(seed,
                                  1:length(sequence),
                                  num_mut_driver,
                                  replace = false)
                sequence, substitution, position_used, mutation_driver =
                               transform_genome(sequence,
                                                pos_rand,
                                                seed,
                                                funct = 1,
                                                mutation_driver=mutation_driver,
                                                position_used=position_used)

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

    return Ref, fasta_samples, Tree_SC, mutation_driver
end



function genomic_evolution_SNV(Seq_f::LongDNASeq,
                               branch_length::AbstractFloat,
                               Model_Selector_matrix::DataFrame,
                               prob_A::Vector{Float64},
                               prob_C::Vector{Float64},
                               prob_G::Vector{Float64},
                               prob_T::Vector{Float64},
                               seed::MersenneTwister,
                               position_used::Vector{Any};
                               mutation_driver::Dict = Dict{}(),
                               num_mut_driver::Int = -1,
                               muts::Vector{Int} = [0])

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

            if num_mut_driver >= 0
                mutation_driver[mut_sub] = string(pos_mutation)
                                           * "_T>"
                                           * string(new_nucleotide)
                push!(position_used, pos_mutation)
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

            if num_mut_driver >= 0
                mutation_driver[mut_sub] = string(pos_mutation)
                                           *"_G>"
                                           *string(new_nucleotide)
                push!(position_used, pos_mutation)

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

            if num_mut_driver >= 0
                mutation_driver[mut_sub] = string(pos_mutation)
                                           * "_C>"
                                           * string(new_nucleotide)
                push!(position_used, pos_mutation)
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

            if num_mut_driver >= 0
                mutation_driver[mut_sub] = string(pos_mutation)
                                           * "_A>"
                                           * string(new_nucleotide)
                push!(position_used, pos_mutation)
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
         return sequence_father, num_mut_driver, position_used, mutation_driver
     else
         return sequence_father
     end
end

function genomic_evolution_INDEL(Seq_f::LongDNASeq,
                                 rate_Indel::AbstractFloat,
                                 size_indel_arr::Vector{Float64},
                                 branch_length::AbstractFloat,
                                 seed::MersenneTwister,
                                 position_used::Vector{Any};
                                 mutation_driver::Dict = Dict{}(),
                                 num_mut_driver::Int = -1,
                                 muts::Vector{Int} = [0])


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

            if num_mut_driver >= 0
                mutation_driver[mut_sub] = "ins $init_pos_ins-"*string(init_pos_ins + length_ins)
                [push!(position_used, p) for p in init_pos:init_pos+length_ins]
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
                mutation_driver[mut_sub] = "del $init_pos_del- $end_pos_del"
            end
        end
    end
    if len_num_mut_driver != -1
         return sequence_father, num_mut_driver, position_used, mutation_driver
     else
         return sequence_father
     end
end
    ### end of file -- SingleCellExperiment.jl
