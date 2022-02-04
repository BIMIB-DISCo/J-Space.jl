### -*- Mode: Julia -*-

### SingleCellExperiment.jl
import BioSequences
using FASTX

"""
Save fasta.
"""
function save_Fasta(Ref::LongDNASeq, fasta_samples::Vector{Any},
                     tree::AbstractMetaGraph, path_save_file::String)
      ## write fasta on files if single_cell is true
      leafs = get_leafs(tree)
      mkpath(path_save_file*"\\Fasta output") # Create folder
      path_for_fasta = path_save_file*"\\Fasta output"
      for i in 1:length(leafs)
          ## Windows-ism!!!!
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
Computes number of mutations for all edges.
"""
function evolution_seq(Tree::AbstractMetaGraph,
                       Neutral_mut_rate::AbstractFloat, # AbstractFloat?
                       length_ROI::Int, seed::MersenneTwister)

    mutations = []
    λ = Neutral_mut_rate * length_ROI

    for e in edges(Tree)
        Tfinal = get_prop(Tree, dst(e), :Time) # Prima era len = ...
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
function transform_genome(genome::LongDNASeq, pos::Vector{Any},
                                                        seed::MersenneTwister)
    transition_matrix =
        DataFrame(A = [0.0, 0.33333333, 0.33333333, 0.33333333],
                  C = [0.33333333, 0.0, 0.33333333, 0.33333333],
                  G = [0.33333333, 0.33333333, 0.0, 0.33333333],
                  T = [0.33333333, 0.33333333, 0.33333333, 0.0])
    for p in pos
        ## println("pos: ",p)
        nucleotide = genome[p]
        ## println("subsequence: ",genome[p - 2 : p + 2])
        ## println("nuclotide: ",nucleotide)

        ## Come nella simulatione, ho una certa prob che avvenga una
        ## simulatione
        prob_cum = cumsum(transition_matrix[!, string(nucleotide)])

        ## println("prob_cum: ",prob_cum)
        k = round( rand(seed), digits = 7) # with 7 decimals
        target_nucl = collect(k .<= prob_cum)
        ## println("target_nucl: ",target_nucl)
        min = findfirst(target_nucl)
        ## println("min: ",min)
        new_nucleotide = collect(names(transition_matrix)[min])[1] # type char
        ## println("new_nucleotide: ",new_nucleotide)
        genome[p] = DNA(new_nucleotide)
    end
    return genome
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
                             len_ROI::Int)
    ## Create reference genome
    g_seq = randdnaseq(seed, len_ROI)
    rec = FASTA.Record("Reference", g_seq)
    w = FASTA.Writer(open("Reference.fasta", "w"))
    write(w, rec)
    close(w)

    ## Compute mutations ∀ node
    for i in [1,2,3]
        check = false
        n_mutations = evolution_seq(Tree, neural_mut_rate, len_ROI, seed)
        check = check_num_path_mutation(Tree, len_ROI)
        if check == true && i >= 3
            ## If it is an error, signal one!
            return "ERROR: mutations are more than lenght ROI -> (file fasta)."
        else
            break
        end
    end
        ## IMHO manca una 'end' qui!!!!

        possible_position = Set(1:len_ROI)
        # g_seq_e = copy(g_seq) #reference g_seq not change

        position_used = []

        # create fasta ∀ nodes
        for e in edges(Tree)
            g_seq_e = LongDNASeq()
            if has_prop(Tree, src(e), :Fasta)
                # println("il source ha un file Fasta")
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
            g_seq_e = transform_genome(g_seq_e, pos_edge, seed)
            set_prop!(Tree, dst(e), :Fasta, g_seq_e)
        end

        # Return fasta of leaf nodes
        leafs = get_leafs(Tree)
        fasta_samples = []
        for l in leafs
            f = get_prop(Tree, l, :Fasta)
            push!(fasta_samples, f)
        end

        return g_seq, fasta_samples, position_used
end

"""
Creates a input for tool ART -> FASTA file (WITH REF FASTA).
"""
function Molecular_evolution_ISA(Tree::AbstractMetaGraph,
                             neural_mut_rate::Float64,
                             seed::MersenneTwister,
                             path::String)
    ## load reference genome
    g_seq = LongDNASeq()
    open(FASTA.Reader, path) do reader
        for record in reader
            g_seq = FASTX.sequence(record)
            #println("g_seq: ", g_seq)
        end
    end
    #println("g_seq: ", g_seq)
    len_ROI = length(g_seq)

    ## Compute mutations ∀ node
    for i in [1,2,3]
        check = false
        n_mutations = evolution_seq(Tree, neural_mut_rate, len_ROI, seed)
        check = check_num_path_mutation(Tree, len_ROI)
        if check == true && i >= 3
            ## If it is an error, signal one!
            return "ERROR: mutations are more than lenght ROI -> (file fasta)."
        else
            break
        end
    end

        possible_position = Set(1:len_ROI)
        # g_seq_e = copy(g_seq) #reference g_seq not change

        position_used = []

        # create fasta ∀ nodes
        for e in edges(Tree)
            g_seq_e = LongDNASeq()
            if has_prop(Tree, src(e), :Fasta)
                # println("il source ha un file Fasta")
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
            g_seq_e = transform_genome(g_seq_e, pos_edge, seed)
            set_prop!(Tree, dst(e), :Fasta, g_seq_e)
        end

        # Return fasta of leaf nodes
        leafs = get_leafs(Tree)
        fasta_samples = []
        for l in leafs
            f = get_prop(Tree, l, :Fasta)
            push!(fasta_samples, f)
        end

        return g_seq, fasta_samples, position_used
end

"""
    comment gemonic evolution
"""
function genomic_evolution(Seq_f::LongDNASeq,
                           Model_Selector::String, #Ploidity::Int,
                           rate_Indel::AbstractFloat,
                           size_indel::Int,
                           branch_length::AbstractFloat,
                           Model_Selector_matrix::DataFrame,
                           prob_A::Vector{Float64},
                           prob_C::Vector{Float64},
                           prob_G::Vector{Float64},
                           prob_T::Vector{Float64},
                           seed::MersenneTwister)


    sequence_father = copy(Seq_f)
    len_father = length(sequence_father)

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
        λ_indel = (len_father+1) * rate_Indel
        λ_tot = sum([λ_A, λ_C, λ_G, λ_T, λ_indel])

        curr_time = curr_time + rand(seed, Exponential(1 / λ_tot), 1)[1]

        prob_vet = vcat(λ_A, λ_C, λ_G, λ_T, λ_indel) ./ λ_tot
        prob_cum = cumsum(prob_vet)
        k = rand(seed)

        mutation = collect(k .<= prob_cum)
        min_mutation = findfirst(mutation)

        if min_mutation == 5 #indel
            e = rand(seed, ["insertion", "deletion"])
            #choose initial position
            init_pos = rand(seed, 1:length(sequence_father))
            if e == "insertion"
                length_ins = rand(seed, Poisson(size_indel), 1)[1] + 1

                if len_father - init_pos < length_ins
                     length_ins = len_father - init_pos
                 end

                 init_pos_ins = rand(seed, 1:length(sequence_father))
                 insertion_sequence =
                                   sequence_father[init_pos:init_pos+length_ins]
                 new_sequence = sequence_father[1:init_pos_ins]
                 #new_sequence = LongSubSeq(sequence_father, 1:init_pos_ins)
                 append!(new_sequence, insertion_sequence)
                 append!(new_sequence, sequence_father[init_pos_ins + 1 : end])
                 sequence_father = new_sequence
                 len_father = len_father + length_ins
                 #update
                 #n_A += count("A" , string(insertion_sequence))
                 #n_C += count("C" , string(insertion_sequence))
                 #n_G += count("G" , string(insertion_sequence))
                 #n_T += count("T" , string(insertion_sequence))
             else #deletion
                 length_del = rand(seed, Poisson(size_indel), 1)[1] + 1
                 init_pos_del = rand(seed, 1:length(sequence_father))
                 end_pos_del = min(len_father, init_pos_del + length_del)
                 new_sequence = sequence_father[1:init_pos_del]

                 if end_pos_del != len_father
                     append!(new_sequence, sequence_father[end_pos_del + 1:end])
                 end

                 deletion_sequence = sequence_father[init_pos_del:end_pos_del]
                 sequence_father = new_sequence
                 len_father = len_father - length_del
                 #update
                 #n_A -= count("A" , string(deletion_sequence))
                 #n_C -= count("C" , string(deletion_sequence))
                 #n_G -= count("G" , string(deletion_sequence))
                 #n_T -= count("T" , string(deletion_sequence))
            end

        elseif min_mutation == 4 #T
            pos_mutation = rand(seed, Ts)
            prob_cum_T = cumsum(prob_T)
            k = rand(seed)
            mutation = collect(k .<= prob_cum_T)
            ff = findfirst(mutation)
            new_nucleotide = collect(names(Model_Selector_matrix)[ff])[1]
            ## println("new_nucleotide: ",new_nucleotide)
            sequence_father[pos_mutation] = DNA(new_nucleotide)

        elseif min_mutation == 3 #G
            pos_mutation = rand(seed, Gs)
            prob_cum_G = cumsum(prob_G)
            k = rand(seed)
            mutation = collect(k .<= prob_cum_G)
            ff = findfirst(mutation)
            new_nucleotide = collect(names(Model_Selector_matrix)[ff])[1]
            sequence_father[pos_mutation] = DNA(new_nucleotide)

        elseif min_mutation == 2 #C
            pos_mutation = rand(seed, Cs)
            prob_cum_C = cumsum(prob_C)
            k = rand(seed)
            mutation = collect(k .<= prob_cum_C)
            ff = findfirst(mutation)
            new_nucleotide = collect(names(Model_Selector_matrix)[ff])[1]
            sequence_father[pos_mutation] = DNA(new_nucleotide)

        elseif min_mutation == 1 #A
            pos_mutation = rand(seed, As)
            prob_cum_A = cumsum(prob_A)
            k = rand(seed)
            mutation = collect(k .<= prob_cum_A)
            ff = findfirst(mutation)
            new_nucleotide = collect(names(Model_Selector_matrix)[ff])[1]
            sequence_father[pos_mutation] = DNA(new_nucleotide)
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
    return sequence_father
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
                          branch_length::AbstractFloat,
                          seed::MersenneTwister)
    Tree_SC = copy(Tree)
    Ref = LongDNASeq()
    ## load reference genome
    open(FASTA.Reader, path) do reader
        for record in reader
            Ref = FASTX.sequence(record)
        end
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
    ##### fino a qui
    for e in edges(Tree_SC)
        g_seq_e = LongDNASeq()
        println(e)
        if has_prop(Tree_SC, src(e), :Fasta)
            # println("il source ha un file Fasta")
            g_seq_e = copy(get_prop(Tree_SC, src(e), :Fasta))
        else
            g_seq_e = copy(Ref)
        end
        sequence = genomic_evolution(g_seq_e,
                                     Selector,
                                     rate_Indel,
                                     size_indel,
                                     branch_length,
                                     Model_Selector_matrix,
                                     prob_A,
                                     prob_C,
                                     prob_G,
                                     prob_T,
                                     seed)
        set_prop!(Tree_SC, dst(e), :Fasta, sequence)
    end

    # Return fasta of leaf nodes
    leafs = get_leafs(Tree_SC)
    fasta_samples = []
    for l in leafs
        f = get_prop(Tree_SC, l, :Fasta)
        push!(fasta_samples, f)
    end


    return Ref, fasta_samples, Tree_SC
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
                          branch_length::AbstractFloat,
                          seed::MersenneTwister)
    Tree_SC = copy(Tree)

    ## Create reference genome
    Ref = randdnaseq(seed, Len)
    rec = FASTA.Record("Reference", Ref)
    w = FASTA.Writer(open("Reference.fasta", "w"))
    write(w, rec)
    close(w)

    #Model_Selector
    Model_Selector_matrix = Q(Selector, params)

    prob_A = Model_Selector_matrix[:,1] ./ sum(Model_Selector_matrix[:,1])
    prob_C = Model_Selector_matrix[:,2] ./ sum(Model_Selector_matrix[:,2])
    prob_G = Model_Selector_matrix[:,3] ./ sum(Model_Selector_matrix[:,3])
    prob_T = Model_Selector_matrix[:,4] ./ sum(Model_Selector_matrix[:,4])
    ##### fino a qui
    for e in edges(Tree_SC)
        g_seq_e = LongDNASeq()
        println(e)
        if has_prop(Tree_SC, src(e), :Fasta)
            # println("il source ha un file Fasta")
            g_seq_e = copy(get_prop(Tree_SC, src(e), :Fasta))
        else
            g_seq_e = copy(Ref)
        end
        sequence = genomic_evolution(g_seq_e,
                                     Selector,
                                     rate_Indel,
                                     size_indel,
                                     branch_length,
                                     Model_Selector_matrix,
                                     prob_A,
                                     prob_C,
                                     prob_G,
                                     prob_T,
                                     seed)
        set_prop!(Tree_SC, dst(e), :Fasta, sequence)
    end

    # Return fasta of leaf nodes
    leafs = get_leafs(Tree_SC)
    fasta_samples = []
    for l in leafs
        f = get_prop(Tree_SC, l, :Fasta)
        push!(fasta_samples, f)
    end

    return Ref, fasta_samples, Tree_SC
end

    ### end of file -- SingleCellExperiment.jl
