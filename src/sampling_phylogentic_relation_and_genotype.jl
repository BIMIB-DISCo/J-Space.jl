#module sampling_phylogenetic_relation_and_genotype
include("JOG_Space.jl")
#using MetaGraphs #questo perchè è qui?
using PhyloNetworks #create newick format
using BioSequences #this library is used for create FASTA file
using FASTX #this labrary is used for read FASTA file
using DataStructures #method counter
#function that return a list of sample
function list_sampling_cell(G::AbstractGraph, Mode::String, L::Int;
                                                                  dist::Int = 0)
    list_cell = []
    list_mut = []
    cs_alive = cells_alive(G)
    if Mode == "Random"#choose random sample -> L = num sample
        set = collect(1:length(cs_alive))
        for i in 1:L
            pos = rand(set)
            push!(list_cell, cs_alive[pos])
            filter!(s -> s != pos, set)
        end
    elseif  Mode == "Neighbourhood"#choose neighbourhood -> L = start pos
        poss_list_cell = neighborhood(G, L, dist)
        list_cell = [c for c in poss_list_cell if has_prop(G, c, :id)]
    else
        return "Error: Choose corret Mode -> Random or Neighbourhood "
    end
    label_cells = []
    #println("list_cell: ")
    #println(list_cell)
    #println("\n\n")
    #println("inizio for")
    for c in list_cell
        #println("c -> ", c)
        #println("label -> ", get_prop(G, c, :id))
        push!(label_cells, get_prop(G, c, :id))
        #println("mutation -> ", get_prop(G, c, :mutation))
        push!(list_mut, get_prop(G, c, :mutation))
    end
    return label_cells, list_mut
end #può essere più efficiente, list_cell potrebe contenere già le labels

#function that create matrix relation of sampling phylogentic relation
function create_matrix_relational(G::AbstractGraph, events_df::DataFrame,
                               sample_list::Vector{Any}, mutations::Vector{Any},
                               set_mut::Vector{Any})

    matrix_R = DataFrame(Father = Any[], Child = Any[], Time = Float64[],
                                         Subpop_Child = Any[], Event = String[])
    #set_mut = unique(get_drivermut(G))
    #println("INIZIO")
    #println("sample_list: ")
    #println(sample_list)
    for i in nrow(events_df):-1:1
        #println(i)
        id_c = events_df[i, :Cell]
        if id_c ∈ sample_list
            #println(events_df[i, :])
            #println("id_c -> ", id_c)
            sample_id = findall(x -> x == id_c, sample_list) #trovo pos
            #println("sample_id: ",sample_id)
            #println("mutations: ",mutations)

            time = events_df[i, :Time]
            father = events_df[i, :Notes][1]
            child = id_c
            event = events_df[i, :Event]
            #println("EVENTO : ",event)
            filter!(s -> s != id_c, sample_list) #tolgo il sample trovato
            #println("AGGIORNO sample_list")
            #println("tolgo ",id_c, " e aggiungo -> ",father)
            push!(sample_list, father) #aggiungo il nuovo sample
            #println(sample_list)
            if events_df[i, :Event] == "Mutation"
                new_driver = events_df[i, :Notes][2]
                #father_mut = filter!(m -> m != events_df[i, :Notes][2],
                #                                    mutations[sample_id][1])
                father_mut = mutations[sample_id][1][1:end-1]
                length(father_mut) == 1 ? push!(mutations, father_mut[1]) :
                                          push!(mutations, father_mut)
                #println("mutations[end]: ", mutations[end])
            else
                #println("mutations[sample_id][1]: ",mutations[sample_id][1])
                push!(mutations, mutations[sample_id][1])
            end
            #println("mutations[end]: ", mutations[end])
            subpop_child = findall(m -> m == mutations[end], set_mut)[1]
            #println("subpop_child: ",subpop_child)
            deleteat!(mutations, sample_id)
            #println("AGGIORNATO mutations: ")
            #println(mutations)
            push!(matrix_R, [father, child, time, subpop_child, event])
            #println("aggiungo alla matrice -> ", matrix_R[end,:])
        end
    end
    return matrix_R
end

#Create and plot tree
function create_tree(matrix::DataFrame; path::String = "")
#!!!!!!!!!!!!||!!! DEVE ESSERE UN GRAFO DIRETTO, altrimenti non funziona il plot
    tree = MetaDiGraph(matrix, :Father, :Child)
    if path != ""
        GLMakie.activate!()
        color = [:blue for i in 1:nv(tree)]
        color[1] = :black
        f, ax, p = graphplot(tree, layout = Buchheim(), node_color=color,
                                nlabels = [string(v) for v in vertices(tree)])
        hidedecorations!(ax)
        hidespines!(ax)
        save(path, f)
        display(f)
    end
    tree_reduce = reduce_tree(tree)
    color = [:blue for i in 1:nv(tree_reduce)]
    color[1] = :black
    f, ax, p = graphplot(tree_reduce, layout = Buchheim(), node_color=color,
                        nlabels = [string(v) for v in vertices(tree_reduce)] )
    save(path, f)
    return tree_reduce
end

#function main of module
function sampling_phylogentic_relation(G::AbstractGraph, Mode::String,
                                           df::DataFrame, L::Int,
                                              set_mut::Vector{Any}; dist::Int=0)
    event_df = copy(df)
    times = df[!, :Time]#df.Time
    tf = times[end]#non so se serve
    event_df = event_df[event_df.Event .!= "Death", :]
    event_df = event_df[event_df.Event .!= "Migration", :]
    list_sample, list_mut = list_sampling_cell(G, Mode, L, dist = dist)
    matrix_R = create_matrix_relational(G, event_df, list_sample, list_mut,
                                                                        set_mut)
    tree_mut = create_tree_mutation_driver(set_mut, event_df)
    #matrix_R_without_migration = matrix_R[matrix_R.Event .!= "Migration", :]
    #return matrix_R_without_migration
    return matrix_R, tree_mut
end

#return all leaf of tree
function get_leafs(tree::AbstractGraph)
    leafs = []
    for v in vertices(tree)
        if indegree(tree, v) == 1 && outdegree(tree, v) == 0
            push!(leafs, v)
        end
    end
    return leafs
end

#return vertex to delete
function get_vertex_middle(tree::AbstractGraph, path::Vector{Int})
    vertix_middle = []
    for v in path
        if indegree(tree,v) == 1 && outdegree(tree,v) == 1
            push!(vertix_middle, v)
        end
    end
    return vertix_middle
end

#function that reduce tree
function reduce_tree(tree::AbstractGraph)
    #check = true
    leafs = get_leafs(tree)
    #println("leafs: ",leafs)
    t_r = SimpleDiGraph()
    tree_reduce = MetaDiGraph(t_r)
    map_nodes = []
    for leaf in leafs
    #    println("leaf: ",leaf)
        yen_k = yen_k_shortest_paths(tree, 1, leaf)
        path = yen_k.paths[1]
        v_mid = get_vertex_middle(tree, path)
    #    println("v_mid: ",v_mid)
        filter!(v -> v ∉ v_mid, path) #tolgo i nodi da eliminare
        for i in 1:length(path)
            if path[i] ∉ map_nodes
                add_vertex!(tree_reduce)
                property = props(tree, path[i])
                set_props!(tree_reduce, vertices(tree_reduce)[end], property)
                push!(map_nodes, path[i])
                if i != 1
                    source = findall(s -> s == path[i-1], map_nodes)[1]
                    dest = findall(d -> d == path[i], map_nodes)[1]
                    add_edge!(tree_reduce, source, dest)
    #                println("ho aggiunto un nodo tra ", source, "e ", dest)
                end
            end
        end
    end
    #set_prop!(tree_reduce, 1, :Time, 0.0)
    return tree_reduce
end

###################################################### Newick mutations drive

#function that return a tree of mutations driver in newick format
function create_tree_mutation_driver(set_mut::Vector{Any}, event_df::DataFrame)
    #initialize graphs
    t_r = SimpleDiGraph()
    tree = MetaDiGraph(t_r)
    add_vertex!(tree)
    set_prop!(tree, 1, :Time, 0)
    set_prop!(tree, 1, :Mutation, 1)
    #times
    df_mutations = event_df[event_df.Event .== "Mutation", :]
    times = df_mutations[!, :Time]
    #set_mut[3]
    for i in 2:length(set_mut)
        muts = set_mut[i]
        if length(muts) == 2 #sicuramente figlio di della prima => 1
            add_vertex!(tree)
            set_prop!(tree, i, :Time, times[i-1])
            set_prop!(tree, i, :Mutation, muts[2])
            add_edge!(tree, 1, i)
        else
            len = length(muts)
            #filtro tutte le muts in set_mut con una lunghezza -> len-1
            possible_father = [m for m in set_mut[1:i] if length(m) == len-1]
            #tra i risultati cerco il modo di verificare un array dentro un array
            father_id = findall(x -> x ⊆ muts, possible_father)[1]
            father = findall(x -> x == possible_father[father_id], set_mut)[1]
            #trovato faccio il solito "gioco" di aggiungere vertice e le props
            add_vertex!(tree)
            set_prop!(tree, i, :Time, times[i-1])
            set_prop!(tree, i, :Mutation, muts[end])
            #aggiungo l'edge
            add_edge!(tree, father, i)
        end
    end
    #richiamo format newick per trasformarlo come tale
    return tree
    #nel mio sampling ho una prob che undica la sottopopolazione, posso usare
    #quello se devo usare solo quei sample
end

###################################################################### ART INPUT

#function that search max shortert path in tree with a root(possibile ̸= 1)
function max_shortest_path(Tree::AbstractMetaGraph, leafs::Vector{Any},
                                                                      root::Int)
    max = 1
    len = 1
    for l in leafs
        yen_k = yen_k_shortest_paths(Tree, root, l)
        path = yen_k.paths[1]
        if len < length(path)
            max = l
            len = length(path)
        end
    end
    return max
end

#return all leaf in tree or subtree with several root
function check_leafs(Tree::AbstractMetaGraph, leafs::Vector{Any}, root::Int)
    new_leafs = []
    for l in leafs
        yen_k = yen_k_shortest_paths(Tree, root, l)
        if yen_k.paths != []
            push!(new_leafs, l)
        end
    end
    return new_leafs
end

#function that trasform from MetaGraph to newick format => String
function create_newick(Tree::AbstractMetaGraph, root::Int64, dict::Dict,
                                                             leafs::Vector{Any})
    if root != 1
        leafs = check_leafs(Tree, leafs, root)
    end
    deepest_leaf = max_shortest_path(Tree, leafs, root)
    parent = dict[deepest_leaf][:in]
    child = dict[parent][:out]
    t_p = get_prop(Tree, parent, :Time)
    t_1 = get_prop(Tree, child[1], :Time) - t_p
    t_2 = get_prop(Tree, child[2], :Time) - t_p
    str = "(" * string(child[1]) * ":" * string(t_1) * "," * string(child[2]) *
                                   ":" * string(t_2) * ")" * string(parent)
    new_parent = dict[parent][:in]
    t_np = t_p - get_prop(Tree, new_parent, :Time)
    str = str * ":" * string(t_np)
    while new_parent != 0 && parent != root
        child = dict[new_parent][:out]
        filter!(v -> v != parent, child)[1]
        t_1 = get_prop(Tree, child[1], :Time)
        t_1 = t_1 - get_prop(Tree, new_parent, :Time)
        if get(dict[child[1]], :out, 0) == 0
            str = "(" * str * "," * string(child[1]) * ":" * string(t_1) *
                  ")" * string(new_parent)
        else
            str_int = create_newick(Tree, child[1], dict, leafs)
            str = "(" * str * "," * str_int * ")" *string(new_parent)
        end
        parent = new_parent
        get(dict[new_parent], :in, 0) == 0 ? new_parent = 0 :
                                             new_parent = dict[new_parent][:in]
        if new_parent == root && length(dict[new_parent][:out]) == 1
            t_p = get_prop(Tree, parent, :Time)
            str = "(" * str * ":" * string(t_p) * ")" * string(new_parent)
            new_parent = 0
        end
        if new_parent != 0
            t_p = get_prop(Tree, parent, :Time)
            if typeof(get_prop(Tree, new_parent, :Time)) != Missing
                t_np = t_p - get_prop(Tree, new_parent, :Time)
                str = str * ":" * string(t_np)
            else
                str = str * ":" * string(t_p)
            end
        end
    end
    return str
end

#create a dict ∀ vertex with "in" and "out" vertex linked
function dict_vertex(Tree::AbstractMetaGraph)
    dict = Dict()
    es = edges(Tree)
    for e in es
        s = src(e)
        d = dst(e)
        #add entry :out
        check_source = get(dict, s, 0)
        if check_source != 0
            if get(get(dict, s, 0), :out, 0) != 0
                push!(dict[s][:out],d)
            else
                dict[s] = Dict(:in => dict[s][:in], :out => [d])
            end
        else
            dict[s] = Dict(:out => [d])
        end
        #add entry :in
        check_dest = get(dict, d, 0)
        if check_dest == 0
            dict[d] = Dict(:in => s)
        end
    end
    return dict
end

#function main for create newick format -> recursive
function format_newick(Tree::AbstractMetaGraph)
    dict = dict_vertex(Tree)
    leafs = get_leafs(Tree)
    s = create_newick(Tree, 1, dict, leafs)*";"
    net = readTopology(s)
    return net
end

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
function create_input_ART(Tree::AbstractMetaGraph, neural_mut_rate::Float64;
                                         path::String = "", len_ROI::Int = 6000,
                                         single_cell::Bool=true, bulk::Bool=true)
    #create format newick
    net = format_newick(Tree)

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
    return net, g_seq, fasta_samples, position_used
end

####################### CALL ART
#function call ART but only tecnology illumina
function call_ART(len_read::Int, mode::String; path = "")
    cd("Fasta output\\") #cambio directory
    for file in readdir() #scorro tutti i file
        f = hcat(split.(file, ".")...)[1, :]
        if length(f)>1 && f[2] == "fasta"
            mkpath(f[1])
            cd(f[1])
            command = `art_illumina -sam -i $file -l $len_read -ss HS25
                                                                -f 10 -o $mode`
            #run(command)
            cd("..\\")
        end
    end
    cd("..\\")
end

########################## Experiments
#function calculate bulk experiment
function create_bulk_groundtruth(G_seq::LongDNASeq, Fasta_sample::Vector{Any},
                                                          Position::Vector{Any})
    #Dataframe -> POS, REF_nucleotide, %A, %C, %G, %T
    #=df = DataFrame(Position = Int[], Reference = DNA[],  A = Any[], C = Any[],
                                     G = Any[], T = Any[])
    positions = []
    n_sample = length(Fasta_sample)
    for i = 1:len_ROI
        nucl_ref = G_seq[i]
        nucleotides = [f[i] for f in Fasta_sample]
        nucl_counter = counter(nucleotides)
        if nucl_counter[nucl_ref] < n_sample
            pos = i
            k_dict = keys(nucl_counter.map)
            values = []
            for n in names(df)[2:end]
                if n != string(nucl_ref)
                    n_dna = DNA(collect(n)[1])
                    push!(values, (nucl_counter[n_dna]*100)/n_sample)
                else
                    push!(values, "ref")
                end
            end
            push!(df, [pos, nucl_ref, values[1], values[2], values[3], values[4]])
        end
    end=#
    df = DataFrame(MUT = Any[], VAF = AbstractFloat[])
    n_sample = length(Fasta_sample)
    Position_sort = sort(Position)
    for pos in Position_sort
        nucl_ref = G_seq[pos]
        nucleotides = [f[pos] for f in Fasta_sample]
        nucl_counter = counter(nucleotides)
        k_dict = keys(nucl_counter.map)
        for k in k_dict
            #println("k :",k,"\t",typeof(k))
            #println("k != nucl_ref -> ",k != nucl_ref)
            if k != nucl_ref
                mut = string(pos) * "_" * string(nucl_ref) * "->" * string(k)
                vaf = nucl_counter[k]/n_sample
                #println("vaf -> ",vaf)
                push!(df, [mut, vaf])
            end
        end
    end
    return df
end

function bulk_with_noise(Coverage, VAF, n_sample, FP, FN, Positions)
    df_noise = DataFrame(MUT = Any[], VAF = AbstractFloat[])
    for i in 1:length(Positions)
        n_reads = rand(Poisson(Coverage), 1)[1]
        println("n_reads: ", n_reads)
        n_read_muts = rand(Binomial(n_reads, VAF[i, :VAF]), 1)[1]
        println("n_read_muts: ", n_read_muts)
        #n_FP = rand(Binomial(n_sample-n_read_muts, FP), 1)[1]
        n_FP = rand(Binomial(n_reads-n_read_muts, FP), 1)[1]
        println("n_FP: ",n_FP)
        n_FN = rand(Binomial(n_read_muts, FN), 1)[1]
        println("n_FN: ",n_FN)
        n_VAF_Noise = (n_read_muts + n_FP - n_FN)/n_sample
        println("n_VAF_Noise: ", n_VAF_Noise)
        push!(df_noise, [VAF[i, :MUT], n_VAF_Noise])
    end
    return df_noise
end

function experiment_bulk(reference::LongDNASeq, fasta_samples::Vector{Any},
                            position_used::Vector{Any}; Noise::Bool=false,
                                coverage::AbstractFloat=0, FP::AbstractFloat=0,
                                                            FN::AbstractFloat=0)
    #calculate bulk experiment if bulk is true
    VAF_GT = create_bulk_groundtruth(reference, fasta_samples, position_used)
    CSV.write("VAF_GT.vcf", VAF_GT, delim=",")
    if Noise
        if coverage == 0 && FP == 0 && FN == 0
            return "Error: insert covarege, FP and FN"
        end
        n_sample = length(fasta_samples)
        VAF_GT_Noise = bulk_with_noise(coverage, VAF_GT, n_sample, FP, FN,
                                                                  position_used)
        CSV.write("VAF_GT_Noise.vcf", VAF_GT_Noise, delim=",")
    end
    return VAF_GT, VAF_GT_Noise#, g_seq
end
