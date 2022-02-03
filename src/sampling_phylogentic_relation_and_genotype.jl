### -*- Mode: Julia -*-

### sampling_phylogenetic_relation_and_genotype.jl

#module sampling_phylogenetic_relation_and_genotype

using PhyloNetworks        # Create newick format
using BioSequences         # This library is used for create FASTA file
using FASTX                # This library is used for read FASTA file
                           # Notes: why is this needed?  Doesn't
                           # BioSequences suffices?
using DataStructures       # Method counter

"""
Returns a list of samples.
"""
function list_sampling_cell(G::AbstractGraph,
                            Mode::String,
                            L::Int,
                            seed::MersenneTwister;
                            dist::Int = 0)
    list_cell = []
    list_mut = []
    cs_alive = cells_alive(G)
    if Mode == "Random"        # Choose random sample -> L = num sample
        set = collect(1:length(cs_alive))
        for i in 1:L
            pos = rand(seed, set)
            push!(list_cell, cs_alive[pos])
            filter!(s -> s != pos, set)
        end
    elseif  Mode == "Neighbourhood" # Choose neighbourhood -> L = start pos
        poss_list_cell = neighborhood(G, L, dist)
        list_cell = [c for c in poss_list_cell if has_prop(G, c, :id)]
    else
        ## AAAAAAAAARRRRRRGGGGGGGGGGH!!!!!!
        return "Error: Choose corret Mode -> Random or Neighbourhood "
    end
    label_cells = []
    ## println("list_cell: ")
    ## println(list_cell)
    ## println("\n\n")
    ## println("inizio for")

    for c in list_cell
        ## println("c -> ", c)
        ## println("label -> ", get_prop(G, c, :id))
        push!(label_cells, get_prop(G, c, :id))
        ## println("mutation -> ", get_prop(G, c, :mutation))
        push!(list_mut, get_prop(G, c, :mutation))
    end
    return label_cells, list_mut
end # Può essere più efficiente, list_cell potrebe contenere già le
    # labels


### Function that create matrix relation of sampling phylogentic
### relation
"""
Creates the matrix of sampling phylogenetic relations.
"""
function create_matrix_relational(G::AbstractGraph,
                                  events_df::DataFrame,
                                  sample_list::Vector{Any},
                                  mutations::Vector{Any},
                                  set_mut::Vector{Any})

    matrix_R = DataFrame(Father = Any[],
                         Child = Any[],
                         Time = Float64[],
                         Subpop_Child = Any[],
                         Event = String[])
    ## set_mut = unique(get_drivermut(G))
    ## println("INIZIO")
    ## println("sample_list: ")
    ## println(sample_list)
    for i in nrow(events_df):-1:1
        ## println(i)
        id_c = events_df[i, :Cell]
        if id_c ∈ sample_list

            ## println(events_df[i, :])
            ## println("id_c -> ", id_c)

            sample_id = findall(x -> x == id_c, sample_list) #trovo pos

            ## println("sample_id: ",sample_id)
            ## println("mutations: ",mutations)

            time = events_df[i, :Time]
            father = events_df[i, :Notes][1]
            child = id_c
            event = events_df[i, :Event]

            ## println("EVENTO : ",event)

            filter!(s -> s != id_c, sample_list) # Tolgo il sample trovato

            ## println("AGGIORNO sample_list")
            ## println("tolgo ",id_c, " e aggiungo -> ",father)
            push!(sample_list, father) # Aggiungo il nuovo sample

            ## println(sample_list)

            if events_df[i, :Event] == "Mutation"
                new_driver = events_df[i, :Notes][2]
                ## father_mut = filter!(m -> m != events_df[i, :Notes][2],
                ##                      mutations[sample_id][1])
                father_mut = mutations[sample_id][1][1:end-1]

                ## Isn't this an 'if'?
                length(father_mut) == 1 ? push!(mutations, father_mut[1]) :
                                          push!(mutations, father_mut)
                ## println("mutations[end]: ", mutations[end])
            else
                ## println("mutations[sample_id][1]: ",mutations[sample_id][1])
                push!(mutations, mutations[sample_id][1])
            end
            ## println("mutations[end]: ", mutations[end])
            subpop_child = findall(m -> m == mutations[end], set_mut)[1]
            ## println("subpop_child: ",subpop_child)
            deleteat!(mutations, sample_id)
            ## println("AGGIORNATO mutations: ")
            ## println(mutations)
            push!(matrix_R, [father, child, time, subpop_child, event])
            ## println("aggiungo alla matrice -> ", matrix_R[end,:])
        end
    end
    return matrix_R
end


"""
Module entry point.
"""
function sampling_phylogentic_relation(G::AbstractGraph,
                                       Mode::String,
                                       df::DataFrame,
                                       L::Int,
                                       set_mut::Vector{Any},
                                       seed::MersenneTwister,
                                       driver_tree::Int;
                                       dist::Int = 0)
    event_df = copy(df)
    times = df[!, :Time]        # df.Time
    tf = times[end]             # Non so se serve
    event_df = event_df[event_df.Event .!= "Death", :]
    event_df = event_df[event_df.Event .!= "Migration", :]
    list_sample, list_mut = list_sampling_cell(G, Mode, L, seed, dist = dist)
    matrix_R =
        create_matrix_relational(G,
                                 event_df,
                                 list_sample,
                                 list_mut,
                                 set_mut)
    if driver_tree == 1
        tree_mut = create_tree_mutation_driver(set_mut, event_df)
        return matrix_R, tree_mut
    else
        return matrix_R
    end
end


"""
Create and plot tree
"""
function create_tree(matrix::DataFrame, newick::Bool; path::String = "")
    # !!!!!!!!!!!!||!!! DEVE ESSERE UN GRAFO DIRETTO, altrimenti non
    # funziona il plot
    tree = MetaDiGraph(matrix, :Father, :Child)
    #=if path != ""
        GLMakie.activate!()
        color = [:blue for i in 1:nv(tree)]
        color[1] = :black
        f, ax, p =
            graphplot(tree,
                      layout = Buchheim(),
                      node_color = color,
                      nlabels = [string(v) for v in vertices(tree)])
        hidedecorations!(ax)
        hidespines!(ax)
        save(path, f)
        display(f)
    end=#

    tree_reduce = reduce_tree(tree)
    ## color = [:blue for i in 1:nv(tree_reduce)]
    ## color[1] = :black
    # f, ax, p = graphplot(tree_reduce, layout = Buchheim(), node_color=color,

    nlabels = [string(v) for v in vertices(tree_reduce)]

    ## save(path, f)

    ## Create format newick
    if newick == true
        net = format_newick(tree_reduce)
        return tree_reduce, net
    end
    return tree_reduce
end


"""
Return all leaves of tree.
"""
function get_leafs(tree::AbstractGraph)
    leafs = []
    for v in vertices(tree)
        if indegree(tree, v) == 1 && outdegree(tree, v) == 0
            push!(leafs, v)
        end
    end
    return leafs
end


"""
Return vertex to delete.
"""
function get_vertex_middle(tree::AbstractGraph, path::Vector{Int})
    vertix_middle = []
    for v in path
        if indegree(tree,v) == 1 && outdegree(tree,v) == 1
            push!(vertix_middle, v)
        end
    end
    return vertix_middle
end


"""
Reduce tree
"""
function reduce_tree(tree::AbstractGraph)
    ## Check = true
    leafs = get_leafs(tree)
    ## println("leafs: ",leafs)
    t_r = SimpleDiGraph()
    tree_reduce = MetaDiGraph(t_r)
    map_nodes = []
    for leaf in leafs
        ## println("leaf: ",leaf)
        yen_k = yen_k_shortest_paths(tree, 1, leaf)
        path = yen_k.paths[1]
        v_mid = get_vertex_middle(tree, path)
        ## println("v_mid: ",v_mid)
        filter!(v -> v ∉ v_mid, path) # Tolgo i nodi da eliminare
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
                    ## println("ho aggiunto un nodo tra ", source, "e ", dest)
                end
            end
        end
    end
    ## set_prop!(tree_reduce, 1, :Time, 0.0)
    return tree_reduce
end

### end of file -- sampling_phylogenetic_relation_and_genotype.jl
