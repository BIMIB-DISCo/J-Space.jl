### -*- Mode: Julia -*-

### sampling_phylogenetic_relation_and_genotype.jl

#module sampling_phylogenetic_relation_and_genotype

using PhyloNetworks        # Create newick format
using BioSequences         # This library is used for create FASTA file
using FASTX                # This library is used for read FASTA file
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
        return "Error: Choose corret Mode -> Random or Neighbourhood "
    end
    label_cells = []

    for c in list_cell
        push!(label_cells, get_prop(G, c, :id))
        push!(list_mut, get_prop(G, c, :mutation))
    end
    return label_cells, list_mut
end


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
    for i in nrow(events_df):-1:1
        id_c = events_df[i, :Cell]
        if id_c ∈ sample_list

            sample_id = findall(x -> x == id_c, sample_list)

            time = events_df[i, :Time]
            father = events_df[i, :Notes][1]
            child = id_c
            event = events_df[i, :Event]

            filter!(s -> s != id_c, sample_list)

            push!(sample_list, father)

            if events_df[i, :Event] == "Mutation"
                new_driver = events_df[i, :Notes][2]
                father_mut = mutations[sample_id][1][1:end-1]

                length(father_mut) == 1 ? push!(mutations, father_mut[1]) :
                                          push!(mutations, father_mut)
            else
                push!(mutations, mutations[sample_id][1])
            end

            subpop_child = findall(m -> m == mutations[end], set_mut)[1]

            deleteat!(mutations, sample_id)

            push!(matrix_R, [father, child, time, subpop_child, event])
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
        return matrix_R, []
    end
end


"""
Create and plot tree
"""
function create_tree(matrix::DataFrame, newick::Bool)

    tree = MetaDiGraph(matrix, :Father, :Child)

    #find root
    t_v = filter_vertices(tree, :Time)
    t_vs = [get_prop(tree, v, :Time) for v in t_v]
    root = findall(t -> t === missing, t_vs)[1]

    if get_prop(tree, root, :Time) === missing
        set_prop!(tree, root, :Time, 0)
        set_prop!(tree, root, :Subpop_Child, 1)
    end

    tree_reduce = reduce_tree(tree, root)

    ## Create format newick
    if newick == true
        println("create format newick....")
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
function reduce_tree(tree::AbstractGraph, root::Int)

    leafs = get_leafs(tree)

    t_r = SimpleDiGraph()
    tree_reduce = MetaDiGraph(t_r)

    map_nodes = []

    for leaf in leafs

        yen_k = yen_k_shortest_paths(tree, root, leaf)
        path = yen_k.paths[1]

        v_mid = get_vertex_middle(tree, path)

        filter!(v -> v ∉ v_mid, path)

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
                end
            end
        end
    end
    return tree_reduce
end

### end of file -- sampling_phylogenetic_relation_and_genotype.jl
