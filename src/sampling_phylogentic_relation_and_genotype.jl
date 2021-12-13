#module sampling_phylogenetic_relation_and_genotype

include("JOG_Space.jl")
using MetaGraphs
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
    for c in list_cell
        push!(label_cells, get_prop(G, c, :id))
        push!(list_mut, get_prop(G, c, :mutation))
    end
    return label_cells, list_mut
end #può essere più efficiente, list_cell potrebe contenere già le labels

#function that create matrix relation of sampling phylogentic relation
function create_matrix_relational(G::AbstractGraph, events_df::DataFrame,
                               sample_list::Vector{Any}, mutations::Vector{Any})

    matrix_R = DataFrame(Father = Any[], Child = Any[], Time = Float64[],
                                         Subpop_Child = Any[], Event = String[])
    set_mut = unique(get_drivermut(G))
    for i in nrow(events_df):-1:1
        #println(i)
        id_c = events_df[i, :Cell]
        if id_c ∈ sample_list
            sample_id = findall(x -> x == id_c, sample_list)
            subpop_child = findall(m->m == mutations[sample_id][1], set_mut)[1]
            time = events_df[i, :Time]
            if events_df[i,:Event] != "Migration"
                father = events_df[i, :Notes][1]
                child = id_c
                event = events_df[i, :Event]
                filter!(s -> s != id_c, sample_list) #tolgo il sample trovato
                push!(sample_list, father) #aggiungo il nuovo sample
                if length(events_df[i, :Event]) == "Mutation"
                    father_mut = filter!(m -> m != events_df[i, :Notes][2],
                                                        mutations[sample_id][1])
                    push!(mutations, father_mut)
                 else
                     push!(mutations, mutations[sample_id][1])
                 end
                deleteat!(mutations, sample_id)
            else
                father = undef
                child = id_c
                event = events_df[i, :Event]
            end
            push!(matrix_R,[father, child, time, subpop_child, event])
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
        f, ax, p = graphplot(tree, layout = Buchheim(), node_color=color, nlabels = [string(v) for v in vertices(tree)])
        hidedecorations!(ax)
        hidespines!(ax)
        save(path, f)
        display(f)
    end
    tree_reduce = reduce_tree(tree)
    return tree_reduce
end

#function main of module
function sampling_phylogentic_relation(G::AbstractGraph, Mode::String,
                                                        df::DataFrame, L::Int;
                                                                    dist::Int=0)
    event_df = copy(df)
    times = df[!, :Time]#df.Time
    tf = times[end]#non so se serve
    event_df = event_df[event_df.Event .!= "Death", :]
    list_sample, list_mut = list_sampling_cell(G, Mode, L, dist = dist)
    matrix_R = create_matrix_relational(G, event_df, list_sample, list_mut)
    matrix_R_without_migration = matrix_R[matrix_R.Event .!= "Migration", :]
    return matrix_R_without_migration
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
    println("leafs: ",leafs)
    t_r = SimpleDiGraph()
    tree_reduce = MetaDiGraph(t_r)
    map_nodes = []
    for leaf in leafs
        println("leaf: ",leaf)
        yen_k = yen_k_shortest_paths(tree, 1, leaf)
        path = yen_k.paths[1]
        v_mid = get_vertex_middle(tree, path)
        println("v_mid: ",v_mid)
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
                    println("ho aggiunto un nodo tra ", source, "e ", dest)
                end
            end
        end
    end
    return tree_reduce
end


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

function format_newick_1(Tree::AbstractMetaGraph, root::Int, dict::Dict,
                                                             leafs::Vector{Any})
    println("dict: ",dict)
    if root != 1
        leafs = check_leafs(Tree, leafs, root)
    end
    deepest_leaf = max_shortest_path(Tree, leafs, root)
    println("deepest_leaf: ",deepest_leaf)
    parent = dict[deepest_leaf][:in]
    println("parent: ",parent)
    child = dict[parent][:out]
    println("child: ",child)
    str = "("*string(child[1])*","*string(child[2])*")"*string(parent)
    println("str -> ",str)
    new_parent = dict[parent][:in]
    println("nuovo parent: ",new_parent)
    while new_parent != 0 && parent != root
        child = dict[new_parent][:out]
        println("child: ",child)
        filter!(v -> v != parent, child)[1]
        println("filtrato child :", child)
        if get(dict[child[1]], :out, 0) == 0
            str = "("*str*","*string(child[1])*")"*string(new_parent)
            println("stringa aggiornata: ", str)
        else
            println("richiamo la funzione con la root: ",child[1])
            str_int = format_newick_1(Tree, child[1], dict, leafs)
            println("str_int: ",str_int)
            str = "("*str*","*str_int*")"*string(new_parent)
            println("stringa aggiornata con pezzo interno: ", str)
        end
        println("new parent: ",new_parent, " e root: ",root)
        #if new_parent == root
        #    println("sono uguali")
        #    return str
        #end
        parent = new_parent
        get(dict[new_parent], :in, 0) == 0 ? new_parent = 0 :
                                             new_parent = dict[new_parent][:in]
        println("nuovo parent: ",new_parent)
        if new_parent == root && length(dict[new_parent][:out]) == 1
            println("sono uguali")
            str = "("*str*")"*string(new_parent)
            new_parent = 0
        end
    end

    return str
end

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
#=
g_meta = spatial_graph(21, 21, dim = 1, n_cell=1)#creo il grafico

df, G, n_cell_alive = simulate_evolution(g_meta, 150.0, 0.2, 0.01,
                                                            0.01, 0.005, 0.4)

#matrix_R = sampling_phylogentic_relation(G, "Random", df, 100)
matrix_R = sampling_phylogentic_relation(G, "Neighbourhood", df, 100, dist = 8)

tree = create_tree(matrix_R, path = "Plot\\tree_relation_original.png")
tree = create_tree(matrix_R)

tree_red = reduce_tree(tree)
=#
#end #module
