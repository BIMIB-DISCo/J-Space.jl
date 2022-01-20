### -*- Mode: Julia -*-

### Format_tree.jl

#module FormatTree

#export create_tree_mutation_driver
#export format_newick

#function that return a tree of mutations driver in newick format

"""
Returns a driver mutation tree in Newick format.
"""
function create_tree_mutation_driver(set_mut::Vector{Any}, event_df::DataFrame)
    ## Initialize graphs
    t_r = SimpleDiGraph()
    tree = MetaDiGraph(t_r)
    add_vertex!(tree)
    set_prop!(tree, 1, :Time, 0)
    set_prop!(tree, 1, :Mutation, 1)
    
    ## Times
    df_mutations = event_df[event_df.Event .== "Mutation", :]
    times = df_mutations[!, :Time]

    # set_mut[3]
    for i in 2:length(set_mut)
        muts = set_mut[i]
        if length(muts) == 2   # sicuramente figlio di della prima => 1
            add_vertex!(tree)
            set_prop!(tree, i, :Time, times[i - 1])
            set_prop!(tree, i, :Mutation, muts[2])
            add_edge!(tree, 1, i)
        else
            len = length(muts)
            
            ## Filtro tutte le muts in set_mut con una
            ## lunghezza -> len-1
            
            possible_father =
                [m for m in set_mut[1:i] if length(m) == len - 1]
            
            ## Tra i risultati cerco il modo di verificare un array
            ## dentro un array
            father_id = findall(x -> x ⊆ muts, possible_father)[1]
            father = findall(x -> x == possible_father[father_id], set_mut)[1]

            ## trovato faccio il solito "gioco" di aggiungere vertice
            ## e le props
            
            add_vertex!(tree)
            set_prop!(tree, i, :Time, times[i-1])
            set_prop!(tree, i, :Mutation, muts[end])

            ## aggiungo l'edge
            add_edge!(tree, father, i)
        end
    end
    ## Richiamo format newick per trasformarlo come tale
    return tree
end


"""
Searches max shortert path in tree with a root(possibile ̸= 1)
"""
function max_shortest_path(Tree::AbstractMetaGraph,
                           leafs::Vector{Any},
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


"""
Returns all leaf in tree or subtree with several root.
"""
function check_leafs(Tree::AbstractMetaGraph,
                     leafs::Vector{Any},
                     root::Int)
    new_leafs = []
    for l in leafs
        yen_k = yen_k_shortest_paths(Tree, root, l)
        if yen_k.paths != []
            push!(new_leafs, l)
        end
    end
    return new_leafs
end


"""
Trasforms from MetaGraph to newick format => String
"""
function create_newick(Tree::AbstractMetaGraph,
                       root::Int64,
                       dict::Dict,
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
    str =
        "(" *
        string(child[1]) *
        ":" *
        string(t_1) *
        "," *
        string(child[2]) *
        ":" *
        string(t_2) *
        ")" *
        string(parent)
    
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

        ## Isn't this an 'if'?
        
        get(dict[new_parent], :in, 0) == 0 ?
            new_parent = 0 :
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


"""
Creates a dict ∀ vertex with "in" and "out" vertex linked.
"""
function dict_vertex(Tree::AbstractMetaGraph)
    dict = Dict()
    es = edges(Tree)
    for e in es
        s = src(e)
        d = dst(e)
        ## add entry :out
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
        ## add entry :in
        check_dest = get(dict, d, 0)
        if check_dest == 0
            dict[d] = Dict(:in => s)
        end
    end
    return dict
end


"""
Main function for creating newick format -> recursive
"""
function format_newick(Tree::AbstractMetaGraph)
    dict = dict_vertex(Tree)
    leafs = get_leafs(Tree)
    s = create_newick(Tree, 1, dict, leafs)*";"
    net = readTopology(s)
    return net
end

#end #module

### end of file -- Format_tree.jl

