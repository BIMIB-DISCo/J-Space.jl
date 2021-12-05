module sampling_phylogenetic_relation_and_genotype
include("JOG_Space.jl")#uso alcune funzioni
include("DataFrameGraphBridge.jl")

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
                push!(sample_list,father) #aggiungo il nuovo sample
                if length(events_df[i, :Event]) == "Mutation"
                    father_mut = filter!(m ->m != events_df[i, :Notes][2],
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
function plot_tree(matrix::DataFrame, name::String)
#!!!!!!!!!!!!||!!! DEVE ESSERE UN GRAFO DIRETTO, altrimenti non funziona il plot
    tree = MetaDiGraph(matrix, :Father, :Child)
    GLMakie.activate!()
    layout = Buchheim()
    color = [:blue for i in 1:nv(tree)]
    color[1] = :black
    f, ax, p = graphplot(tree, layout = Buchheim(), node_color=color)
    hidedecorations!(ax)
    hidespines!(ax)
    save(name, f)
    display(f)
    return tree
end

#function main of module
function sampling_phylogentic_relation(G::AbstractGraph, Mode::String,
                                                        df::DataFrame, L::Int;
                                                                    dist::Int=0)
    event_df = copy(df)
    times = df[!, :Time]#df.Time
    tf = times[end]#non so se serve
    event_df = event_df[event_df.Event .!= "Death", :]
    list_sample, list_mut = list_sampling_cell(G, Mode, L, dist=dist)
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
            push!(vertix_middle,v)
        end
    end
    return vertix_middle
end

#return index path for delete nodes/vertex
function get_start_finish(path::Vector{Int},indexs::Vector{Int})
    start = []
    finish = []
    push!(start,indexs[1]-1)
    for i in 2:length(indexs)
        #println("indexs[i]: ", indexs[i])
        #println("indexs[i-1]: ", indexs[i-1])
        if indexs[i] == indexs[i-1] + 2
            #println("è la fine di una serie")
            push!(finish, indexs[i] - 1)
            #println("aggiorno finish con ->", indexs[i]-1)
            #println("Finish: ", finish)
            #println("i, length(indexs): ",i,"\t",length(indexs))
            if i <= length(indexs)
                #println("non sono arrivato alla fine")
                push!(start, indexs[i] - 1)
                #println("ho agggiornato lo start con -> ",indexs[i]-1)
                #println(start)
            end
        elseif indexs[i] > indexs[i-1] + 2
            #println("ho un salto maggiore, 2 o + nodi consecutivi che non devo cancellare")
            push!(finish, indexs[i-1] + 1)
            #println("aggiorno finish con ->",indexs[i-1] + 1)
            #println("Finish: ", finish)
            push!(start, indexs[i] - 1)
            #println("ho agggiornato lo start con -> ",indexs[i]-1)
            #println(start)
        end
    end
    push!(finish,indexs[end]+1)

    return start, finish
end

#function that return leaf with vertices to be eliminated
function control_leaf(tree::AbstractGraph, leafs::Vector{Any})
    leaf = 0
    path = []
    v_mid = []
    for l in leafs
        yen_k = yen_k_shortest_paths(tree, 1, l)
        path = yen_k.paths[1]
        v_mid = get_vertex_middle(tree, path)
        if length(v_mid) != 0
            leaf = l
            break
        end
    end
    return leaf, path, v_mid
end

#function that reduce tree
function reduce_tree(tree::AbstractGraph)
    check = true
    nplot = 1 #da cancellare
    #for i in 1:length(leafs)
    while check
        leafs = get_leafs(tree)
        leaf, path, v_mid = control_leaf(tree, leafs)
        if leaf == 0
            break
        end
        println("ora guardo: ",leaf)
        println("path: ",path)
        println("v_mid: ",v_mid)
        # color
        color = [:black for i in 1:nv(tree)]
        color[v_mid] .= :red
        color[path[1]] = :green
        color[leaf] = :blue
        f, ax, p = graphplot(tree, layout = Buchheim(), node_color=color,
        nlabels=[string(v) for v in vertices(tree)])
        save("Plot\\"*string(nplot)*".png", f)
        #fine color
        indexs = [findall(id -> id == v, path)[1] for v in v_mid]
        println("indexs: ",indexs)
        println("inizia il while")
        n = 0
        while length(indexs) != 0
            start, finish = get_start_finish(path,indexs)
            println("start[i]: ",start[end])
            println("finish[i]: ",finish[end])
            vx_delete = path[start[end]+1:finish[end]-1]
            println("vx_delete: ",vx_delete)
            #println("num vertex: ",nv(tree))
            add_edge!(tree, path[start[end]], path[finish[end]])#ricollego il grafo
            println("ho aggiunto l'edge tra :",path[start[end]]," e ",path[finish[end]])
            new_leaf = false
            #if leaf in nv(tree)-length(vx_delete):nv(tree)
            if leaf in (nv(tree)-length(vx_delete)+1):nv(tree)
                println("n vertex on tree : ", nv(tree))
                println("leaf: ", leaf)
                #nuova leaf: num vertici - foglia - vertici maggiori della foglia
                # -> perchè la delete dei nodi funziona cosi, la label dell'ultimo
                #nodo
                println("leaf .< vx_delete -> ", leaf .< vx_delete )
                println("sum(leaf .< vx_delete) + 1", sum(leaf .< vx_delete))
                id = nv(tree) - leaf - sum(leaf .< vx_delete) + 1
                new_leaf = true
                println("ho cambiato leaf")
                println("id :", id)
            end
            #elimino vertici dall'albero
            rem_vertices!(tree.graph, vx_delete)
            println("ho cancellato il/i nodo/i")
            sort!(vx_delete)
            println("vx_delete ordinato ->", vx_delete)
            new_leaf ? yen_k = yen_k_shortest_paths(tree, 1, vx_delete[id]) :
                       yen_k = yen_k_shortest_paths(tree, 1, leaf)
            path = yen_k.paths[1]
            println("Nuovo path: ",path)
            leaf = path[end]
            v_mid = get_vertex_middle(tree, path)
            println("v_mid: ",v_mid)
            indexs = [convert(Int64,findall(id -> id == v, path)[1]) for v in v_mid]
            println("indexs: ",indexs)
            color = [:black for i in 1:nv(tree)]
            color[v_mid] .= :red
            color[path[1]] = :green
            color[leaf] = :blue
            f, ax, p = graphplot(tree, layout = Buchheim(), node_color=color,
                                    nlabels=[string(v) for v in vertices(tree)])
            println("salvo il plot ->",nplot,"_",n)
            save("Plot\\"*string(nplot)*"_"*string(n)*".png", f)
            n = n + 1
        end
        nplot = nplot + 1
    end
    color = [:black for i in 1:nv(tree)]
    color[1] = :red
    f, ax, p = graphplot(tree, layout = Buchheim(), node_color=color)
    save("Plot\\fine1.png", f)
    return tree
end

g_meta = spatial_graph(21, 21, dim = 1, n_cell=1)#creo il grafico

df, G, n_cell_alive = simulate_evolution(g_meta, 150.0, 0.2, 0.01,
                                                            0.01, 0.005, 0.4)

matrix_R = sampling_phylogentic_relation(G, "Random", df, 100)
#matrix_R = sampling_phylogentic_relation(G, "Neighbourhood", df, 100, dist = 8)

tree = plot_tree(matrix_R, "Plot\\tree_relation_original.png")


matrix_R[matrix_R.Father .== get_prop(tree,1,:name),:]#da controllare per vedere
                                                  #se ho una root in comune o no


tree_backup = copy(tree)
tree = copy(tree_backup)


tree_red = reduce_tree(tree)

end #module
