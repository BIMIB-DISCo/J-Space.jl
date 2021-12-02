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


g_meta = spatial_graph(21, 21, dim = 1, n_cell=1)#creo il grafico

df, G, n_cell_alive = simulate_evolution(g_meta, 150.0, 0.2, 0.01,
                                                            0.01, 0.005, 0.4)

matrix_R = sampling_phylogentic_relation(G, "Random", df, 11)
#matrix_R = sampling_phylogentic_relation(G, "Neighbourhood", df, 100, dist = 8)

tree = plot_tree(matrix_R, "Plot\\tree_relation_original.png")


matrix_R[matrix_R.Father .== get_prop(tree,1,:name),:]#da controllare per vedere
                                                  #se ho una root in comune o no

leafs = []
for v in vertices(tree)
    if indegree(tree,v) == 1 && outdegree(tree,v) == 0
        push!(leafs,v)
    end
end
leafs

color = [:blue for i in 1:nv(tree)]
color[1] = :black
color[leafs] .= :green
f, ax, p = graphplot(tree, layout = Buchheim(), node_color=color)
hidedecorations!(ax)
hidespines!(ax)
save("Plot\\tree_relation_leafs.png", f)
display(f)

path = yen_k_shortest_paths(tree,1,leafs[2])
path.paths[1]

color = [:blue for i in 1:nv(tree)]
color[1] = :black
color[path.paths[1]] .= :yellow
f, ax, p = graphplot(tree, layout = Buchheim(), node_color=color)
hidedecorations!(ax)
hidespines!(ax)
save("Plot\\tree_relation_pathleaf.png", f)
display(f)

function get_vertex_middle(tree::AbstractGraph, path::Vector{Int})
    vertix_middle = []
    for v in path
        if indegree(tree,v) == 1 && outdegree(tree,v) == 1
            push!(vertix_middle,v)
        end
    end
    return vertix_middle
end

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
tree2 = copy(tree)
tree = copy(tree2)

n_delete_nodes = 0
sort(leafs)
for leaf in leafs
    leaf = leaf - n_delete_nodes
    println("\n",leaf)
    yen_k = yen_k_shortest_paths(tree, 1, leaf)
    println("yen_k: ",yen_k)
    path = yen_k.paths[1]
    v_mid = get_vertex_middle(tree, path)
    #n_delete_nodes = n_delete_nodes + sum(leaf .> v_mid)
    # color
    color = [:black for i in 1:nv(tree)]
    color[v_mid] .= :red
    color[path[1]] = :green
    color[leaf] = :blue
    f, ax, p = graphplot(tree, layout = Buchheim(), node_color=color, nlabels=[string(v) for v in vertices(tree)])
    save("Plot\\"*string(leaf)*".png", f)
    #fine color
    println("path: ",path)
    println("v_mid: ",v_mid)
    indexs = [findall(id -> id == v, path)[1] for v in v_mid]
    println("indexs: ",indexs)
    start, finish = get_start_finish(path,indexs)
    println("start: ",start)
    println("finish: ",finish)

    for i in length(start):-1:1
        println("start[i]: ",start[i])
        println("finish[i]: ",finish[i])
        vx_delete = path[start[i]+1:finish[i]-1]
        println("vx_delete: ",vx_delete)
        add_edge!(tree, path[start[i]], path[finish[i]])#ricollego il grafo
        println("ho aggiunto l'edge")
        #elimino vertici dall'albero
        length(vx_delete) > 1 ? rem_vertices!(tree.graph, vx_delete) :
                                rem_vertex!(tree.graph, vx_delete[1])
        println("ho cancellato il/i nodo/i")
    end

    #=if v ∈ v_mid
        println("ho trovato un nodo da eliminare")
        #println("tot nodi rimasti: ",nv(tree))
        start = path[1] #inizio taglio
        println("il mio inizio: ",start)
        child = outneighbors(tree, v)[1]
        println("il figlio: ", child)
        vx_delete = [v] #lista nodi da cancellare
        println("vs_delete: ",vx_delete)
        while child ∈ vertix_middle
            push!(vx_delete,child)#aggiorno
            child = outneighbors(tree,child)[1]
        end
        println("finito while: ",vx_delete)
        finish = child#fine taglio
        println("fine: ",finish)
        color = [:black for i in 1:nv(tree)]
        color[vx_delete] .= :red
        color[start] = :green
        color[finish] = :blue
        f, ax, p = graphplot(tree, layout = Buchheim(), node_color=color)
        save("Plot\\"*string(v)*"_prima.png", f)
        #elimino i vertici dall'albero
        length(vx_delete) > 1 ? rem_vertices!(tree.graph, vx_delete) :
                                rem_vertex!(tree.graph, vx_delete[1])
        println("ho cancellato il/i nodo/i")
        add_edge!(tree, start, finish)#ricollego il grafo
        println("ho aggiunto l'edge")
        #color = [:black for i in 1:nv(tree)]
        color[start] = :green
        color[finish] = :blue
        f, ax, p = graphplot(tree, layout = Buchheim())#, node_color=color)
        save("Plot\\"*string(v)*"_zdopo.png", f)
        #tolgo i valori da vertix_middle
        println("filtro i valori in vertex_middle: ",vertix_middle)
        filter!(v -> v ∉ vx_delete, vertix_middle)
        println("ho tolto -> ",vx_delete)
        println(vertix_middle)
    end
    println("v ∉ vertex_middle -> quindi non faccio niente")
    if continuo == false
        break
    end=#
end

color = [:black for i in 1:nv(tree)]
color[66] = :red
f, ax, p = graphplot(tree, layout = Buchheim(), node_color=color)
end #module
