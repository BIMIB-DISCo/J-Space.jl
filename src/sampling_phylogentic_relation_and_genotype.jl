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
            #trovo id nella mia lista sample
            #println("\nHO TROVATO IL SAMPLE")
            sample_id = findall(x -> x == id_c, sample_list)
            #println("id_c: ",id_c)
            #println("sample_id: ", sample_id)
            #println("sample_list: ",sample_list)
            #println("set_mut: ",set_mut)
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

function sampling_phylogentic_relation(G::AbstractGraph, Mode::String,
                                                        df::DataFrame, L::Int;
                                                                    dist::Int=0)
    event_df = copy(df)
    times = df[!, :Time]#df.Time
    tf = times[end]#non so se serve
    event_df = event_df[event_df.Event .!= "Death", :]
    list_sample, list_mut = list_sampling_cell(G, Mode, L, dist=dist)
    matrix_R = create_matrix_relational(G, event_df, list_sample, list_mut)
    return matrix_R
end


g_meta = spatial_graph(21, 21, dim = 1, n_cell=1)#creo il grafico

df, G, n_cell_alive = simulate_evolution(g_meta, 150.0, 0.2, 0.01,
                                                            0.01, 0.005, 0.4)

matrix_R = sampling_phylogentic_relation(G,"Random",df,100)

matrix_R_without_migration = matrix_R[matrix_R.Event .!= "Migration", :]

#!!!!!!!!!!!!!!!!!!!!! DEVE ESSERE UN GRAFO DIRETTO, altrimenti non funziona il
#plot
tree = MetaDiGraph(matrix_R_without_migration, :Father, :Child)
layout = Buchheim()
f, ax, p = graphplot(tree, layout=Buchheim())
end #module
