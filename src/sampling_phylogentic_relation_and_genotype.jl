include("project.jl")#uso alcune funzioni

#function that return a list of sample
function list_sampling_cell(G::AbstractMetaGraph, Mode::String, L::Int;
                                                                  dist::Int = 0)
    list_cell = []
    list_mut = []
    cs_alive = cells_alive(G)
    #println("cs_alive: ",cs_alive)
    if Mode == "Random"#choose random sample -> L = num sample
        set = collect(1:length(cs_alive))
        #println("set: ",set)
        for i in 1:L
            #println("i: ",i)
            pos = rand(set)
            #println("pos: ",pos)
            push!(list_cell, cs_alive[pos])
            #println("ho aggiornato list_cell: ",list_cell)
            filter!(s -> s != pos, set)
            #println("ho aggiornato- set: ",set)
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
    #println("set_mut: ",set_mut)
    #println("sample_list: ",sample_list)
    #println("mutation: ",mutations)
    #println("\n")
    for i in nrow(events_df):-1:1
        println(i)
        #events_df[i, :Event] == "duplicate" ?  id_c = events_df[i, :Nodes][1] :
        #                                     id_c = events_df[i, :Nodes]
        find_sample = false
        if events_df[i, :Event] == "duplicate"
            if events_df[i, :Nodes][1] ∈ sample_list
                id_c = events_df[i, :Nodes][1]
                find_sample = true
            elseif events_df[i, :Nodes][2] ∈ sample_list
                id_c = events_df[i, :Nodes][2]
                find_sample = true
            end
        elseif events_df[i, :Nodes] ∈ sample_list
            id_c = events_df[i, :Nodes]
            find_sample = true
        end
        #println("Sto guardando questo evento: ", events_df[i, :Event], "\t" ,events_df[i, :Nodes],"\t" , events_df[i, :Notes])
        #println(sample_list,"\nqui c'è? ", id_c)
        #if id_c ∈ sample_list
        if find_sample == true
            println("Ho trovato un sample")
            #trovo id nella mia lista sample
            sample_id = findall(x -> x == id_c, sample_list)
            #println("sample_id: ",sample_id)
            #trovo la sottopopolazione di appartenteza
            println("set_mut: ",set_mut)
            println("mutations[sample_id]: ",mutations[sample_id])
            subpop_child = findall(m->m == mutations[sample_id][1], set_mut)[1]
            println("subpop_child: ",subpop_child)
            time = events_df[i, :Time]
            if events_df[i,:Event] == "duplicate"
                father = events_df[i, :Notes][1]
            #    println("father: ",father)
                child = id_c
            #    println("length(events_df[i, :Notes]): ",length(events_df[i, :Notes]))
                length(events_df[i, :Notes]) == 1 ? event=events_df[i,:Event] :
                                                    event = "Mutation"
            #    println("event: ",event)
                #aggiorno sample_list e mutations
            #    println("\nDA CONTROLLARE\n")
            #    println(sample_list)
                filter!(s -> s != id_c, sample_list)
            #    println("Ho aggiornato sample_list, ho tolto -> ",id_c)
            #    println(sample_list)
                push!(sample_list,father)
            #    println("Ho aggiornato sample_list, ho aggiunto -> ",father)
            #    println(sample_list)
            #    println(events_df[i, :Notes]," ",typeof(events_df[i, :Notes]))
                if length(events_df[i, :Notes]) > 1
                    println("voglio ottenere la mutazione associata alla cellula padre,quindi senza mutazione")
                    println("events_df[i, :Notes]: ",events_df[i, :Notes])
                    println("events_df[i, :Notes][2]: ",events_df[i, :Notes][2])
                    println("mutations[sample_id]: ", mutations[sample_id])
                    println("old mut: ",filter!(m ->m != events_df[i, :Notes][2],mutations[sample_id][1]))
                end
            #    println("length(events_df[i, :Notes]): ",length(events_df[i, :Notes]))
                length(events_df[i, :Notes]) > 1 ? push!(mutations,
                                       filter!(m ->m != events_df[i, :Notes][2],
                                                       mutations[sample_id][1])) :
                                       push!(mutations, mutations[sample_id][1])
                println("Aggiornato mutations, aggiunto: ",mutations)
            #    println("\n\n")
                println(mutations)
                deleteat!(mutations, sample_id)
                println("Ho aggiornato mutations, ho tolto -> ",sample_id)
                println(mutations)
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
                                                            0.01, 0.05, 0.4)

matrix_R = sampling_phylogentic_relation(G,"Random",df,10)
