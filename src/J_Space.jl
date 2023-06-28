### -*- Mode: Julia -*-

### J_Space.jl

module J_Space
using Plots
using MetaGraphs
using Graphs
using NetworkLayout
using GraphMakie
using CairoMakie
using DataFrames
using Random
using UUIDs
using StatsBase
using Distributions
using CSV, Tables
using TOML
using DelimitedFiles
#using Distributed
using LinearAlgebra
### function to exports
export
    ## Simulation
    spatial_graph, simulate_evolution, Start_J_Space,
    ## Sampling
    sampling_phylogentic_relation, create_tree,
    ## Experiment
    experiment_ISA, experiment_bulk, experiment_noISA, experiment_noISA_sign,
    save_Fasta, Q,
    ## ART
    call_ART,
    ## plot
    plot_lattice, plot_tree


"""
Creates a lattice with parameters row and col(optional dim), and it
calls function for initialize lattice with p as probability that
there are a cell have a driver mutation
"""
#### create graph
function spatial_graph(row::Int,
                       col::Int,
                       seed::MersenneTwister;
                       dim::Int = 2,
                       n_cell::Int = 1)
    ###3D
    if dim == 3
        d = convert(Int, trunc((row * col) ^ (1/3))) + 1
        G = Graphs.grid((d, d, d))
    ###2D
    elseif dim == 2
        G = Graphs.grid((row, col, 1))
    ###dim > 3
    else
        println("WARNING -> dim not corret, it set 3")
        d = convert(Int, trunc((row * col) ^ (1/3))) + 1
        G = Graphs.grid((d, d, d))
    end

    G_meta = initialize_graph(G, row, col, n_cell, seed)
    return G_meta
end

####Load graphs
function spatial_graph(path::String, seed::MersenneTwister; n_cell::Int = 1)
    mat = readdlm(path, Int)
    n = size(mat, 1) # Numero di nodi
    graph = SimpleGraph(n)
    for i in 1:n, j in 1:i
        if mat[i,j] != 0
            add_edge!(graph, i, j)
        end
    end
    abstract_g = Graph(graph)
    G_meta = initialize_graph(abstract_g, n_cell, seed)
    return G_meta
end

####Load Graphs from DataFrames
function spatial_graph(path_dataframe_edges::String, path_dataframe_labels::String)
    df_edges_file = CSV.File(path_dataframe_edges)
    df_edges = DataFrame(df_edges_file)

    df_labels_file = CSV.File(path_dataframe_labels)
    df_labels = DataFrame(df_labels_file)
    G_meta = DataFrameGraphBridge.metagraph_from_dataframe_JHistint(MetaGraph, df_edges, :origin, :destination, :weight, :weight, df_labels, :label)
    return G_meta
end

"""
Initialize graph with un single cell in the middle of the
lattice/square
"""
###Lattice not regular, initialized graph at node with the most closeness
function initialize_graph(G::AbstractGraph, n_cell::Int,  seed::MersenneTwister)
    G_meta = MetaGraphs.MetaGraph(G)

    if n_cell == 1
        closeness = closeness_centrality(G)
        max_value=maximum(closeness)
        pos_value = findall(p -> p == max_value , closeness)
        v₀ = rand(seed, pos_value)
        id = uuid1(seed)
        for i in 1:nv(G_meta.graph)
            ## Isn't this an 'if'?
            i == v₀ ?
                set_props!(G_meta,
                           i,
                           Dict(:mutation => (1),
                                :name => "node$i",
                                :id => id)) :
                set_props!(G_meta, i, Dict(:name => "node$i"))
        end

    elseif n_cell == 0
        for i in 1:nv(G_meta.graph)
             set_props!(G_meta, i, Dict(:name=> "node$i"))
         end

    else
        positions = Set(1:nv(G))

        for i in 1:n_cell
            id = uuid1(seed)
            pos = rand(seed, positions)
            delete!(positions, pos)
            set_props!(G_meta, pos, Dict(:mutation => (1), :id => id))
        end

        for i in 1:nv(G_meta.graph)
             set_props!(G_meta, i, Dict(:name=> "node$i"))
        end
    end
    return G_meta
end

#SquareGrid
function initialize_graph(G::AbstractGraph,
                          row::Int,
                          col::Int,
                          n_cell::Int,
                          seed::MersenneTwister)

    G_meta = MetaGraphs.MetaGraph(G)

    if n_cell == 1
        nv(G_meta) % 2 == 0 ?
            v₀ = (col * floor(row / 2)) + (floor(col / 2) + 1) :
            v₀ = floor(nv(G) / 2) + 1
        id = uuid1(seed)
        for i in 1:nv(G_meta.graph)
            i == v₀ ?
                set_props!(G_meta,
                           i,
                           Dict(:mutation => (1),
                                :name => "node$i",
                                :id => id)) :
                set_props!(G_meta, i, Dict(:name => "node$i"))
        end

    elseif n_cell == 0
        for i in 1:nv(G_meta.graph)
             set_props!(G_meta, i, Dict(:name=> "node$i"))
        end

    else
        positions = Set(1:nv(G))
        for i in 1:n_cell
            id = uuid1(seed)
            pos = rand(seed, positions)
            delete!(positions, pos)
            set_props!(G_meta, pos, Dict(:mutation => (1), :id => id))
        end
        for i in 1:nv(G_meta.graph)
             set_props!(G_meta, i, Dict(:name=> "node$i"))
        end
    end
    return G_meta
end


"""
Plots a graph in 2D.
"""
function plot_lattice(G::MetaGraph, Set_mut::Vector{Any}; dim::Int=2)
    driver_mut, labels, colors = get_drivermut_name_colors(G, Set_mut)
    mylayout = NetworkLayout.SquareGrid(cols=:auto)
    f, ax, p = graphplot(G,
                         layout = mylayout,
                         node_size = repeat([5], nv(G)),
                         node_color = colors)
    hidedecorations!(ax)
    hidespines!(ax)
    return f, ax, p, colors
end

function plot_lattice_JHistint(G::MetaGraph, Set_mut::Vector{Any}; dim::Int=3)
    driver_mut, labels, colors = get_drivermut_name_colors(G, Set_mut)
    mylayout = NetworkLayout.Spectral(dim=3)
    f, ax, p = graphplot(G,
                         layout = mylayout,
                         node_size = repeat([5], nv(G)),
                         edge_width=1.0,
                         node_color = colors)
    return f, ax, p, colors
end

function plot_lattice_metagraph(G::MetaGraph; dim::Int=3)
    mylayout = NetworkLayout.Spectral(dim=3)
    f, ax, p = graphplot(G,
                         layout = mylayout,
                         node_size = repeat([5], nv(G)),
                         edge_width=1.0)
    return f, ax, p, colors
end


"""
Creates a palette based on the number of driver mutations present.
"""
function color_index(driver_mut::Vector{Any}, Set_mut::Vector{Any})
    colors = []
    n_color = 20
    #if there are more 20 mutation in simulation
    if length(Set_mut) > 20
        n_color = length(Set_mut)
    end
    #create palette color
    nodecolor_range =
        distinguishable_colors(n_color,
                               [RGB(1, 1, 1), RGB(0, 0, 0)],
                               dropseed = true)

    for i in 1:length(driver_mut)
        if driver_mut[i] != []
            idx = findall(x -> x == driver_mut[i], Set_mut)[1]
            push!(colors, nodecolor_range[idx])
        else
            #push!(colors, colorant"white")
            push!(colors, RGB(1, 1, 1))
        end
    end

    return colors
end


"""
Returns drivers mutations and names for each cell.
"""
function get_drivermut_name_colors(G::MetaGraph, Set_mut::Vector{Any})
    driver_muts = []
    names = []
    for i in 1:nv(G)
        has_prop(G, i, :mutation) ?
            push!(driver_muts, get_prop(G, i, :mutation)) :
            push!(driver_muts, [])
        push!(names, get_prop(G, i, :name))
    end
    colors = color_index(driver_muts, Set_mut)
    return driver_muts, names, colors
end


"""
Plots tree.
"""
function plot_tree(Tree::AbstractMetaGraph)
    color = [:black for i in 1:nv(Tree)]
    color[1] = :red
    f, ax, p =
    graphplot(Tree,
              layout = Buchheim(),
              node_color = color,
              nlabels = [string(v) for v in vertices(Tree)])
    hidedecorations!(ax)
    hidespines!(ax)
    return f, ax, p
end


"""
Returns only drivers mutations of each cell
"""
function get_drivermut(G::MetaGraph)
    driver_muts = []
    for i in 1:nv(G)
        if has_prop(G, i, :mutation)
            push!(driver_muts, get_prop(G, i, :mutation))
        end
    end
    return driver_muts
end


"""
Returns all drive mutations.
"""
function get_all_mut(G::MetaGraph)
    all_muts = []
    for i in 1:nv(G)
        has_prop(G, i,:mutation) ?
            push!(all_muts, get_prop(G, i, :mutation)) :
            push!(all_muts, [])
    end
    return all_muts
end


"""
Converts a graph into a dataframe, where each row correspond to space
of lattice.

## type of event || time of event || cell involved || Notes ->
## In notes we put
## 1)Parent, Mut if Event is duplicate with new driver
## 2)Parent if Event is duplicate without new driver
## 3)undef if Event is death
"""
function Graph_to_Dataframe(G::AbstractGraph)
    df = DataFrame(Event = String[],
                   Time = Float64[],
                   Cell = UUID[],
                   Notes = Any[])
    return df
end


"""
Returns a list of cells -> node with a cell
"""
function cells_alive(G::AbstractGraph)
    f_v = filter_vertices(G, :id)
    f_vs = [v for v in f_v]
    return f_vs
end


"""
Returns dict with {id => [neighbors]}
Note: id is identifier for the node, not cell
"""
function cells_neighbors(G::AbstractGraph)
    c_neighbors = Dict{Int64,Any}()
    for c in vertices(G)
        c_n = all_neighbors(G, c)
        c_neighbors[c] = c_n
    end
    return c_neighbors
end


"""
Update df and graph when there is a birth.
"""
function cell_birth(G::AbstractGraph,
                    cell::Int,
                    pos::Int,
                    df::DataFrame,
                    μ_dri::AbstractFloat,
                    set_mut::Vector{Any},
                    α::Vector{Float64},
                    driv_average_advantage::AbstractFloat,
                    driv_std_advantage::AbstractFloat,
                    time::AbstractFloat,
                    ca_subpop::Vector{Any},
                    idx::Int,
                    seed::MersenneTwister)

    #delete old cell if use h_voter or voter
    if has_prop(G, pos, :id)
        filter!(e -> e != pos, ca_subpop[get_prop(G, pos, :Subpop)])
        cell_death(G, pos, df, time)
    end

    id = uuid1(seed)
    id2 = uuid1(seed)
    parent = get_prop(G, cell, :id)
    muts = get_prop(G, cell, :mutation)
    subpop = get_prop(G, cell, :Subpop)
    fitness = get_prop(G, cell, :Fit)

    r = rand(seed)
    if r > μ_dri     # NO driver mut

        ## Update Dataframe
        push!(df, ["Duplicate", time, id, [parent]])
        push!(df, ["Duplicate", time, id2, [parent]])

        ## Update  Graph
        set_props!(G, cell, Dict(:mutation => muts,
                                 :id => id,
                                 :Subpop => subpop,
                                 :Fit => fitness))
        set_props!(G, pos, Dict(:mutation => muts,
                                :id => id2,
                                :Subpop => subpop,
                                :Fit => fitness))
        push!(ca_subpop[idx], pos)

    else     # NEW Driver mut
        possible_mut = Set(1:1000) # Possible mutations
        ms = unique(collect(Iterators.flatten(set_mut))) # mut already present
        [delete!(possible_mut, m) for m in ms][1] #clear

        # Create new "cluster" at mut
        new_drive = rand(seed, possible_mut)
        old_muts = collect(Iterators.flatten(muts))
        new_mut = push!(old_muts, new_drive)

        ## compute new α
        new_alpha_driver(G,
                         cell,
                         set_mut,
                         α,
                         driv_average_advantage,
                         driv_std_advantage,
                         seed)

        push!(set_mut, new_mut)

        #Update DataFrame
        push!(df, ["Mutation", time, id2, [parent, new_drive]])
        push!(df, ["Duplicate", time, id, [parent]])

        ## Update Graph and Subpop
        if rand(seed) < 0.5 # Prob mut is in pos or new pos
            set_props!(G, cell, Dict(:mutation => muts,
                                     :id => id,
                                     :Subpop => subpop,
                                     :Fit => fitness))
            set_props!(G, pos, Dict(:mutation => new_mut,
                                    :id => id2,
                                    :Subpop => length(set_mut),
                                    :Fit => α[end]))
            push!(ca_subpop, [pos])
        else
            set_props!(G, pos, Dict(:mutation => muts,
                                    :id => id,
                                    :Subpop => subpop,
                                    :Fit => fitness))
            set_props!(G, cell, Dict(:mutation => new_mut,
                                     :id => id2,
                                     :Subpop => length(set_mut),
                                     :Fit => α[end]))
            push!(ca_subpop[idx], pos)
            filter!(e -> e != cell, ca_subpop[idx])
            push!(ca_subpop, [cell])
        end
    end
    return length(set_mut)
end

###second method
function cell_birth(G::AbstractGraph,
                    cell::Int,
                    pos::Int,
                    df::DataFrame,
                    μ_dri::AbstractFloat,
                    set_mut::Vector{Any},
                    α::Vector{Float64},
                    edge_list::Matrix{String},
                    driv_adv::Matrix{String},
                    time::AbstractFloat,
                    ca_subpop::Vector{Any},
                    idx::Int,
                    seed::MersenneTwister)

    #delete old cell if use h_voter or voter
    if has_prop(G, pos, :id)
        filter!(e -> e != pos, ca_subpop[get_prop(G, pos, :Subpop)])
        cell_death(G, pos, df, time)
    end

    id = uuid1(seed)
    id2 = uuid1(seed)
    parent = get_prop(G, cell, :id)
    muts = get_prop(G, cell, :mutation)
    subpop = get_prop(G, cell, :Subpop)
    fitness = get_prop(G, cell, :Fit)

    drivers_child = ""
    r = rand(seed)
    if r <= μ_dri
        #controllo che posso eseguire una mutazione
        if typeof(muts) == String
            drivers_father = findall(x -> x == muts, edge_list[:,1])
        else
            drivers_father = findall(x -> x == muts[end], edge_list[:,1])
        end
        if drivers_father != []
            for d_f in drivers_father
                drivers_child = edge_list[d_f, :][2]
                new_mut = []
                if typeof(muts) == String
                    push!(new_mut, muts)
                    push!(new_mut, drivers_child)
                else
                    new_mut = copy(muts)
                    push!(new_mut, drivers_child)
                end
                if new_mut ∉ set_mut
                    break
                else
                    drivers_child = ""
                end

            end
        end
    end

    if drivers_child == ""# NO driver mut
        ## Update Dataframe
        push!(df, ["Duplicate", time, id, [parent]])
        push!(df, ["Duplicate", time, id2, [parent]])

        ## Update  Graph
        set_props!(G, cell, Dict(:mutation => muts,
                                 :id => id,
                                 :Subpop => subpop,
                                 :Fit => fitness))
        set_props!(G, pos, Dict(:mutation => muts,
                                :id => id2,
                                :Subpop => subpop,
                                :Fit => fitness))
        push!(ca_subpop[idx], pos)

    else # NEW Driver mut
        new_mut = []
        if typeof(muts) == String
            push!(new_mut, muts)
            push!(new_mut, drivers_child)
        else
            new_mut = copy(muts)
            push!(new_mut, drivers_child)
        end

        ## compute new α
        new_α_id = findall(x -> x == drivers_child, driv_adv[:,1])[1]
        new_α = parse(Float64, driv_adv[new_α_id, 2])
        push!(α, new_α)
        push!(set_mut, new_mut)
        #Update DataFrame
        push!(df, ["Mutation", time, id2, [parent, new_mut]])
        push!(df, ["Duplicate", time, id, [parent]])

        ## Update Graph and Subpop
        if rand(seed) < 0.5 # Prob mut is in pos or new pos
            set_props!(G, cell, Dict(:mutation => muts,
                                     :id => id,
                                     :Subpop => subpop,
                                     :Fit => fitness))
            set_props!(G, pos, Dict(:mutation => new_mut,
                                    :id => id2,
                                    :Subpop => length(set_mut),
                                    :Fit => α[end]))
            push!(ca_subpop, [pos])
        else
            set_props!(G, pos, Dict(:mutation => muts,
                                    :id => id,
                                    :Subpop => subpop,
                                    :Fit => fitness))
            set_props!(G, cell, Dict(:mutation => new_mut,
                                     :id => id2,
                                     :Subpop => length(set_mut),
                                     :Fit => α[end]))
            push!(ca_subpop[idx], pos)
            filter!(e -> e != cell, ca_subpop[idx])
            push!(ca_subpop, [cell])
        end
    end
    return length(set_mut)
end

"""
Returns alive cells in subpopulation.
"""
function cells_alive_subpop(G::AbstractMetaGraph, Set_mut::Vector{Any})
    num_of_pop = []
    for m in Set_mut
        vs = filter_vertices(G, :mutation, m)
        subpop = [v for v in vs]
        push!(num_of_pop, subpop)
    end
    return num_of_pop
end


"""
Compute new rate birth for new driver.
"""
function new_alpha_driver(G::AbstractMetaGraph,
                          cell::Int,
                          Set_mut::Vector{Any},
                          α::Vector{Float64},
                          driv_average_advantage::AbstractFloat,
                          driv_std_advantage::AbstractFloat,
                          seed::MersenneTwister)
    mut = get_prop(G, cell, :mutation)
    idx = findall(x -> x == mut, Set_mut)[1]
    norm = Normal(driv_average_advantage, driv_std_advantage)
    new_α = α[idx] + abs(rand(seed, norm, 1)[1])
    push!(α, new_α)
end


"""
Cell death.
"""
function cell_death(G::AbstractGraph,
                    cell::Int,
                    df::DataFrame,
                    time::AbstractFloat)
    id = get_prop(G, cell, :id)
    push!(df, ["Death", time, id, undef]) # Update Dataframe

    ## Update Graph
    rem_prop!(G, cell, :id)
    rem_prop!(G, cell, :mutation)
    rem_prop!(G, cell, :Subpop)
    rem_prop!(G, cell, :Fit)
end


"""
Cell migration.
"""
function migration_cell(G::AbstractGraph,
                        cell::Int,
                        pos::Int,
                        df::DataFrame,
                        time::AbstractFloat)
    id = get_prop(G, cell, :id)
    mut = get_prop(G, cell, :mutation)
    subpop = get_prop(G, cell, :Subpop)
    fitness = get_prop(G, cell, :Fit)
    # Update Graph
    rem_prop!(G, cell, :id)
    rem_prop!(G, cell, :mutation)
    rem_prop!(G, cell, :Subpop)
    rem_prop!(G, cell, :Fit)

    set_props!(G, pos, Dict(:mutation => mut,
                            :id => id,
                            :Subpop => subpop,
                            :Fit => fitness))
    # Update Dataframe
    push!(df, ["Migration", time, id, undef])
end


"""
Returns the average on axis y.
"""
function average_axis_y(xs, time_tot, n_cell_alive_tot)
    average_nca = []            # Media num cell alive
    for x in xs
        sum = 0
        for i in 1:length(times_tot)
            p = findall(x .>= times_tot[i])
            s = n_cell_alive_tot[i][p][end]
            sum += s
        end
        push!(average_nca, sum / length(times_tot))
    end
    m = findmax(average_nca)
    return m[1]
end


"""
Evolution simulation.
"""
function simulate_evolution(G::AbstractGraph,
                            Tf::AbstractFloat,
                            rate_birth::AbstractFloat,
                            rate_death::AbstractFloat,
                            rate_migration::AbstractFloat,
                            μ_dri::AbstractFloat,
                            driv_average_advantage::AbstractFloat,
                            driv_std_advantage::AbstractFloat,
                            model::String,
                            seed::MersenneTwister;
                            Time_of_sampling = [],
                            t_bottleneck = [],
                            ratio_bottleneck = [])
    ## Prepare values for simulation
    Time_index = 1
    Gs_plot = []                           #Graphs to plot
    CA_Alive_TOT = []
    cs_alive = cells_alive(G)              # Position cells alive = nodes busy
    list_len_node_occ = []                 # List num cells busy for each event
    n_cs_alive = length(cs_alive)          # Num cells alive
    push!(list_len_node_occ, n_cs_alive)   # Update list nodes

    set_mut_pop = unique(get_drivermut(G)) # All mutation
    df = Graph_to_Dataframe(G)             # create df foreach event
    t_curr = 0.0                           # initial time
    # create array alpha
    α = [rate_birth]
    num_pop = length(set_mut_pop)          # Num mutation tot
    # initialize bottleneck
    if t_bottleneck != []
        id_bottleneck = 1
    else
        id_bottleneck = 0
    end

    #initialize metadata
    for alive in cs_alive
        mutation = get_prop(G, alive, :mutation)
        subpop = findall(m -> m == mutation, set_mut_pop)[1]
        set_prop!(G, alive, :Subpop, subpop)
        set_prop!(G, alive, :Fit, α[subpop])
    end

    cs_neighbors = cells_neighbors(G)
    ca_subpop = cells_alive_subpop(G, set_mut_pop) # Num cells ∀ subpop
    α_subpop_f = []
    #simulation
    while (t_curr < Tf && n_cs_alive > 0) ||
          (rate_death == 0 && n_cs_alive == nv(G))

        α_subpop = []           # alpha for each subpop

        for i in 1:length(α)
            push!(α_subpop, α[i] * length(ca_subpop[i]))
        end
        birth = sum(α_subpop)           # Tot prob birth
        death = rate_death * n_cs_alive # Tot prob death
        M = rate_migration * n_cs_alive # Tot prob migration
        λ = birth + death + M
        t_event = rand(seed, Exponential(1 / λ), 1)[1]
        t_old = copy(t_curr)
        t_curr += t_event

        Aₙ = α_subpop ./ λ
        Bₙ = death / λ
        Mₙ = M / λ

        #plot
        if Time_of_sampling != [] &&
           Time_index <= length(Time_of_sampling) &&
           t_curr > Time_of_sampling[Time_index]
            push!(Gs_plot, copy(G))
            push!(CA_Alive_TOT, length.(ca_subpop))
            Time_index += 1
        end

        ## Probability vector
        prob_vet = vcat(Aₙ, Bₙ, Mₙ)
        prob_cum = cumsum(prob_vet)

        ## Choose event
        k = rand(seed)
        target_subpop = collect(k .<= prob_cum)
        min = findfirst(target_subpop)

        if min == num_pop + 2   # Migration Event
            # Choose random cell
            x = rand(seed, (1:n_cs_alive))
            cell = cs_alive[x]
            # Choose random neighbor
            pos = rand(seed, cs_neighbors[cell])

            if pos ∉ cs_alive #phandom event?
                mut = get_prop(G, cell, :mutation)
                idx = findall(x -> x == mut, set_mut_pop)[1]
                migration_cell(G, cell, pos, df, t_curr)
                push!(cs_alive, pos)
                push!(ca_subpop[idx], pos)
                filter!(e -> e != cell, cs_alive)
                filter!(e -> e != cell, ca_subpop[idx])
                push!(list_len_node_occ, n_cs_alive)
                push!(CA_Alive_TOT, length.(ca_subpop))
            end

        elseif min == num_pop + 1    # Death Event
            # Choose random cell
            x = rand(seed, (1:n_cs_alive))
            cell = cs_alive[x]

            mut = get_prop(G, cell, :mutation)
            idx = findall(x -> x == mut, set_mut_pop)[1]
            cell_death(G, cell, df, t_curr)

            #Update list
            filter!(e -> e != cell, cs_alive)
            filter!(e -> e != cell, ca_subpop[idx])

            n_cs_alive -= 1
            push!(list_len_node_occ, n_cs_alive)
            push!(CA_Alive_TOT, length.(ca_subpop))

        else     # Birth event
            x = rand(seed, 1:length(ca_subpop[min]))
            cell = ca_subpop[min][x]
            pos = rand(seed, cs_neighbors[cell])

            #check model for birth
            rule = true
            if model == "contact"
                if pos ∈ cs_alive
                    rule = false
                end

            elseif  model == "voter"
                if pos ∈ cs_alive
                    if get_prop(G, cell, :Subpop) == get_prop(G, pos, :Subpop)
                        rule = false
                    end
                end

            elseif model == "h_voter" || model == "hvoter"
                if pos ∈ cs_alive
                    if get_prop(G, cell, :Fit) <= get_prop(G, pos, :Fit)
                        rule = false
                    end
                end
            end

            if rule
                num_pop = cell_birth(G,
                                     cell,
                                     pos,
                                     df,
                                     μ_dri,
                                     set_mut_pop,
                                     α,
                                     driv_average_advantage,
                                     driv_std_advantage,
                                     t_curr,
                                     ca_subpop,
                                     min,
                                     seed)

                if pos ∉ cs_alive
                    push!(cs_alive, pos)
                    n_cs_alive += 1
                end
                push!(list_len_node_occ, n_cs_alive)
                push!(CA_Alive_TOT, length.(ca_subpop))
            end
        end

        #bottleneck
        if id_bottleneck != 0
            if t_bottleneck[id_bottleneck] > t_old &&
               t_bottleneck[id_bottleneck] <= t_curr

               ratio = ratio_bottleneck[id_bottleneck]
               #compute the number of cells to be take
               n_cells_taked = n_cs_alive * ratio
               cells_taked = sample(seed,
                                    cs_alive,
                                    convert(Int, trunc(n_cells_taked)),
                                    replace = false)
                for ct in cells_taked
                    #foreach cell execute death event
                    mut = get_prop(G, ct, :mutation)
                    idx = findall(x -> x == mut, set_mut_pop)[1]
                    cell_death(G, ct, df, t_bottleneck[id_bottleneck])
                    #Update list
                    filter!(e -> e != ct, cs_alive)
                    filter!(e -> e != ct, ca_subpop[idx])
                    n_cs_alive -= 1
                end
                #aggiorno tutti i parametri della morte
                push!(list_len_node_occ, n_cs_alive)
                push!(CA_Alive_TOT, length.(ca_subpop))
           end
       end

    end
    return df, G, list_len_node_occ, set_mut_pop, Gs_plot, CA_Alive_TOT, α
end


"""
Evolution simulation #2 given mutation driver tree
"""
function simulate_evolution(G::AbstractGraph,
                            Tf::AbstractFloat,
                            rate_death::AbstractFloat,
                            rate_migration::AbstractFloat,
                            μ_dri::AbstractFloat,
                            model::String,
                            edge_list_path::String,
                            driv_adv_path::String,
                            seed::MersenneTwister;
                            Time_of_sampling = [],
                            t_bottleneck = [],
                            ratio_bottleneck = [])

    ## Prepare values for simulation
    Time_index = 1
    Gs_plot = []                           #Graphs to plot
    CA_Alive_TOT = []
    cs_alive = cells_alive(G)              # Position cells alive = nodes busy
    list_len_node_occ = []                 # List num cells busy for each event
    n_cs_alive = length(cs_alive)          # Num cells alive
    push!(list_len_node_occ, n_cs_alive)   # Update list nodes


    edge_list = readdlm(edge_list_path, String) #read edge_list_driver
    driv_adv = readdlm(driv_adv_path, String) #read advantage_driver
    root = setdiff!(edge_list[:,1], edge_list[:,2])[1]
    for alive in cs_alive
        set_prop!(G, alive, :mutation, root) #change mutation
    end

    set_mut_pop = unique(get_drivermut(G)) # All mutation
    df = Graph_to_Dataframe(G)             # create df foreach event
    t_curr = 0.0                           # initial time
    α = [parse(Float64, driv_adv[1,2])]    # create array alpha
    num_pop = length(set_mut_pop)          # Num mutation tot

    # initialize bottleneck
    if t_bottleneck != []
        id_bottleneck = 1
    else
        id_bottleneck = 0
    end

    #initialize metadata
    for alive in cs_alive
        mutation = get_prop(G, alive, :mutation)
        subpop = findall(m -> m == mutation, set_mut_pop)[1]
        set_prop!(G, alive, :Subpop, subpop)
        set_prop!(G, alive, :Fit, α[subpop])
    end

    cs_neighbors = cells_neighbors(G)
    ca_subpop = cells_alive_subpop(G, set_mut_pop) # Num cells ∀ subpop
    α_subpop_f = []

    #simulation
    while (t_curr < Tf && n_cs_alive > 0) ||
          (rate_death == 0 && n_cs_alive == nv(G))

        α_subpop = []           # alpha for each subpop

        for i in 1:length(α)
            push!(α_subpop, α[i] * length(ca_subpop[i]))
        end

        birth = sum(α_subpop)           # Tot prob birth
        death = rate_death * n_cs_alive # Tot prob death
        M = rate_migration * n_cs_alive # Tot prob migration
        λ = birth + death + M
        t_event = rand(seed, Exponential(1 / λ), 1)[1]
        t_old = copy(t_curr)
        t_curr += t_event

        Aₙ = α_subpop ./ λ
        Bₙ = death / λ
        Mₙ = M / λ

        #plot
        if Time_of_sampling != [] &&
           Time_index <= length(Time_of_sampling) &&
           t_curr > Time_of_sampling[Time_index]
            push!(Gs_plot, copy(G))
            push!(CA_Alive_TOT, length.(ca_subpop))
            Time_index += 1
        end

        ## Probability vector
        prob_vet = vcat(Aₙ, Bₙ, Mₙ)
        prob_cum = cumsum(prob_vet)

        ## Choose event
        k = rand(seed)
        target_subpop = collect(k .<= prob_cum)
        min = findfirst(target_subpop)

        if min == num_pop + 2   # Migration Event
            # Choose random cell
            x = rand(seed, (1:n_cs_alive))
            cell = cs_alive[x]
            # Choose random neighbor
            pos = rand(seed, cs_neighbors[cell])

            if pos ∉ cs_alive #phandom event?
                mut = get_prop(G, cell, :mutation)
                idx = findall(x -> x == mut, set_mut_pop)[1]
                migration_cell(G, cell, pos, df, t_curr)
                push!(cs_alive, pos)
                push!(ca_subpop[idx], pos)
                filter!(e -> e != cell, cs_alive)
                filter!(e -> e != cell, ca_subpop[idx])
                push!(list_len_node_occ, n_cs_alive)
                push!(CA_Alive_TOT, length.(ca_subpop))
                color[cell] = 0
                color[pos] = idx
                push!(colors, copy(color)) # Aggiorno il colore
            end

        elseif min == num_pop + 1    # Death Event
            # Choose random cell
            x = rand(seed, (1:n_cs_alive))
            cell = cs_alive[x]

            mut = get_prop(G, cell, :mutation)
            idx = findall(x -> x == mut, set_mut_pop)[1]
            cell_death(G, cell, df, t_curr)

            #Update list
            filter!(e -> e != cell, cs_alive)
            filter!(e -> e != cell, ca_subpop[idx])

            n_cs_alive -= 1
            push!(list_len_node_occ, n_cs_alive)
            push!(CA_Alive_TOT, length.(ca_subpop))

            color[cell] = 0
            push!(colors, copy(color)) # Aggiorno il colore

        else     # Birth event
            x = rand(seed, 1:length(ca_subpop[min]))
            cell = ca_subpop[min][x]
            pos = rand(seed, cs_neighbors[cell])

            #check model for birth
            rule = true
            if model == "contact"
                if pos ∈ cs_alive
                    rule = false
                end

            elseif  model == "voter"
                if pos ∈ cs_alive
                    if get_prop(G, cell, :Subpop) == get_prop(G, pos, :Subpop)
                        rule = false
                    end
                end

            elseif model == "h_voter" || model == "hvoter"
                if pos ∈ cs_alive
                    if get_prop(G, cell, :Fit) <= get_prop(G, pos, :Fit)
                        rule = false
                    end
                end
            end

            if rule
                num_pop = cell_birth(G,
                                     cell,
                                     pos,
                                     df,
                                     μ_dri,
                                     set_mut_pop,
                                     α,
                                     edge_list,
                                     driv_adv,
                                     t_curr,
                                     ca_subpop,
                                     min,
                                     seed)

                if pos ∉ cs_alive
                    push!(cs_alive, pos)
                    n_cs_alive += 1
                end
                push!(list_len_node_occ, n_cs_alive)
                push!(CA_Alive_TOT, length.(ca_subpop))
            end
        end
        #bottleneck
        if id_bottleneck != 0
            if t_bottleneck[id_bottleneck] > t_old &&
               t_bottleneck[id_bottleneck] <= t_curr

               ratio = ratio_bottleneck[id_bottleneck]
               #compute the number of cells to be take
               n_cells_taked = n_cs_alive * ratio
               cells_taked = sample(seed,
                                    cs_alive,
                                    convert(Int, trunc(n_cells_taked)),
                                    replace = false)
                #foreach cell execute death event
                for ct in cells_taked
                    mut = get_prop(G, ct, :mutation)
                    idx = findall(x -> x == mut, set_mut_pop)[1]
                    cell_death(G, ct, df, t_bottleneck[id_bottleneck])
                    #Update list
                    filter!(e -> e != ct, cs_alive)
                    filter!(e -> e != ct, ca_subpop[idx])
                    n_cs_alive -= 1
                end
                #update parameters
                push!(list_len_node_occ, n_cs_alive)
                push!(CA_Alive_TOT, length.(ca_subpop))

           end
           if id_bottleneck < length(t_bottleneck)
               id_bottleneck += 1
           end
       end #end bottleneck
    end
    if length(set_mut_pop) != length(driv_adv[:,1])
        println("WARNING: This simulation did not generate every subpopulation
                                        specified in the driver mutation tree.")
        println("\t in particolar, follow driver: ")
        all_driver = []
        for d in set_mut_pop
            if typeof(d) == String
                push!(all_driver, d)
            else
                push!(all_driver, d[end])
            end
        end
        driver_print = setdiff!(driv_adv[:,1], all_driver)
        for d in driver_print
            println("\t \t \t \t ",d)
        end
    end
    if 0 in CA_Alive_TOT[end]
        leafs = setdiff!(edge_list[:,2], edge_list[:,1])
        id_zero = findall(x -> x == 0, CA_Alive_TOT[end])
        subpop_extinct = [set_mut_pop[id][end] for id in id_zero]
        driver_print = []
        for leaf in leafs
            if leaf in subpop_extinct
                push!(driver_print, leaf)
            end
        end
        if driver_print != []
            println("WARNING: In this simulation did not present every leaf in the
                        driver mutation tree, because they have become extinct")
            println("\t in particolar, follow driver: ")
            for d in driver_print
                println("\t \t \t \t ",d)
            end
        end
    end
    return df, G, list_len_node_occ, set_mut_pop, Gs_plot, CA_Alive_TOT, α
end

include("sampling_phylogentic_relation_and_genotype.jl")
include("Format_tree.jl")
include("SingleCellExperiment.jl")
include("Bulk_Experiment.jl")
include("CallART.jl")
include("DataFrameGraphBridge.jl")
include("Models.jl")
include("Start.jl")
include("DataFrameAPI.jl")
end # module


### end of file -- J_Space.jl
