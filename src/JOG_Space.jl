### -*- Mode: Julia -*-

### JOG_Space.jl

module JOG_Space

using MetaGraphs                # Add property to graph
using Graphs                    # Create graph
using Plots                     # Library as support
using NetworkLayout             # Layout for lattice
using GraphMakie                # Plot graph
using CairoMakie                # Library for heatmap
using JSServe                   # Create page for plot
using WGLMakie                  # Library for 2D/3D plot on web page
using GLMakie                   # Library for 2D/3D plot
using DataFrames                # Library for the main struct
using Makie
using Makie.MakieLayout
using Random
using UUIDs                     # Library for unique id
#using BenchmarkTools # Library for checking time and allocations of
                     # the function
using Distributions  # Library for calculate normal distributions
using CSV

### File da esportare
export
    ## Simulation
    spatial_graph, simulate_evolution, simulate_MC,
    ## Sampling
    sampling_phylogentic_relation, create_tree,
    ## Experiment
    Molecular_evolution, experiment_bulk, singlecell_NoISA,
    ## non utili, servono per i plot
    plot_lattice_3D_web, animation_2D, create_heatmap


### Inizio progetto
### ===============
###'''/!\ATTENZIONE/!\: per i plot interattivi, bisogna rendere
### disabled il plot automatico dato che Atom e jupyter non supportano
### Makie '''

### create graph, wiht parameters row,col

"""
Creates a lattice with parameters row and col(optional dim), and it
calls function for initialize lattice with p as probability that
there are a cell have a driver mutation
"""
function spatial_graph(row::Int, col::Int, seed::MersenneTwister;
                                                dim::Int = 1, n_cell::Int = 1)
    G = Graphs.grid((row, col, dim))
    G_meta = initialize_graph(G, row, col, n_cell, seed)
    return G_meta
end

function spatial_graph(path::String, seed::MersenneTwister; n_cell::Int = 1)
    G = Graphs.grid((row, col, dim))
    G_meta = initialize_graph(G, row,col, n_cell)
    return G_meta
end


"""
Initialize graph with un single cell in the middle of the
lattice/square
"""
function initialize_graph(G::AbstractGraph, row::Int, col::Int, n_cell::Int,
                                                          seed::MersenneTwister)
    G_meta = MetaGraphs.MetaGraph(G)
    if n_cell == 1
        ## Isn't this an 'if'?
        nv(G_meta) % 2 == 0 ?
            v₀ = (col * floor(row / 2)) + (floor(col / 2) + 1) :
            v₀ = floor(nv(G) / 2) + 1
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
    end
    return G_meta
end


"""
Plots a graph in 2D.
"""
function plot_lattice(G::MetaGraph, Set_mut::Vector{Any}; dim::Int=2)
    driver_mut, labels, colors = get_drivermut_name_colors(G, Set_mut)
    GLMakie.activate!()
    mylayout = NetworkLayout.SquareGrid(cols=:auto)
    f, ax, p = graphplot(G,
                         layout = mylayout,
                         node_size = [10 for i in 1:nv(G)],
                         node_color = colors)
    return f, ax, p, colors
end


"""
Creates a palette based on the number of driver mutations present.
"""
function color_index(driver_mut::Vector{Any}, Set_mut::Vector{Any})
    colors = []
    nodecolor_range =
        distinguishable_colors(20,
                               [RGB(1, 1, 1), RGB(0, 0, 0)],
                               dropseed = true)
    for i in 1:length(driver_mut)
        if driver_mut[i] != []
            idx =
                findall(x -> x == driver_mut[i], Set_mut)[1] # Ne esiste solo 1
            push!(colors, nodecolor_range[idx])
        else
            push!(colors, colorant"white")
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
        ## Isn't this an 'if'?
        has_prop(G, i, :mutation) ?
            push!(driver_muts, get_prop(G,i,:mutation)) :
            push!(driver_muts, [])
        push!(names, get_prop(G, i, :name))
    end
    colors = color_index(driver_muts, Set_mut)
    return driver_muts, names, colors
end


"""
Plots graph in 3D.
"""
function plot_lattice_3D_web(G::MetaGraph, set_mut::Vector{Any})
    driver_mut, labels, colors = get_drivermut_name_colors(G, set_mut)
    Page(exportable = true, offline = true) # use JSServe
    WGLMakie.activate!()
    set_theme!(resolution=(800, 600)) # use WGLMakie
    f, ax, p = graphplot(G,
                         layout = Spring(dim = 3),
                         node_color = colors,
                         node_size = 100)
    display(f)
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
        ## Isn't this an 'if'?
        has_prop(G, i,:mutation) ?
            push!(all_muts, get_prop(G, i, :mutation)) :
            push!(all_muts, [])
    end
    ## Even better; are you sure you cannot just map/filter over the
    ## vertices?
    return all_muts
end


"""
Converts a graph into a dataframe, where each row correspond to space
of lattice.
"""
function Graph_to_Dataframe(G::AbstractGraph)
    df = DataFrame(Event = String[],
                   Time = Float64[],
                   Cell = UUID[],
                   Notes = Any[])
    return df
    ## type of event || time of event || cell involved || Notes ->
    ## In notes we put
    ## 1)Parent, Mut if Event is duplicate with new driver
    ## 2)Parent if Event is duplicate without new driver
    ## 3)undef if Event is death
end


"""
Returns a list of cells -> node with a cell
"""
#='''Ho creato una funzione a parte perchè non so se mi interessa solo la
posizione o anche ID della cellula'''=#
function cells_alive(G::AbstractGraph)
    f_v = filter_vertices(G, :id)
    f_vs = [v for v in f_v]
    return f_vs
end


"""
Returns dict with {id => [neighbors]}
""" #/!\id è del nodo non della cell
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
    id = uuid1(seed)
    id2 = uuid1(seed)
    parent = get_prop(G, cell, :id)     # Trovo il vecchio id
    muts = get_prop(G, cell, :mutation) # Prendo la mutazione
    subpop = get_prop(G, cell, :Subpop) # Prendo la sottopopolazione
    fitness = get_prop(G, cell, :Fit)   # Prendo la fitness
    r = rand(seed)                  # Controlle se ho una nuova mut driver
    if r > μ_dri                # NO driver mut
        ## Aggiorno il Dataframe
        push!(df, ["Duplicate", time, id, [parent]])
        push!(df, ["Duplicate", time, id2, [parent]])

        ## Aggiorno il Grafo
        set_props!(G, cell, Dict(:mutation => muts, :id => id,
                                            :Subpop => subpop, :Fit => fitness))
        set_props!(G, pos, Dict(:mutation => muts, :id => id2,
                                            :Subpop => subpop, :Fit => fitness))
        push!(ca_subpop[idx], pos) # Aggiorno lista dei nodi delle subpop
    else                           # NEW Driver mut
        possible_mut = Set(1:1000) # Setto un insieme di possibili mutazioni
        ms = unique(collect(Iterators.flatten(set_mut))) # Mutazioni già presenti
        [delete!(possible_mut, m) for m in ms][1] # Tolgo le mut già presenti
        new_drive = rand(seed, possible_mut) # Seleziono una mutazione a caso
        old_muts = collect(Iterators.flatten(muts)) # Flat mut for insert new_drive
        new_mut = push!(old_muts, new_drive) # Create new "cluster" at mut

        ## Calcolo il nuovo alpha
        new_alpha_driver(G, cell, set_mut, α, driv_average_advantage,
                                                       driv_std_advantage, seed)
        push!(set_mut, new_mut)
        push!(df, ["Mutation", time, id2, [parent, new_drive]])
        push!(df, ["Duplicate", time, id, [parent]])

        ## Aggiorno il Grafo e la lista delle sottopopolazioni
        if rand(seed) < 0.5 # Prob casuale che la mutazione sia esterna o interna
            set_props!(G, cell, Dict(:mutation => muts, :id => id,
                                            :Subpop => subpop, :Fit => fitness))
            set_props!(G, pos, Dict(:mutation => new_mut, :id => id2,
                                    :Subpop => length(set_mut), :Fit => α[end]))
            push!(ca_subpop, [pos]) # Aggiungo cell alla nuova subpop
        else
            set_props!(G, pos, Dict(:mutation => muts, :id => id,
                                            :Subpop => subpop, :Fit => fitness))
            set_props!(G, cell, Dict(:mutation => new_mut, :id => id2,
                                    :Subpop => length(set_mut), :Fit => α[end]))
            push!(ca_subpop[idx], pos) # Aggiungo pos alla vecchia subpop
            filter!(e -> e != cell, ca_subpop[idx]) # delete elm of list of subpop
            push!(ca_subpop, [cell])
        end
    end
    return length(set_mut)
end


"""
Returns alive cells in subpopulation.
"""
function cells_alive_subpop(G::AbstractMetaGraph, Set_mut::Vector{Any})
    # num_of_pop = Dict{Any,Int64}()
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
    idx = findall(x -> x == mut, Set_mut)[1] # Ne esiste solo 1
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
    id = get_prop(G, cell, :id)           # Trovo l'id della cellula
    push!(df, ["Death", time, id, undef]) # Creo nuovo record

    ## Aggiorno il grafo
    rem_prop!(G, cell, :id)
    rem_prop!(G, cell, :mutation)
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
    rem_prop!(G, cell, :id)
    rem_prop!(G, cell, :mutation)
    rem_prop!(G, cell, :Subpop)
    rem_prop!(G, cell, :Fit)
    set_props!(G, pos, Dict(:mutation => mut, :id => id, :Subpop => subpop,
                                                               :Fit => fitness))
    push!(df, ["Migration", time, id, undef])
end


"""
Return animation for simulation.
"""
function animation_2D(row::Int,
                      col::Int,
                      colors::Vector{Any},
                      name_file::String)
    c = Observable(1)
    color_observable = @lift(colors[$c])
    p1 = Graphs.grid((row, col))
    mylayout = NetworkLayout.SquareGrid(cols=:auto)
    f,ax,p = graphplot(p1,
                       layout = mylayout,
                       node_color = color_observable,
                       node_size = [10 for i in 1:nv(G)])
    hidedecorations!(ax)
    hidespines!(ax)
    l = 1:length(colors)
    record(f, name_file, l; framerate = 24) do i
        c[] = i
    end
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
Return value of 'axis z' -> value of area
"""
function heatmap_value(xs, ys, times_tot, n_cell_alive_tot)
    z = []
    y_old = 0
    for x in xs, y in ys        # Intervallo del tempo e num cell
        sum = 0
        for i in 1:length(times_tot)
            p2 = findall(x .>= times_tot[i])
            nca = n_cell_alive_tot[i][p2][end]
            if nca <= y && nca > y_old
                sum += 1
            elseif y == ys[end] && nca > y
                sum +=1
            end
        end
        y != ys[end] ? y_old = y : y_old = 0
        push!(z, sum)
    end
    return z
end


"""
Creates the heatmap.
"""
function create_heatmap(Tf::AbstractFloat,
                        times_tot::Vector{Any},
                        n_cell_alive_tot::Vector{Any};
                        name::String = "PlotHeatmap.png",
                        bin::Int=50)
    xs = range(0, Tf, length = bin) # Tempo

    ## Calcolo la media ad ogni tempo
    m = average_axis_y(xs, times_tot, n_cell_alive_tot)
    ys = range(1, trunc(Int, m), length = bin) # Num medio cellule
    zs = heatmap_value(xs ,ys ,times_tot, n_cell_alive_tot)
    zs = [convert(AbstractFloat, z) for z in zs]
    zs1 = reshape(zs, bin, bin)
    zs2 = permutedims(zs1)
    fig,ax,p =CairoMakie.heatmap(xs, ys, zs2)
    CairoMakie.Colorbar(fig[1, 2], p)
    save(name, fig, px_per_unit = 1)
    return zs2
end

#=
"""
Evolution coloring function.
"""
function simulate_evolution_color(G::AbstractGraph,
                                  Tf::Float64,
                                  rate_birth::Float64,
                                  rate_death::Float64,
                                  rate_migration::Float64,
                                  μ_dri::Float64,
                                  driv_average_advantage::Float64,
                                  seed::MersenneTwister;
                                  n_save_graph::Int64 = 1)
    colors = []
    set_mut_pop = unique(get_drivermut(G)) # Insieme delle mutazioni
    gad = get_all_mut(G)
    color = color_index(gad, set_mut_pop)  # Ritorno i colori
    push!(colors, color)
    df = Graph_to_Dataframe(G)    # Creo il dataframe
    t_curr = 0.0                  # Tempo iniziale
    α = [rate_birth]              # Creo array di alpha
    num_pop = length(set_mut_pop) # Num delle mutazioni tot
    cs_alive=cells_alive(G)       # Pos cellule vive = nodi occupati
    list_len_node_occ = [] # Lista del numero di cell occupate ad ogni evento
    n_cs_alive= length(cs_alive)         # Num cellule vive
    push!(list_len_node_occ, n_cs_alive) # Update list nodes occ al tempo t

    ## Ritorna la lista dei vicini di ogni nodo(sia vuoto che pieno)
    cs_neighbors = cells_neighbors(G)
    ca_subpop = cells_alive_subpop(G, set_mut_pop) # List delle cell ∀ subpop
    while t_curr < Tf && n_cs_alive > 0
        α_subpop = []           # alpha for each subpop
        for i in 1:length(α)
            push!(α_subpop, α[i] * length(ca_subpop[i]))
        end
        birth = sum(α_subpop)           # Tot prob nascita
        death = rate_death * n_cs_alive # Tot prob di morte
        M = rate_migration * n_cs_alive # Tot prob di migrazione
        λ = birth + death + M
        t_event = rand(seed, Exponential(1 / λ), 1)[1]
        t_curr += t_event
        Aₙ = α_subpop ./ λ
        Bₙ = death / λ
        Mₙ = M / λ

        ## Probability vector
        prob_vet = vcat(Aₙ, Bₙ, Mₙ)
        prob_cum = cumsum(prob_vet)
        k = rand(seed)
        target_subpop = collect(k .<= prob_cum)
        min = findfirst(target_subpop)
        if min == num_pop + 2              # Evento migrazione
            x = rand(seed, (1:n_cs_alive))       # Scelgo una cellula a caso
            cell = cs_alive[x]             # Prendo la cellula
            pos = rand(seed, cs_neighbors[cell]) # Scelgo una posizione a caso
            if pos ∉ cs_alive   # Controllo che non è un phandom event
                mut = get_prop(G, cell, :mutation) # Recupero mutazione
                idx = findall(x -> x == mut, set_mut_pop)[1]
                migration_cell(G, cell, pos, df, t_curr)
                color = color_index(get_all_mut(G), set_mut_pop)
                push!(colors, color) # Aggiorno il colore
                push!(cs_alive, pos) # Aggiorno la lista dei nodi occupati
                push!(ca_subpop[idx], pos) # Update list of subpop nodes
                filter!(e -> e != cell, cs_alive) # Tolgo la vecchia
                                                  # cell occupata
                filter!(e -> e != cell, ca_subpop[idx]) # Update list of subpop
                push!(list_len_node_occ, n_cs_alive)
            end
        elseif min == num_pop + 1    # Evento morte
            ## println("Morte!!")
            x = rand(seed, (1:n_cs_alive)) # Scelgo una cellula a caso
            cell = cs_alive[x]       # Prendo la cellula
            mut = get_prop(G, cell, :mutation) # Recupero mutazione
            idx = findall(x -> x == mut, set_mut_pop)[1]
            cell_death(G, cell, df, t_curr) # Funzione morte cellula
            color = color_index(get_all_mut(G), set_mut_pop) # Get colors
            push!(colors, color)              # Aggiorno il colore
            filter!(e -> e != cell, cs_alive) # Tolgo la vecchia cell occupata
            filter!(e -> e != cell, ca_subpop[idx]) # Update list of subpop
            n_cs_alive -= 1
            push!(list_len_node_occ, n_cs_alive)
        else
            ## println("Nascita!!")
            x = rand(seed, (1:n_cs_alive))       # Scelgo una cellula a caso
            cell = cs_alive[x]             # Prendo la cellula
            pos = rand(seed, cs_neighbors[cell]) # Scelgo una posizione a caso
            if pos ∉ cs_alive   # Controllo che non sia un phantom event

                ## funzione che duplica una cellula
                num_pop = cell_birth(G,
                                     cell,
                                     pos,
                                     df,
                                     μ_dri,
                                     set_mut_pop,
                                     α,
                                     driv_average_advantage,
                                     t_curr,
                                     ca_subpop,
                                     min,
                                     seed)
                ## Return colors
                color = color_index(get_all_mut(G), set_mut_pop)
                push!(colors, color) # Aggiorno il colore
                push!(cs_alive, pos) # Aggiorno la lista dei nodi occupati
                n_cs_alive += 1      # Number of cell alive on lattice

                ## Update list nodes occ al tempo t
                push!(list_len_node_occ, n_cs_alive)
            end
        end
    end
    return df, G, colors, list_len_node_occ
end


"""
Montecarlo simulation algorithm.
"""
function simulate_MC(n_sim::Int64,
                     rows::Int64,
                     cols::Int64,
                     Tf::Float64,
                     rate_birth::Float64,
                     rate_death::Float64,
                     rate_migration::AbstractFloat,
                     rate_new_mut::Float64,
                     driv_average_advantage::Float64,
                     seed::MersenneTwister)
    times_tot = []
    n_cell_alive_tot = []
    for i in 1:n_sim
        ## Print?  Is this debugging?
        println(i)
        g_meta = spatial_graph(rows, cols, dim = 1, n_cell = 1)
        df, G, n_cell_alive =
            simulate_evolution(g_meta,
                               Tf,
                               rate_birth,
                               rate_death,
                               rate_migration,
                               rate_new_mut,
                               driv_average_advantage,
                               seed)
        push!(times_tot, df.Time)
        push!(n_cell_alive_tot, n_cell_alive)
    end
    return times_tot, n_cell_alive_tot
end
=#

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
                            n_save_graph::Int = 1)
    set_mut_pop = unique(get_drivermut(G)) # Insieme delle mutazioni
    df = Graph_to_Dataframe(G)             # Creo il dataframe
    t_curr = 0.0                           # Tempo iniziale
    α = [rate_birth]                       # Creo array di alpha
    num_pop = length(set_mut_pop)          # Num delle mutazioni tot
    cs_alive = cells_alive(G)   # Pos cellule vive = nodi occupati
    list_len_node_occ = [] # Lista del numero di cell occupate ad ogni evento
    n_cs_alive = length(cs_alive)        # Num cellule vive
    push!(list_len_node_occ, n_cs_alive) # Update list nodes occ al tempo t
    #initialize metadata
    for alive in cs_alive
        #println("alive: ", alive)
        mutation = get_prop(G, alive, :mutation)
        #println("mut -> ", mutation)
        subpop = findall(m -> m == mutation, set_mut_pop)[1]
        #println("subpop: ", subpop)
        #println("fit: ", α[subpop])
        set_prop!(G, alive, :Subpop, subpop)
        set_prop!(G, alive, :Fit, α[subpop])
    end
    ## Ritorna la lista dei vicini di ogni nodo(sia vuoto che pieno)
    cs_neighbors = cells_neighbors(G)
    ca_subpop = cells_alive_subpop(G, set_mut_pop) # Num di cell ∀ subpop
    while t_curr < Tf && n_cs_alive > 0
        α_subpop = []           # alpha for each subpop
        for i in 1:length(α)
            push!(α_subpop, α[i] * length(ca_subpop[i]))
        end
        birth = sum(α_subpop)           # Tot prob nascita
        death = rate_death * n_cs_alive # Tot prob di morte
        M = rate_migration * n_cs_alive # Tot prob di migrazione
        λ = birth + death + M
        t_event = rand(seed, Exponential(1 / λ), 1)[1]
        t_curr += t_event
        Aₙ = α_subpop ./ λ
        Bₙ = death / λ
        Mₙ = M / λ

        ## Probability vector
        prob_vet = vcat(Aₙ, Bₙ, Mₙ)
        prob_cum = cumsum(prob_vet)
        k = rand(seed)

        ## Choose event
        target_subpop = collect(k .<= prob_cum)
        min = findfirst(target_subpop)

        if min == num_pop + 2              # Evento migrazione
            x = rand(seed, (1:n_cs_alive))       # Scelgo una cellula a caso
            cell = cs_alive[x]             # Prendo la cellula
            pos = rand(seed, cs_neighbors[cell]) # Scelgo una posizione a caso
            if pos ∉ cs_alive   # Controllo che non sia un phantom event
                mut = get_prop(G, cell, :mutation) # Recupero mutazione
                idx = findall(x -> x == mut, set_mut_pop)[1]
                migration_cell(G, cell, pos, df, t_curr) # Evento migrazione
                push!(cs_alive, pos) # Aggiorno la lista dei nodi occupati
                push!(ca_subpop[idx], pos) # Update list of subpop nodes
                filter!(e -> e != cell, cs_alive) # Tolgo la vecchia cell occupata
                filter!(e -> e != cell, ca_subpop[idx]) # Update list of subpop
                push!(list_len_node_occ, n_cs_alive)
            end
        elseif min == num_pop + 1    # Evento morte
            x = rand(seed, (1:n_cs_alive)) # Scelgo una cellula a caso
            cell = cs_alive[x]       # Prendo la cellula
            mut = get_prop(G, cell, :mutation) # Recupero mutazione
            idx = findall(x -> x == mut, set_mut_pop)[1]
            cell_death(G, cell, df, t_curr) # Funzione morte cellula
            filter!(e -> e != cell, cs_alive) # Tolgo la vecchia cell occupata
            filter!(e -> e != cell, ca_subpop[idx]) # Update list of subpop
            n_cs_alive -= 1
            push!(list_len_node_occ, n_cs_alive)
        else                                   # Evento nascita
            x = rand(seed, 1:length(ca_subpop[min])) # Scelgo cell a caso
            cell = ca_subpop[min][x]           # Prendo la cell a caso
            pos = rand(seed, cs_neighbors[cell]) # Scelgo una posizione a caso
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
            elseif model == "h_voter"
                if pos ∈ cs_alive
                    if get_prop(G, cell, :Fit) < get_prop(G, pos, :Fit)
                        rule = false
                    end
                end
            end
            if rule   # Controllo che non è un phandom event

                ## Funzione che duplica una cellula
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
                push!(cs_alive, pos) # Aggiorno la lista dei nodi occupati
                n_cs_alive += 1      # Number of cell alive on lattice
                push!(list_len_node_occ, n_cs_alive) # Update list nodes occ a t
            end
        end
    end
    return df, G, list_len_node_occ, set_mut_pop
end

### Include?  Non è una cosa che si fa a livello di progetto?

include("sampling_phylogentic_relation_and_genotype.jl")
include("Format_tree.jl")
include("SingleCellExperiment.jl")
include("Bulk_Experiment.jl")
include("CallART.jl")
include("DataFrameGraphBridge.jl")

end # module


### end of file -- JOG_Space.jl
