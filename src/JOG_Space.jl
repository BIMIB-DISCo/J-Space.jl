module JOG_Space

using MetaGraphs #add property to graph
using Graphs #crete graph
using Plots # library as support
using NetworkLayout #layout for lattice
using GraphMakie #for plot graph
using CairoMakie #library for heatmap
using JSServe #Create page for plot
using WGLMakie #library for 2D/3D plot on web page
using GLMakie #library for 2D/3D plot
using DataFrames #library for the main struct
using Makie
using Makie.MakieLayout
using Random
using UUIDs #library for unique id
#using BenchmarkTools #library for check time and allocations of the function
using Distributions #library for calculate normal distributions
using CSV

#file da esportare
export
#Simulation
spatial_graph, simulate_evolution, simulate_MC,
#Sampling
sampling_phylogentic_relation, create_tree,
#experiment
SC_experiment, experiment_bulk,
#non utili, servono per i plot
plot_lattice_3D_web, animation_2D, create_heatmap
####################### Inizio progetto #########################
#'''/!\ATTENZIONE/!\: per i plot interattivi, bisogna rendere disabled il plot
#                    automatico dato che Atom e jupyter non supportano Makie '''

#create graph, wiht paramaters row,col
#='this function create a lattice with parameters row and col(optional dim),
 and it calls function for initialize lattice with p as probability that
 there are a cell have a driver mutation'=#
function spatial_graph(row::Int, col::Int; dim::Int=1, n_cell::Int=1)
    G = Graphs.grid((row, col,dim))
    G_meta = initialize_graph(G, row,col, n_cell)
    return G_meta
end

function spatial_graph(path::String; n_cell::Int=1)

    G = Graphs.grid((row, col,dim))
    G_meta = initialize_graph(G, row,col, n_cell)
    return G_meta
end

#initialize graph with un single cell in the middle of the lattice/square
function initialize_graph(G::AbstractGraph, row::Int, col::Int, n_cell::Int)
    G_meta = MetaGraphs.MetaGraph(G)
    rng = MersenneTwister(1234)
    if n_cell == 1
        nv(G_meta)%2 == 0  ? v₀ = (col * floor(row/2)) + (floor(col/2) + 1) :
                             v₀ = floor(nv(G)/2) + 1
        id = uuid1(rng)
        for i in 1:nv(G_meta.graph)
            i == v₀ ? set_props!(G_meta, i, Dict(:mutation => (1),
                                                 :name=> "node$i", :id => id)) :
                      set_props!(G_meta, i, Dict(:name=> "node$i"))
        end
    elseif n_cell == 0
        for i in 1:nv(G_meta.graph)
             set_props!(G_meta, i, Dict(:name=> "node$i"))
        end
    else
        positions = Set(1:nv(G))
        for i in 1:n_cell
            id = uuid1(rng)
            pos = rand(positions)
            delete!(positions, pos)
            set_props!(G_meta, pos, Dict(:mutation => (1), :id => id))
        end
    end
    return G_meta
end

#function for plot graph into 2D
function plot_lattice(G::MetaGraph, Set_mut::Vector{Any}; dim::Int=2)
    driver_mut, labels, colors = get_drivermut_name_colors(G, Set_mut)
    GLMakie.activate!()
    mylayout = NetworkLayout.SquareGrid(cols=:auto)
    f, ax, p = graphplot(G, layout=mylayout,
                         node_size = [10 for i in 1:nv(G)],
                         node_color = colors)
    return f, ax, p, colors
end

#color_node, it create a pallete based on driver mutation number there are
function color_index(driver_mut::Vector{Any}, Set_mut::Vector{Any})
    colors = []
    nodecolor_range = distinguishable_colors(20, [RGB(1,1,1), RGB(0,0,0)],
                                                                dropseed=true)
    for i in 1:length(driver_mut)
        if driver_mut[i] != []
            idx = findall(x->x == driver_mut[i],Set_mut)[1] #ne esiste solo 1
            push!(colors, nodecolor_range[idx])
        else
            push!(colors, colorant"white")
        end
    end
    return colors
end

#return drivers mutations and names for each cell
function get_drivermut_name_colors(G::MetaGraph, Set_mut::Vector{Any})
    driver_muts = []
    names = []
    for i in 1:nv(G)
        has_prop(G, i,:mutation) ? push!(driver_muts, get_prop(G,i,:mutation)) :
                                   push!(driver_muts, [])
        push!(names, get_prop(G, i, :name))
    end
    colors = color_index(driver_muts, Set_mut)
    return driver_muts, names, colors
end

#fuction for plot graph into 3D
function plot_lattice_3D_web(G::MetaGraph, set_mut::Vector{Any})
    driver_mut, labels, colors = get_drivermut_name_colors(G, set_mut)
    Page(exportable=true, offline=true)#use JSServe
    WGLMakie.activate!()
    set_theme!(resolution=(800, 600))#use WGLMakie
    f, ax, p = graphplot(G, layout=Spring(dim=3),node_color = colors,
                                                                node_size=100)
    display(f)
end

#return only drivers mutations of each cell
function get_drivermut(G::MetaGraph)
    driver_muts = []
    for i in 1:nv(G)
        if has_prop(G, i, :mutation)
            push!(driver_muts, get_prop(G, i, :mutation))
        end
    end
    return driver_muts
end

#function that return all drive mut
function get_all_mut(G::MetaGraph)
    all_muts = []
    for i in 1:nv(G)
        has_prop(G, i,:mutation) ? push!(all_muts, get_prop(G,i,:mutation)) :
                                   push!(all_muts, [])
    end
    return all_muts
end

#convert a graph into dataframe,in which each row correspond to space of lattice
function Graph_to_Dataframe(G::AbstractGraph)
    df = DataFrame(Event = String[], Time = Float64[], Cell = UUID[],
                                                        Notes = Any[])
    return df
    #type of event || time of event || cell involved || Notes ->
    #In notes we put 1)Parent,Mut if Event is duplicate with new driver
                    #2)Parent if Event is duplicate without new driver
                    #3)undef if Event is death
end

#return a list of cells -> node with a cell
#='''Ho creato una funzione a parte perchè non so se mi interessa solo la
posizione o anche ID della cellula'''=#
function cells_alive(G::AbstractGraph)
    f_v = filter_vertices(G, :id)
    f_vs = [v for v in f_v]
    return f_vs
end

#return dict with {id => [neighbors]} #/!\id è del nodo non della cell
function cells_neighbors(G::AbstractGraph)
    c_neighbors = Dict{Int64,Any}()
    for c in vertices(G)
        c_n = all_neighbors(G, c)
        c_neighbors[c] = c_n
    end
    return c_neighbors
end

#Function that update df and graph when there is a birth
function cell_birth(G::AbstractGraph, cell::Int, pos::Int, df::DataFrame,
                      μ_dri::AbstractFloat, set_mut::Vector{Any},
                      α::Vector{Float64}, driv_average_advantage::AbstractFloat,
                      time::AbstractFloat, ca_subpop::Vector{Any}, idx::Int,
                      rng)
    id = uuid1(rng)
    id2 = uuid1(rng)
    parent = get_prop(G, cell, :id) #trovo il vecchio id
    muts = get_prop(G, cell, :mutation) #prendo la mutazione
    r = rand()#controlle se ho una nuova mut driver
    if r > μ_dri #NO driver mut
        #aggiorno il Dataframe
        push!(df, ["Duplicate", time, id, [parent]])
        push!(df,["Duplicate", time, id2, [parent]])
        #aggiorno il Grafo
        set_props!(G, cell, Dict(:mutation => muts, :id => id))
        set_props!(G, pos, Dict(:mutation => muts, :id => id2))
        push!(ca_subpop[idx], pos)#aggiorno lista dei nodi delle subpop
    else #NEW Driver mut
        possible_mut = Set(1:1000)#setto un insieme di possibili mutazioni
        ms = unique(collect(Iterators.flatten(set_mut)))#mutazioni già presenti
        [delete!(possible_mut, m) for m in ms][1]#tolgo le mut già presenti
        new_drive = rand(possible_mut) #seleziono una mutazione a caso
        old_muts = collect(Iterators.flatten(muts))#flat mut for insert new_drive
        new_mut = push!(old_muts, new_drive)#create new "cluster" at mut
        #calcolo il nuovo alpha
        new_alpha_driver(G, cell, set_mut, α, driv_average_advantage)
        if new_mut == [1]
            println("ms: ",ms)
            println("new_drive: ",new_drive)
            println("old_muts: ",old_muts)
            println("new_mut: ",new_mut)
        end
        push!(set_mut, new_mut)
        push!(df, ["Mutation", time, id2, [parent, new_drive]])
        push!(df, ["Duplicate", time, id, [parent]])
        #aggiorno il Grafo e la lista delle sottopopolazioni
        if rand() < 0.5 #prob casuale che la mutazione sia esterna o interna
            set_props!(G, cell, Dict(:mutation => muts, :id => id))
            set_props!(G, pos, Dict(:mutation => new_mut, :id => id2))
            push!(ca_subpop, [pos]) #aggiungo cell alla nuova subpop
        else
            set_props!(G, pos, Dict(:mutation => muts, :id => id))
            set_props!(G, cell, Dict(:mutation => new_mut, :id => id2))
            push!(ca_subpop[idx], pos) #aggiungo pos alla vecchia subpop
            filter!(e -> e != cell, ca_subpop[idx])# del elm of list of subpop
            push!(ca_subpop, [cell])
        end
    end
    return length(set_mut)
end

#function that return cell for subpop
function cells_alive_subpop(G::AbstractMetaGraph, Set_mut::Vector{Any})
    #num_of_pop = Dict{Any,Int64}()
    num_of_pop = []
    for m in Set_mut
        vs = filter_vertices(G, :mutation, m)
        subpop = [v for v in vs]
        push!(num_of_pop, subpop)
    end
    return num_of_pop
end

#computation new rate birth for new driver
function new_alpha_driver(G::AbstractMetaGraph, cell::Int,
                            Set_mut::Vector{Any}, α::Vector{Float64},
                            driv_average_advantage::AbstractFloat)
    mut = get_prop(G, cell, :mutation)
    idx = findall(x -> x == mut, Set_mut)[1] #ne esiste solo 1
    norm = Normal(driv_average_advantage, driv_average_advantage/2)
    new_α = α[idx] + abs(rand(norm, 1)[1])
    push!(α, new_α)
end

#function death cell
function cell_death(G::AbstractGraph, cell::Int, df::DataFrame,
                                                            time::AbstractFloat)
    id = get_prop(G, cell, :id) #trovo l'id della cellula
    push!(df, ["Death", time, id, undef]) #creo nuovo record
    #aggiorno il grafo
    rem_prop!(G, cell, :id)
    rem_prop!(G, cell, :mutation)
end

function migration_cell(G::AbstractGraph, cell::Int, pos::Int, df::DataFrame,
                                                            time::AbstractFloat)
    id = get_prop(G, cell, :id)
    mut = get_prop(G, cell, :mutation)
    rem_prop!(G, cell, :id)
    rem_prop!(G, cell, :mutation)
    set_props!(G, pos, Dict(:mutation => mut, :id => id))
    push!(df, ["Migration", time, id, undef])
end

#return animation for simulate
function animation_2D(row::Int, col::Int, colors::Vector{Any},
                                                              name_file::String)
    c = Observable(1)
    color_observable = @lift(colors[$c])
    p1 = Graphs.grid((row, col))
    mylayout = NetworkLayout.SquareGrid(cols=:auto)
    f,ax,p = graphplot(p1, layout = mylayout, node_color = color_observable,
                                              node_size = [10 for i in 1:nv(G)])
    hidedecorations!(ax)
    hidespines!(ax)
    l = 1:length(colors)
    record(f, name_file, l; framerate = 24) do i
                c[] = i
            end
end

#function that return value of axis y
function average_axis_y(xs, time_tot, n_cell_alive_tot)
    average_nca = [] #media num cell alive
    for x in xs
        sum = 0
        for i in 1:length(times_tot)
            p = findall(x .>= times_tot[i])
            s = n_cell_alive_tot[i][p][end]
            sum += s
        end
        push!(average_nca, sum/length(times_tot))
    end
    m = findmax(average_nca)
    return m[1]
end

#function that return value of "axis z"-> value of area
function heatmap_value(xs, ys, times_tot, n_cell_alive_tot)
    z = []
    y_old = 0
    for x in xs, y in ys #intervallo del tempo e num cell
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

function create_heatmap(Tf::AbstractFloat, times_tot::Vector{Any},
                            n_cell_alive_tot::Vector{Any};
                                name::String="PlotHeatmap.png", bin::Int=50)
    xs = range(0, Tf, length = bin)#tempo
    #calcolo la media ad ogni tempo
    m = average_axis_y(xs, times_tot, n_cell_alive_tot)
    ys = range(1, trunc(Int, m), length = bin)#num medio cellule
    zs = heatmap_value(xs ,ys ,times_tot, n_cell_alive_tot)
    zs = [convert(AbstractFloat, z) for z in zs]
    zs1 = reshape(zs, bin, bin)
    zs2 = permutedims(zs1)
    fig,ax,p =CairoMakie.heatmap(xs, ys, zs2)
    CairoMakie.Colorbar(fig[1, 2], p)
    save(name, fig, px_per_unit = 1)
    return zs2
end


#function that aim to simulate cancer evolution
function simulate_evolution_color(G::AbstractGraph, Tf::Float64,
                                    rate_birth::Float64, rate_death::Float64,
                                       rate_migration::Float64, μ_dri::Float64,
                                            driv_average_advantage::Float64,
                                               seed::Int; n_save_graph::Int64=1)
    rng = MersenneTwister(seed)
    colors = []
    set_mut_pop = unique(get_drivermut(G)) #Insieme delle mutazioni
    gad = get_all_mut(G)
    color = color_index(gad, set_mut_pop)#ritorno i colori
    push!(colors, color)
    df = Graph_to_Dataframe(G)#creo il dataframe
    t_curr = 0.0 #tempo iniziale
    α = [rate_birth] #creo array di alpha
    num_pop = length(set_mut_pop)#num delle mutazioni tot
    cs_alive=cells_alive(G)#pos cellule vive = nodi occupati
    list_len_node_occ = []#lista del numero di cell occupate ad ogni evento
    n_cs_alive= length(cs_alive)#num cellule vive
    push!(list_len_node_occ, n_cs_alive)#update list nodes occ al tempo t
    #ritorna la lista dei vicini di ogni nodo(sia vuoto che pieno)
    cs_neighbors = cells_neighbors(G)
    ca_subpop = cells_alive_subpop(G, set_mut_pop)#list delle cell ∀ subpop
    while t_curr < Tf && n_cs_alive > 0
        α_subpop = [] #alpha for each subpop
        for i in 1:length(α)
            push!(α_subpop, α[i] * length(ca_subpop[i]))
        end
        birth = sum(α_subpop) #tot prob nascita
        death = rate_death * n_cs_alive #tot prob di morte
        M = rate_migration * n_cs_alive #tot prob di migrazione
        λ = birth + death + M
        t_event = rand(Exponential(1/λ), 1)[1]
        t_curr += t_event
        Aₙ = α_subpop ./ λ
        Bₙ = death / λ
        Mₙ = M / λ
        # probability vector
        prob_vet = vcat(Aₙ, Bₙ, Mₙ)
        prob_cum = cumsum(prob_vet)
        k = rand()
        target_subpop = collect(k .<= prob_cum)
        min = findfirst(target_subpop)
        if min == num_pop + 2#evento migrazione
            x = rand((1:n_cs_alive))#scelgo una cellula a caso
            cell = cs_alive[x] #prendo la cellula
            pos = rand(cs_neighbors[cell])#scelgo una posizione a caso
            if pos ∉ cs_alive #controllo che non è un phandom event
                mut = get_prop(G, cell, :mutation) #recupero mutazione
                idx = findall(x -> x == mut, set_mut_pop)[1]
                migration_cell(G, cell, pos, df, t_curr)
                color = color_index(get_all_mut(G), set_mut_pop)
                push!(colors, color)#aggiorno il colore
                push!(cs_alive, pos) #aggiorno la lista dei nodi occupati
                push!(ca_subpop[idx], pos)#update list of subpop nodes
                filter!(e -> e != cell, cs_alive)#tolgo la vecchia cell occupata
                filter!(e -> e != cell, ca_subpop[idx])# update list of subpop
                push!(list_len_node_occ, n_cs_alive)
            end
        elseif min == num_pop + 1 #evento morte
            #println("Morte!!")
            x = rand((1:n_cs_alive))#scelgo una cellula a caso
            cell = cs_alive[x] #prendo la cellula
            mut = get_prop(G, cell, :mutation) #recupero mutazione
            idx = findall(x -> x == mut, set_mut_pop)[1]
            cell_death(G, cell, df, t_curr)#funzione morte cellula
            color = color_index(get_all_mut(G), set_mut_pop)#get colors
            push!(colors, color)#aggiorno il colore
            filter!(e -> e != cell, cs_alive)#tolgo la vecchia cell occupata
            filter!(e -> e != cell, ca_subpop[idx])# update list of subpop
            n_cs_alive -= 1
            push!(list_len_node_occ, n_cs_alive)
        else #println("Nascita!!")
            x = rand((1:n_cs_alive)) #scelgo una cellula a caso
            cell = cs_alive[x] #prendo la cellula
            pos = rand(cs_neighbors[cell])#scelgo una posizione a caso
            if pos ∉ cs_alive #controllo che non è un phandom event
                #funzione che duplica una cellula
                num_pop = cell_birth(G, cell, pos, df, μ_dri, set_mut_pop,
                                            α, driv_average_advantage, t_curr,
                                                            ca_subpop, min, rng)
                #return colors
                color = color_index(get_all_mut(G), set_mut_pop)
                push!(colors, color)#aggiorno il colore
                push!(cs_alive, pos) #aggiorno la lista dei nodi occupati
                n_cs_alive += 1 #number of cell alive on lattice
                #update list nodes occ al tempo t
                push!(list_len_node_occ, n_cs_alive)
            end
        end
    end

    return df, G, colors, list_len_node_occ
end

#simulate montecarlo algorithm
function simulate_MC(n_sim::Int64, rows::Int64, cols::Int64, Tf::Float64,
                        rate_birth::Float64, rate_death::Float64,
                           rate_migration::AbstractFloat, rate_new_mut::Float64,
                              driv_average_advantage::Float64, seed::Int)
    times_tot = []
    n_cell_alive_tot = []
    for i in 1:n_sim
        println(i)
        g_meta = spatial_graph(rows, cols, dim = 1, n_cell=1)
        df, G, n_cell_alive = simulate_evolution(g_meta, Tf,
                                                rate_birth, rate_death,
                                                rate_migration, rate_new_mut,
                                                driv_average_advantage, seed)
        push!(times_tot, df.Time)
        push!(n_cell_alive_tot, n_cell_alive)
    end
    return times_tot, n_cell_alive_tot
end

#function that aim to simulate cancer evolution
function simulate_evolution(G::AbstractGraph, Tf::AbstractFloat,
                        rate_birth::AbstractFloat, rate_death::AbstractFloat,
                            rate_migration::AbstractFloat, μ_dri::AbstractFloat,
                               driv_average_advantage::AbstractFloat, seed::Int;
                                    n_save_graph::Int=1)
    rng = MersenneTwister(seed)
    e_f = 0
    set_mut_pop = unique(get_drivermut(G)) #Insieme delle mutazioni
    df = Graph_to_Dataframe(G)#creo il dataframe
    t_curr = 0.0 #tempo iniziale
    α = [rate_birth] #creo array di alpha
    num_pop = length(set_mut_pop) #num delle mutazioni tot
    cs_alive=cells_alive(G)#pos cellule vive = nodi occupati
    list_len_node_occ = []#lista del numero di cell occupate ad ogni evento
    n_cs_alive= length(cs_alive)#num cellule vive
    push!(list_len_node_occ, n_cs_alive)#update list nodes occ al tempo t
    #ritorna la lista dei vicini di ogni nodo(sia vuoto che pieno)
    cs_neighbors = cells_neighbors(G)
    ca_subpop = cells_alive_subpop(G, set_mut_pop)#num di cell ∀ subpop
    while t_curr < Tf && n_cs_alive > 0
        α_subpop = [] #alpha for each subpop
        for i in 1:length(α)
            push!(α_subpop, α[i] * length(ca_subpop[i]))
        end
        birth = sum(α_subpop) #tot prob nascita
        death = rate_death * n_cs_alive #tot prob di morte
        M = rate_migration * n_cs_alive #tot prob di migrazione
        λ = birth + death + M
        t_event = rand(Exponential(1/λ), 1)[1]
        t_curr += t_event
        Aₙ = α_subpop ./ λ
        Bₙ = death / λ
        Mₙ = M / λ
        # probability vector
        prob_vet = vcat(Aₙ, Bₙ, Mₙ)
        prob_cum = cumsum(prob_vet)
        k = rand()
        #choose event
        target_subpop = collect(k .<= prob_cum)
        min = findfirst(target_subpop)
        if min == num_pop + 2#evento migrazione
            x = rand((1:n_cs_alive))#scelgo una cellula a caso
            cell = cs_alive[x] #prendo la cellula
            pos = rand(cs_neighbors[cell])#scelgo una posizione a caso
            if pos ∉ cs_alive #controllo che non è un phandom event
                mut = get_prop(G, cell, :mutation) #recupero mutazione
                idx = findall(x -> x == mut, set_mut_pop)[1]
                migration_cell(G, cell, pos, df, t_curr)#eveno migrazione
                push!(cs_alive, pos) #aggiorno la lista dei nodi occupati
                push!(ca_subpop[idx], pos)#update list of subpop nodes
                filter!(e -> e != cell, cs_alive)#tolgo la vecchia cell occupata
                filter!(e -> e != cell, ca_subpop[idx])# update list of subpop
                push!(list_len_node_occ, n_cs_alive)
            else
                e_f +=1
            end
        elseif min == num_pop + 1 #evento morte
            x = rand((1:n_cs_alive))#scelgo una cellula a caso
            cell = cs_alive[x] #prendo la cellula
            mut = get_prop(G, cell, :mutation) #recupero mutazione
            idx = findall(x -> x == mut, set_mut_pop)[1]
            cell_death(G, cell, df, t_curr)#funzione morte cellula
            filter!(e -> e != cell, cs_alive)#tolgo la vecchia cell occupata
            filter!(e -> e != cell, ca_subpop[idx])# update list of subpop
            n_cs_alive -= 1
            push!(list_len_node_occ, n_cs_alive)
        else #evento nascita
            x = rand(1:length(ca_subpop[min]))#scelgo cell a caso
            cell = ca_subpop[min][x]#prendo la cell a caso
            pos = rand(cs_neighbors[cell])#scelgo una posizione a caso
            if pos ∉ cs_alive #controllo che non è un phandom event
                #funzione che duplica una cellula
                num_pop = cell_birth(G, cell, pos, df, μ_dri, set_mut_pop,
                                            α, driv_average_advantage, t_curr,
                                            ca_subpop, min, rng)
                push!(cs_alive, pos) #aggiorno la lista dei nodi occupati
                n_cs_alive += 1 #number of cell alive on lattice
                push!(list_len_node_occ, n_cs_alive)#update list nodes occ a t
            end
        end
    end

    return df, G, list_len_node_occ, set_mut_pop
end

include("sampling_phylogentic_relation_and_genotype.jl")
include("Format_tree.jl")
include("SingleCellExperiment.jl")
include("Bulk_Experiment.jl")
include("CallART.jl")
include("DataFrameGraphBridge.jl")

end # module
