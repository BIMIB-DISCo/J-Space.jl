### Inizio progetto
### ===============
###'''/!\ATTENZIONE/!\: per i plot interattivi, bisogna rendere
### disabled il plot automatico dato che Atom e jupyter non supportano
### Makie '''

"""
Plots a graph in 2D.
"""
function plot_lattice(G::MetaGraph, Set_mut::Vector{Any}; dim::Int=2)
    driver_mut, labels, colors = get_drivermut_name_colors(G, Set_mut)
    mylayout = NetworkLayout.SquareGrid(cols=:auto)
    f, ax, p = graphplot(G,
                         layout = mylayout,
                         node_size = repeat([5], nv(G)),#[7 for i in 1:nv(G)],
                         #edge_width = repeat([1], ne(G)),
                         #edge_color = repeat([:white], ne(G)),
                         node_color = colors)
    hidedecorations!(ax)
    hidespines!(ax)
    return f, ax, p, colors
end


"""
Plots graph in 3D.
"""
function plot_lattice_3D_web(G::MetaGraph, set_mut::Vector{Any})
    driver_mut, labels, colors = get_drivermut_name_colors(G, set_mut)
    Page(exportable = true, offline = true) # use JSServe
    #WGLMakie.activate!()
    set_theme!(resolution=(800, 600)) # use WGLMakie
    f, ax, p = graphplot(G,
                         layout = Spring(dim = 3),
                         node_color = colors,
                         node_size = 100)
    display(f)
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

function create_bulk_groundtruth(G_seq::LongDNASeq,
                                 Fasta_sample::Vector{Any},
                                 Position::Vector{Any})
    # Dataframe -> POS, REF_nucleotide, %A, %C, %G, %T
    #=
    df = DataFrame(Position = Int[], Reference = DNA[],  A = Any[], C = Any[],
                   G = Any[], T = Any[])
    positions = []
    n_sample = length(Fasta_sample)
    for i = 1:len_ROI
        nucl_ref = G_seq[i]
        nucleotides = [f[i] for f in Fasta_sample]
        nucl_counter = counter(nucleotides)
        if nucl_counter[nucl_ref] < n_sample
            pos = i
            k_dict = keys(nucl_counter.map)
            values = []
            for n in names(df)[2:end]
                if n != string(nucl_ref)
                    n_dna = DNA(collect(n)[1])
                    push!(values, (nucl_counter[n_dna]*100)/n_sample)
                else
                    push!(values, "ref")
                end
            end
            push!(df, [pos, nucl_ref, values[1], values[2], values[3], values[4]])
        end
    end
    =#

    df = DataFrame(MUT = Any[], VAF = AbstractFloat[])
    n_sample = length(Fasta_sample)
    Position_sort = sort(Position)
    for pos in Position_sort
        nucl_ref = G_seq[pos]
        nucleotides = [f[pos] for f in Fasta_sample]
        nucl_counter = counter(nucleotides)
        k_dict = keys(nucl_counter.map)
        for k in k_dict
            # println("k :",k,"\t",typeof(k))
            # println("k != nucl_ref -> ",k != nucl_ref)
            if k != nucl_ref
                mut = string(pos) * "_" * string(nucl_ref) * "->" * string(k)
                vaf = nucl_counter[k]/n_sample
                # println("vaf -> ",vaf)
                push!(df, [mut, vaf])
            end
        end
    end
    return df
end
