#example per lanciare JOG_Space
using JOG_Space
using BenchmarkTools
using Random
Start("Parameters.toml", "Config.toml")
@btime Start("Parameters.toml", "Config.toml")
pwd()
cd("..")
seed = MersenneTwister(1234)

g_meta = spatial_graph(71, 71, seed, dim =2 , n_cell=1)
g_meta, t... = @timed spatial_graph(71, 71, seed, dim =2 , n_cell=1)
t.time
@time df, G, n_cell_alive, set_mut, Gs_conf, CA_subpop, α_subpop = simulate_evolution(
                                                        g_meta,
                                                        200.0,
                                                        0.4,
                                                        0.01,
                                                        0.0,
                                                        0.0001,
                                                        0.2,
                                                        0.1,
                                                        "contact",
                                                        seed)#,

@time simulate_evolution(
                                                        g_meta,
                                                        200.0,
                                                        0.4,
                                                        0.01,
                                                        0.0,
                                                        0.0001,
                                                        0.2,
                                                        0.1,
                                                        "hvoter",
                                                        seed)#
                                                        #Time_of_sampling = [10,20])
props(G, 100)
using Graphs, MetaGraphs
typeof(LGFormat())
#savegraph("prova", G, GraphIO.GML.GMLFormat())
file = ".\\Final_conf_0.0001-5.lgz"
G = loadgraph("Tree_Fil_0.0001-8.lgz","graph" ,MGFormat())
#using GraphIO.GML
#using ParserCombinator
df = CSV.read("Dinamica_0.0001-9.csv", DataFrame)
df2 = df[df.Event .!= "Duplicate", :]
df3 = df2[df2.Event .!= "Death", :]
t = @timed matrix_R, tree_mut  = sampling_phylogentic_relation(G,
                                                                         "Random",
                                                                         df,
                                                                         100,
                                                                         set_mut,
                                                                         seed,
                                                                         1)


tree_red, net = create_tree(matrix_R, true)
props(tree_red, 1)
writeTopology(net)
params = IdDict( "alpha" => 0.5, "beta" => 0.3)
params = IdDict( "alpha" => 0.5)
g_seq, fastaX, Tree_SC, mutation_driver = Molecular_evolution_NoISA(tree_red,
                                                       10000,
                                                       "JC69",
                                                       params,
                                                       0.000000001,
                                                       100,
                                                       seed,
                                                       set_mut,
                                                       0.5,
                                                       1)#da cambiare in 0 o 1



g_seq, fastaX, position_used,t... =@timed Molecular_evolution_ISA(tree_red,
                                                       0.00000001,
                                                       seed,
                                                       "reference.fasta",
                                                       set_mut)

VAF_GT = JOG_Space.create_bulk_groundtruth(g_seq, fastaX, position_used)
f = histogram(V[:, :VAF])
V =  VAF_GT[VAF_GT.VAF .>= 0.05, :]
VAF_GT[:, :VAF]


adj_matrix= rand(seed, (0,1), (100,100))
using LinearAlgebra
Supper = Symmetric(adj_matrix)


sum = 0
for e in edges(tree_red)
    n_m = get_prop(tree_red, e, :N_Mutations)
    sum += n_m
end

call_ART("HS25",
                            "./fileoutput",
                            150,
                            10,
                            "example",
                            false,
                            seed,
                            sam = false,
                            ef = true,
                            mate_pair = false,
                            mean_fragsize = 200,
                            std_fragsize = 10,
                            no_ALN = false)

l
using Random
seed = MersenneTwister(1234)
call_ART("HS25",
                                   "./ART/",
                                   150,
                                   20,
                                   "exp_paired",
                                   true,
                                   seed,
                                   sam = false,
                                   ef = true,
                                   mate_pair = false,
                                   mean_fragsize = 200,
                                   std_fragsize = 10,
                                   no_ALN = false)
