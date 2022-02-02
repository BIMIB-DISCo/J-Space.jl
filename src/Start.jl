#include(".\\JOG_Space.jl")
#include("sampling_phylogentic_relation_and_genotype.jl")
#include("DataFrameGraphBridge.jl")
#using JOG_Space
#using Random
#using TOML
################### function START

function Start(paramaters::String, config::String)
      println("START...")
      #read paramaters
      Par_dict = TOML.parsefile(paramaters)
      #read config
      Conf_dict = TOML.parsefile(config)

      #set general config
      seed_n = Conf_dict["Config"][1]["seed"]
      seed = MersenneTwister(seed_n)
      path_save_file = Conf_dict["Config"][1]["path_to_save_files"]
      path_save_plot = Conf_dict["Config"][1]["path_to_save_plot"]

      ###GRAPH
      println("CREATE/LOAD GRAPH....")
      #check if graph already exist
      if Conf_dict["Config"][1]["Generate_graph"] == 0
            #probabilmente qui cambiare formato!!!!!!!!!
            matrix_adjacency =  Conf_dict["Config"][1]["Path_to_Graph"]
            ##leggere in base a che file è
      else
            row = Par_dict["Graph"][1]["row"]
            col = Par_dict["Graph"][1]["col"]
            dim = Par_dict["Graph"][1]["dim"]
            start_cell = Par_dict["Graph"][1]["N_starting_cells"]
      end
      #create graph
      g_meta = spatial_graph(row, col, seed, dim = dim, n_cell=start_cell)

      ###DYNAMIC
      println("RUN DYNAMIC....")
      Dynamic_dict = Par_dict["Dynamic"][1]
      Model = Dynamic_dict["Model"]
      Time = Dynamic_dict["Max_time"]
      rate_birth = Dynamic_dict["rate_birth"]
      rate_death = Dynamic_dict["rate_death"]
      rate_migration = Dynamic_dict["rate_migration"]
      driv_mut_rate = Dynamic_dict["drive_mut_rate"]
      avg_driv_mut_rate = Dynamic_dict["average_driver_mut_rate"]
      std_driv_mut_rate = Dynamic_dict["std_driver_mut_rate"]
      #run simulation
      df, G, n_cell_alive, set_mut = simulate_evolution(g_meta,
                                                        Time,
                                                        rate_birth,
                                                        rate_death,
                                                        rate_migration,
                                                        driv_mut_rate,
                                                        avg_driv_mut_rate,
                                                        std_driv_mut_rate,
                                                        Model,
                                                        seed)
      #save final configuration
      if Conf_dict["OutputGT"][1]["Final_configuration"] == 1
            f, ax, p, colors = plot_lattice(G, set_mut)
            save(path_save_plot*"\\Final_conf.png", f)
      end

      #save driver list
      if Conf_dict["OutputGT"][1]["Driver_list"] == 1
            CSV.write(path_save_file*"/DriverList.csv",
                      Tables.table(set_mut),
                      header=false)
      end

      ##SAMPLING
      println("SAMPLING....")
      driver_tree = Conf_dict["OutputGT"][1]["Driver_Tree"]
      Sampling_dict = Par_dict["Sampling"][1]
      if Sampling_dict["Random_sampling"] == 1
            num_cell = Sampling_dict["num_cell"]
            matrix_R, tree_mut... = sampling_phylogentic_relation(G,
                                                                  "Random",
                                                                  df,
                                                                  num_cell,
                                                                  set_mut,
                                                                  seed,
                                                                  driver_tree)
      else
            pos_center = Sampling_dict["pos_center"]
            rad_sampl = Sampling_dict["radius_sampling"]
            matrix_R, tree_mut... = sampling_phylogentic_relation(
                                                               G,
                                                               "Neighbourhood",
                                                               df,
                                                               pos_center,
                                                               set_mut,
                                                               seed,
                                                               driver_tree,
                                                               dist = rad_sampl)
      end
      if driver_tree == 1
            tree_mut = tree_mut[1]
            plot_tree(tree_mut, path_save_plot, "\\driver_tree")
      end
      if Conf_dict["OutputGT"][1]["Tree_Newick"] == 1
            tree_red, net = create_tree(matrix_R, true)
            path_complete = path_save_file*"\\formatNewick"#inserire un formato?
            writeTopology(net, path_complete)
      else
            tree_red, net = create_tree(matrix_R, false)
      end

      ##MolecularEvolution
      println("MOLECULAR EVOLUTION....")
      MolEvo_dict = Par_dict["MolecularEvolution"][1]
      is_ISA = MolEvo_dict["type_isa"]
      #controllo se ho la ref
      if Conf_dict["Config"][1]["generate_reference"] == 0
            ref = Conf_dict["Config"][1]["path_reference"]
      else
            ref =  MolEvo_dict["length_genome"]
      end
      if is_ISA == 1
            neu_mut_rate = MolEvo_dict["neut_mut_rate"]
            g_seq, fastaX, position_used = Molecular_evolution(tree_red,
                                                               neu_mut_rate,
                                                               seed,
                                                               ref)
      else
            submodel = MolEvo_dict["sub_model"]
            indel_size = MolEvo_dict["indel_size"]
            indel_rate = MolEvo_dict["indel_rate"]
            branch_length = MolEvo_dict["branch_length"]
            params = IdDict(MolEvo_dict["params"])
            println("params: ",params)
            println("type: ",typeof(params))
            g_seq, fastaX, Tree_SC = singlecell_NoISA(tree_red,
                                                      ref,
                                                      submodel,
                                                      params,
                                                      indel_rate,
                                                      branch_length,
                                                      seed)
      end

      #save fasta
      if Conf_dict["FileOutputExperiments"][1]["Single_cell_fasta"] == 1
            save_fasta(g_seq, fastaX, tree_red, path_save_file)
      end

      ##BulkExperiment
      if Conf_dict["FileOutputExperiments"][1]["VAF_GT"] == 1
            println("BULK EXPERIMENT....")
            BulkExp_dict = Par_dict["BulkExperiment"][1]
            if Conf_dict["FileOutputExperiments"][1]["Bulk_noise"] == 1
                  coverage = BulkExp_dict["coverage"]
                  FP = BulkExp_dict["FP"]
                  FN = BulkExp_dict["FN"]
                  df_bulk, df_bulk_noise... = experiment_bulk(g_seq,
                                                           fastaX,
                                                           position_used,
                                                           path_save_file,
                                                           seed;
                                                           Noise = 1,
                                                           coverage = coverage,
                                                           FP = FP,
                                                           FN = FN)
                  if typeof(df_bulk) == String
                        return "Error -> incorrect parameters for bulk noise"
                  else
                        df_bulk_noise = df_bulk_noise[1]
                  end
            else
                  df_bulk = experiment_bulk(g_seq,
                                            fastaX,
                                            position_used,
                                            path_save_file,
                                            seed)
            end
      end

      ##ART
      println("CALL ART....")
      ART_dict = Par_dict["ART"][1]
      if ART_dict["command"] != ``
            call_Art(ART_dict["command"])
      else
            profile = ART_dict["profile"]
            len_read = ART_dict["len_read"]
            tot_num_reads = ART_dict["tot_num_reads"]
            outfile_prefix = ART_dict["outfile_prefix"]
            if ART_dict["ef"] == 1
                  ef = true
            else
                  ef = false
            end
            if ART_dict["no_ALN"] == 1
                  no_ALN = true
            else
                  no_ALN = false
            end
            if ART_dict["paired_end"] == 1
                  paired_end = true
                  mate_pair = false
            elseif ART_dict["mate_pair"] == 1
                  mate_pair = true
                  paired_end = false
            else
                  mate_pair = false
                  paired_end = false
            end
            mean_fragsize = ART_dict["mean_fragsize"]
            std_fragsize = ART_dict["std_fragsize"]
      end
      #controllare il path_ref
      call_ART(profile, path_ref, #cambiare se ho già la ref o no
                        len_read,
                        tot_num_reads,
                        outfile_prefix,
                        paired_end,
                        seed,
                        ef = ef,
                        mate_pair = mate_pair,
                        mean_fragsize = mean_fragsize,
                        std_fragsize = std_fragsize,
                        no_ALN = no_ALN)
end
#=Start("Parameters.toml", "Config.toml")
################### START, single function
seed = MersenneTwister(1234)
g_meta = spatial_graph(21, 21, seed, dim = 1, n_cell=1)#creo il grafico
df, G, n_cell_alive, set_mut = simulate_evolution(g_meta, 150.0, 0.2, 0.01,
                                                          0.01, 0.05, 0.4, 0.1,
                                                                "contact", seed)
#G,T_final,birth,death,rate_migration,new_drive,driv_average_advantage,
#                                             std_driver_mut_rate, model, seed

times_tot, n_cell_alive_tot = simulate_MC(100, 21, 21, 150.0, 0.2, 0.01,
                                                         0.01, 0.005, 0.4, seed)
#n_simulation, rows, cols, time_f, rate_birth, rate_death, rate_migration,
#                                           new_mut,driv_average_advantage, seed

matrix_R, tree_mut = sampling_phylogentic_relation(G, "Random", df, 10, set_mut,
                                                                        seed, 1)
#matrix_R = sampling_phylogentic_relation(G, "Neighbourhood", df, 100, set_mut,
                                                                      #dist = 8)

tree_red, net = create_tree(matrix_R, true)
writeTopology(net) #mostro la stringa

#25*10⁻⁸ too small
g_seq, fastaX, position_used = Molecular_evolution(tree_red, 0.000025, seed,
                                      6000)


df_bulk, df_bulk_noise = experiment_bulk(g_seq, fastaX, position_used, seed;
                                            Noise=true, coverage=10.0,
                                           FP = 0.0000005, FN = 0.0000002)

#un pò lento
g_seq, fastaX, Tree_SC = singlecell_NoISA(tree_red, g_seq, "isa", 0.0005, 300, 10.0, seed)

################### PLOT
plot_lattice_3D_web(g_meta)#stampo il grafico in 3D #100x100x2 => 20'000 nodi

animation_2D(11, 11, colors, "Video\\prova.mp4")

zs2 = create_heatmap(100.0, times_tot, n_cell_alive_tot,
                        name = "Plot\\plot_heatmap_11x11_100_b20.png", bin = 20)


################### Save DATA
df = DataFrame(times = times_tot, cell = n_cell_alive_tot)
CSV.write("df_11x11_100.csv", df, delim = ",")
df2 = DataFrame(CSV.File("Data\\df_11x11_100.csv"))


#################### phylogentic
#Tree, Ref, Selector, rate_Indel, size_indel,branch_length, seed)
#sequence = genomic_evolution(g_seq, "isa", 0.0005, 300, 10.0, Model_Selector_matrix)
for v in vertices(tree_red)
      println(v)
      println(props(tree_red, v))
end




################### HDF5
using HDF5
fid = h5open("prova", "w") #create file -> qui metto path_save_file prima di prova
close(fid) #close file
=#
