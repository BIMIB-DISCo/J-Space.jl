
################### function START

function Start_J_Space(paramaters::String, config::String)
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
      #check signature time
      vector_activities = Par_dict["MolecularEvolution"][1]["vector_activities"]
      if vector_activities != []
            vector_activities = reduce(vcat,transpose.(vector_activities))
            for row in eachrow(vector_activities)
                  if sum(row) <= 0.999
                        return println("Error -> ", row,
                                               " : its sum does not equal 1.0")
                  end
            end
      end
      ###GRAPH
      println("CREATE/LOAD GRAPH....")

      #check if graph already exist
      start_cell = Par_dict["Graph"][1]["N_starting_cells"]
      if Conf_dict["Config"][1]["Generate_graph"] == 0
            path_adj_matrix =  Conf_dict["Config"][1]["Path_to_Graph"]
            g_meta = spatial_graph(path_adj_matrix, seed, n_cell = start_cell)
      else
            row = Par_dict["Graph"][1]["row"]
            col = Par_dict["Graph"][1]["col"]
            dim = Par_dict["Graph"][1]["dim"]
            #create graph
            g_meta = spatial_graph(row,
                                   col,
                                   seed,
                                   dim = dim,
                                   n_cell = start_cell)
      end

      ###DYNAMIC
      println("RUN DYNAMIC....")
      #Dyn_Clon_genotype = Conf_dict["FileOutputPlot"][1]["Dynamic_Clonal_genotype"]
      Graph_configuration = Conf_dict["FileOutputPlot"][1]["Graph_configuration"]
      if Graph_configuration == 1
            Time_of_sampling =  Conf_dict["FileOutputPlot"][1]["Time_of_sampling"]
      else
            Time_of_sampling = []
      end

      Dynamic_dict = Par_dict["Dynamic"][1]
      Model = Dynamic_dict["Model"]
      Time = Dynamic_dict["Max_time"]
      rate_death = Dynamic_dict["rate_death"]
      rate_migration = Dynamic_dict["rate_migration"]
      driv_mut_rate = Dynamic_dict["drive_mut_rate"]
      t_bottleneck = Dynamic_dict["t_bottleneck"]
      ratio_bottleneck = Dynamic_dict["ratio_bottleneck"]
      #run simulation
      if Conf_dict["Config"][1]["Tree_Driver_Configure"] != 1
            rate_birth = Dynamic_dict["rate_birth"]
            avg_driv_mut_rate = Dynamic_dict["average_driver_mut_rate"]
            std_driv_mut_rate = Dynamic_dict["std_driver_mut_rate"]
            df, G, n_cell_alive, set_mut, Gs_conf, CA_subpop, α_subpop =
                         simulate_evolution(g_meta,
                                            Time,
                                            rate_birth,
                                            rate_death,
                                            rate_migration,
                                            driv_mut_rate,
                                            avg_driv_mut_rate,
                                            std_driv_mut_rate,
                                            Model,
                                            seed,
                                            Time_of_sampling = Time_of_sampling,
                                            t_bottleneck = t_bottleneck,
                                            ratio_bottleneck =ratio_bottleneck)
      else
            edge_list_path = Conf_dict["Config"][1]["edgelist_treedriver"]
            driv_adv_path = Conf_dict["Config"][1]["driver_birth_rates"]
            df, G, n_cell_alive, set_mut, Gs_conf, CA_subpop, α_subpop =
                         simulate_evolution(g_meta,
                                            Time,
                                            rate_death,
                                            rate_migration,
                                            driv_mut_rate,
                                            Model,
                                            edge_list_path,
                                            driv_adv_path,
                                            seed,
                                            Time_of_sampling = Time_of_sampling,
                                            t_bottleneck = t_bottleneck,
                                            ratio_bottleneck =ratio_bottleneck)
      end
      CSV.write(path_save_file * "Dinamica.csv",
                df,
                delim = ",")

      CSV.write(path_save_file * "n_cell_alive.csv",
                Tables.table(n_cell_alive),
                header=false)

      List_driver = DataFrame(Driver = set_mut, Fitness = α_subpop)
      CSV.write(path_save_file * "DriverList.csv",
                List_driver,
                delim = ",")

      CSV.write(path_save_file * "CA_subpop.csv",
                Tables.table(CA_subpop),
                header=false)

      #save configuration of time at specific time
      if Graph_configuration == 1
            println("save plot...")
            for i in 1:length(Gs_conf)
                  f, ax, p, colors = plot_lattice(Gs_conf[i], set_mut)
                  if Sys.iswindows()
                        save(path_save_plot
                        * "\\Conf_t_"
                        * string(Time_of_sampling[i])
                        *".png",
                        f)
                  elseif Sys.islinux()
                        save(path_save_plot
                        * "/Conf_t_"
                        * string(Time_of_sampling[i])
                        * ".png",
                        f)
                  end
            end
      end

      if Conf_dict["OutputGT"][1]["Final_configuration"] == 1

            f, ax, p, colors = plot_lattice(G, set_mut)
            if Sys.iswindows()
                  save(path_save_plot * "\\Final_conf.png", f)
            elseif Sys.islinux()
                  save(path_save_plot * "/Final_conf.png", f)
            end
      end
      #save driver list
      if Conf_dict["OutputGT"][1]["Driver_list"] == 1
            List_driver = DataFrame(Driver = set_mut, Fitness = α_subpop)
            if Sys.iswindows()
                  CSV.write(path_save_file * "\\DriverList.csv",
                            List_driver,
                            header=false)

            elseif Sys.islinux()
                  CSV.write(path_save_file * "/DriverList.csv",
                            List_driver,
                            header=false)
            end
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
            f, ax, p = plot_tree(tree_mut)
            if Sys.iswindows()
                  save(path_save_plot * "\\driver_tree.png", f)
            elseif Sys.islinux()
                  save(path_save_plot * "/driver_tree.png", f)
            end
      end

      if Conf_dict["OutputGT"][1]["Tree_Newick"] == 1
            tree_red, net = create_tree(matrix_R, true, Time)
            if Sys.iswindows()
                  path_complete = path_save_file * "\\formatNewick"
                  writeTopology(net, path_complete)
            elseif Sys.islinux()
                  path_complete = path_save_file * "/formatNewick"
                  writeTopology(net, path_complete)
            end
      else
            tree_red = create_tree(matrix_R, false)
      end

      ##MolecularEvolution
      println("MOLECULAR EVOLUTION....")
      MolEvo_dict = Par_dict["MolecularEvolution"][1]
      is_ISA = MolEvo_dict["type_isa"]
      #check if ref is present
      if Conf_dict["Config"][1]["generate_reference"] == 0
            ref = Conf_dict["Config"][1]["path_reference"]
            prob_base = []
      else
            ref =  MolEvo_dict["length_genome"]
            prob_base = MolEvo_dict["prob_base"]
      end

      if is_ISA == 1
            neu_mut_rate = MolEvo_dict["neut_mut_rate"]
            g_seq, fastaX, position_used, mutations_tot =
                                    experiment_ISA(tree_red,
                                                   neu_mut_rate,
                                                   seed,
                                                   ref,
                                                   set_mut,
                                                   frequency_dna = prob_base)
      else
            indel_size = MolEvo_dict["indel_size"]
            lavalette = MolEvo_dict["lavalette_par"]
            indel_rate = MolEvo_dict["indel_rate"]
            submodel = MolEvo_dict["sub_model"]
            if submodel in ["SBS-37","SBS-38"]
                  mut_rate_avg = MolEvo_dict["mut_rate_avg"]
                  used_sign = MolEvo_dict["used_sign"]
                  vector_change_points = MolEvo_dict["vector_change_points"]
                  vector_activities = MolEvo_dict["vector_activities"]
                  vector_activities = reduce(vcat,transpose.(vector_activities))
                  ratio_bg_signature = MolEvo_dict["ratio_background_signature"]
                  g_seq, fastaX, Tree_SC, mutations_tot =
                                   experiment_noISA_sign(tree_red,
                                                    ref,
                                                    submodel,
                                                    mut_rate_avg,
                                                    indel_rate,
                                                    indel_size,
                                                    seed,
                                                    set_mut,
                                                    lavalette,
                                                    used_sign,
                                                    vector_change_points,
                                                    vector_activities,
                                                    ratio_bg_signature,
                                                    frequency_dna = prob_base)

            else
                  params = IdDict(MolEvo_dict["params"][1])
                  approx_snv_indel = MolEvo_dict["approx_snv_indel"]
                  g_seq, fastaX, Tree_SC, mutations_tot =
                                   experiment_noISA(tree_red,
                                                    ref,
                                                    submodel,
                                                    params,
                                                    indel_rate,
                                                    indel_size,
                                                    seed,
                                                    set_mut,
                                                    lavalette,
                                                    approx_snv_indel,
                                                    frequency_dna = prob_base)
            end

      end

      if fastaX == []
            return "Fasta are empty, check your input"
      end
      #save mutations_tot
      if Sys.iswindows()
            CSV.write(path_save_file * "\\Mutations_tot.csv",
                      mutations_tot,
                      header=false)

      else
            CSV.write(path_save_file * "/Mutations_tot.csv",
                      mutations_tot,
                      header=false)
      end

      #save fasta
      if Conf_dict["FileOutputExperiments"][1]["Single_cell_fasta"] == 1
            if Sys.iswindows()
                  save_Fasta_W(g_seq, fastaX, tree_red, path_save_file)
            elseif Sys.islinux()
                  save_Fasta_L(g_seq, fastaX, tree_red, path_save_file)
            end

      end

      ##BulkExperiment
      if Conf_dict["FileOutputExperiments"][1]["VAF_GT"] == 1 && is_ISA == 1
            println("BULK EXPERIMENT....")
            BulkExp_dict = Par_dict["BulkExperiment"][1]
            if Conf_dict["FileOutputExperiments"][1]["Bulk_noise"] == 1
                  coverage = BulkExp_dict["coverage"]
                  FP = BulkExp_dict["FP"]
                  FN = BulkExp_dict["FN"]
                  df_bulk, df_bulk_noise... =
                                           experiment_bulk(g_seq,
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
      SC_noise = Conf_dict["FileOutputExperiments"][1]["Single_cell_noise"]
      if SC_noise != 0
            println("CALL ART....")
            if Sys.iswindows()
                  path_fasta = path_save_file*"\\Fasta output\\"
            else
                  path_fasta = path_save_file*"/Fasta output/"
            end
            if Conf_dict["FileOutputExperiments"][1]["Single_cell_fasta"] == 0
                  println("WARNING -> ART must need file Fasta...")
                  println("I am save the files")
                  save_Fasta(g_seq, fastaX, tree_red, path_save_file)
            end

            ART_dict = Par_dict["ART"][1]

            if ART_dict["command"] != ""
                  call_ART(ART_dict["command"], path_fasta)
            else
                  profile = ART_dict["profile"]
                  len_read = ART_dict["len_read"]
                  tot_num_reads = ART_dict["tot_num_reads"]
                  outfile_prefix = ART_dict["outfile_prefix"]
                  if SC_noise == 1
                        ef = false
                        sam = true
                  else
                        ef = true
                        sam = false
                  end
                  if Conf_dict["FileOutputExperiments"][1]["Alignment"] == 0
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
                  call_ART(profile,
                           path_fasta,
                           path_save_file,
                           len_read,
                           tot_num_reads,
                           outfile_prefix,
                           paired_end,
                           seed,
                           sam = sam,
                           ef = ef,
                           mate_pair = mate_pair,
                           mean_fragsize = mean_fragsize,
                           std_fragsize = std_fragsize,
                           no_ALN = no_ALN)
            end
      end
end


# Interface with JHistint , START J-Space

function Start_J_Space(filepath_reference::AbstractString,
                        filepath_matrix::AbstractString,
                        filepath_file::AbstractString,
                        filpeath_plot::AbstractString,
                        slide_id::AbstractString,
                        filepath_dataframe_edges::AbstractString,
                        filepath_dataframe_labels::AbstractString)

      println("J-SPACE: START... ($slide_id)")
      paramaters =  joinpath(@__DIR__, "..", "Parameters.toml")
      config =  joinpath(@__DIR__, "..", "Config.toml")

      #read paramaters
      Par_dict = TOML.parsefile(paramaters)
      #read config
      Conf_dict = TOML.parsefile(config)

      #set general config
      seed_n = Conf_dict["Config"][1]["seed"]
      seed = MersenneTwister(seed_n)
      # path_save_file = Conf_dict["Config"][1]["path_to_save_files"]
      path_save_file = filepath_file
      # path_save_plot = Conf_dict["Config"][1]["path_to_save_plot"]
      path_save_plot = filpeath_plot
      #check signature time
      vector_activities = Par_dict["MolecularEvolution"][1]["vector_activities"]
      if vector_activities != []
            vector_activities = reduce(vcat,transpose.(vector_activities))
            for row in eachrow(vector_activities)
                  if sum(row) <= 0.999
                        return println("Error -> ", row,
                                               " : its sum does not equal 1.0")
                  end
            end
      end
      ###GRAPH
      println("J-SPACE: CREATE/LOAD GRAPH... ($slide_id)")

      #check if graph already exist
      start_cell = Par_dict["Graph"][1]["N_starting_cells"]
      if Conf_dict["Config"][1]["Generate_graph"] == 0
            path_adj_matrix = filepath_matrix
            # path_adj_matrix =  Conf_dict["Config"][1]["Path_to_Graph"]
            g_meta = spatial_graph(path_adj_matrix, seed, n_cell = start_cell)
            # g_meta = spatial_graph(filepath_dataframe_edges, filepath_dataframe_labels)
      else
            row = Par_dict["Graph"][1]["row"]
            col = Par_dict["Graph"][1]["col"]
            dim = Par_dict["Graph"][1]["dim"]
            #create graph
            g_meta = spatial_graph(row,
                                   col,
                                   seed,
                                   dim = dim,
                                   n_cell = start_cell)
      end

      ###DYNAMIC
      println("J-SPACE: RUN DYNAMIC... ($slide_id)")
      #Dyn_Clon_genotype = Conf_dict["FileOutputPlot"][1]["Dynamic_Clonal_genotype"]
      Graph_configuration = Conf_dict["FileOutputPlot"][1]["Graph_configuration"]
      if Graph_configuration == 1
            Time_of_sampling =  Conf_dict["FileOutputPlot"][1]["Time_of_sampling"]
      else
            Time_of_sampling = []
      end

      Dynamic_dict = Par_dict["Dynamic"][1]
      Model = Dynamic_dict["Model"]
      Time = Dynamic_dict["Max_time"]
      rate_death = Dynamic_dict["rate_death"]
      rate_migration = Dynamic_dict["rate_migration"]
      driv_mut_rate = Dynamic_dict["drive_mut_rate"]
      t_bottleneck = Dynamic_dict["t_bottleneck"]
      ratio_bottleneck = Dynamic_dict["ratio_bottleneck"]
      #run simulation
      if Conf_dict["Config"][1]["Tree_Driver_Configure"] != 1
            rate_birth = Dynamic_dict["rate_birth"]
            avg_driv_mut_rate = Dynamic_dict["average_driver_mut_rate"]
            std_driv_mut_rate = Dynamic_dict["std_driver_mut_rate"]

            # SIMULATION BLOCKED ON METAGRAPH, ISSUE ON GITHUB
            df, G, n_cell_alive, set_mut, Gs_conf, CA_subpop, α_subpop =
                         simulate_evolution(g_meta,
                                            Time,
                                            rate_birth,
                                            rate_death,
                                            rate_migration,
                                            driv_mut_rate,
                                            avg_driv_mut_rate,
                                            std_driv_mut_rate,
                                            Model,
                                            seed,
                                            Time_of_sampling = Time_of_sampling,
                                            t_bottleneck = t_bottleneck,
                                            ratio_bottleneck =ratio_bottleneck)
      else
            edge_list_path = Conf_dict["Config"][1]["edgelist_treedriver"]
            driv_adv_path = Conf_dict["Config"][1]["driver_birth_rates"]
            df, G, n_cell_alive, set_mut, Gs_conf, CA_subpop, α_subpop =
                         simulate_evolution(g_meta,
                                            Time,
                                            rate_death,
                                            rate_migration,
                                            driv_mut_rate,
                                            Model,
                                            edge_list_path,
                                            driv_adv_path,
                                            seed,
                                            Time_of_sampling = Time_of_sampling,
                                            t_bottleneck = t_bottleneck,
                                            ratio_bottleneck =ratio_bottleneck)
      end
      CSV.write(path_save_file * "/Dinamica.csv",
                df,
                delim = ",")

      CSV.write(path_save_file * "/n_cell_alive.csv",
                Tables.table(n_cell_alive),
                header=false)

      List_driver = DataFrame(Driver = set_mut, Fitness = α_subpop)
      CSV.write(path_save_file * "/DriverList.csv",
                List_driver,
                delim = ",")

      CSV.write(path_save_file * "/CA_subpop.csv",
                Tables.table(CA_subpop),
                header=false)

      #save configuration of time at specific time
      if Graph_configuration == 1
            println("J-SPACE: SAVE PLOT... ($slide_id)")
            for i in 1:length(Gs_conf)
                  f, ax, p, colors = plot_lattice_JHistint(Gs_conf[i], set_mut)
                  if Sys.iswindows()
                        save(path_save_plot
                        * "\\Conf_t_"
                        * string(Time_of_sampling[i])
                        *".png",
                        f)
                  elseif Sys.islinux()
                        save(path_save_plot
                        * "/Conf_t_"
                        * string(Time_of_sampling[i])
                        * ".png",
                        f)
                  end
            end
      end

      if Conf_dict["OutputGT"][1]["Final_configuration"] == 1

            f, ax, p, colors = plot_lattice_JHistint(G, set_mut)
            if Sys.iswindows()
                  save(path_save_plot * "\\Final_conf.png", f)
            elseif Sys.islinux()
                  save(path_save_plot * "/Final_conf.png", f)
            end
      end
      #save driver list
      if Conf_dict["OutputGT"][1]["Driver_list"] == 1
            List_driver = DataFrame(Driver = set_mut, Fitness = α_subpop)
            if Sys.iswindows()
                  CSV.write(path_save_file * "\\DriverList.csv",
                            List_driver,
                            header=false)

            elseif Sys.islinux()
                  CSV.write(path_save_file * "/DriverList.csv",
                            List_driver,
                            header=false)
            end
      end

      ##SAMPLING
      println("J-SPACE: SAMPLING.... ($slide_id)")
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
            f, ax, p = plot_tree(tree_mut)
            if Sys.iswindows()
                  save(path_save_plot * "\\driver_tree.png", f)
            elseif Sys.islinux()
                  save(path_save_plot * "/driver_tree.png", f)
            end
      end

      if Conf_dict["OutputGT"][1]["Tree_Newick"] == 1
            tree_red, net = create_tree(matrix_R, true, Time)
            if Sys.iswindows()
                  path_complete = path_save_file * "\\formatNewick"
                  writeTopology(net, path_complete)
            elseif Sys.islinux()
                  path_complete = path_save_file * "/formatNewick"
                  writeTopology(net, path_complete)
            end
      else
            tree_red = create_tree(matrix_R, false)
      end

      ##MolecularEvolution
      println("J-SPACE: MOLECULAR EVOLUTION... ($slide_id)")
      MolEvo_dict = Par_dict["MolecularEvolution"][1]
      is_ISA = MolEvo_dict["type_isa"]
      #check if ref is present
      if Conf_dict["Config"][1]["generate_reference"] == 0
            # ref = Conf_dict["Config"][1]["path_reference"]
            ref = filepath_reference
            prob_base = []
      else
            ref =  MolEvo_dict["length_genome"]
            prob_base = MolEvo_dict["prob_base"]
      end

      if is_ISA == 1
            neu_mut_rate = MolEvo_dict["neut_mut_rate"]
            g_seq, fastaX, position_used, mutations_tot =
                                    experiment_ISA(tree_red,
                                                   neu_mut_rate,
                                                   seed,
                                                   ref,
                                                   set_mut,
                                                   frequency_dna = prob_base)
      else
            indel_size = MolEvo_dict["indel_size"]
            lavalette = MolEvo_dict["lavalette_par"]
            indel_rate = MolEvo_dict["indel_rate"]
            submodel = MolEvo_dict["sub_model"]
            if submodel in ["SBS-37","SBS-38"]
                  mut_rate_avg = MolEvo_dict["mut_rate_avg"]
                  used_sign = MolEvo_dict["used_sign"]
                  vector_change_points = MolEvo_dict["vector_change_points"]
                  vector_activities = MolEvo_dict["vector_activities"]
                  vector_activities = reduce(vcat,transpose.(vector_activities))
                  ratio_bg_signature = MolEvo_dict["ratio_background_signature"]
                  g_seq, fastaX, Tree_SC, mutations_tot =
                                   experiment_noISA_sign(tree_red,
                                                    ref,
                                                    submodel,
                                                    mut_rate_avg,
                                                    indel_rate,
                                                    indel_size,
                                                    seed,
                                                    set_mut,
                                                    lavalette,
                                                    used_sign,
                                                    vector_change_points,
                                                    vector_activities,
                                                    ratio_bg_signature,
                                                    frequency_dna = prob_base)

            else
                  params = IdDict(MolEvo_dict["params"][1])
                  approx_snv_indel = MolEvo_dict["approx_snv_indel"]
                  g_seq, fastaX, Tree_SC, mutations_tot =
                                   experiment_noISA(tree_red,
                                                    ref,
                                                    submodel,
                                                    params,
                                                    indel_rate,
                                                    indel_size,
                                                    seed,
                                                    set_mut,
                                                    lavalette,
                                                    approx_snv_indel,
                                                    frequency_dna = prob_base)
            end

      end

      if fastaX == []
            return "Fasta are empty, check your input"
      end
      #save mutations_tot
      if Sys.iswindows()
            CSV.write(path_save_file * "\\Mutations_tot.csv",
                      mutations_tot,
                      header=false)

      else
            CSV.write(path_save_file * "/Mutations_tot.csv",
                      mutations_tot,
                      header=false)
      end

      #save fasta
      if Conf_dict["FileOutputExperiments"][1]["Single_cell_fasta"] == 1
            if Sys.iswindows()
                  save_Fasta_W(g_seq, fastaX, tree_red, path_save_file)
            elseif Sys.islinux()
                  save_Fasta_L(g_seq, fastaX, tree_red, path_save_file)
            end

      end

      ##BulkExperiment
      if Conf_dict["FileOutputExperiments"][1]["VAF_GT"] == 1 && is_ISA == 1
            println("J-SPACE: BULK EXPERIMENT... ($slide_id)")
            BulkExp_dict = Par_dict["BulkExperiment"][1]
            if Conf_dict["FileOutputExperiments"][1]["Bulk_noise"] == 1
                  coverage = BulkExp_dict["coverage"]
                  FP = BulkExp_dict["FP"]
                  FN = BulkExp_dict["FN"]
                  df_bulk, df_bulk_noise... =
                                           experiment_bulk(g_seq,
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

       SC_noise = Conf_dict["FileOutputExperiments"][1]["Single_cell_noise"]
       if SC_noise != 0
           """
        println("J-SPACE: CALL ART... ($slide_id)")
            if Sys.iswindows()
                 path_fasta = path_save_file*"\\Fasta output\\"
            else
                  path_fasta = path_save_file*"/Fasta output/"
            end
            if Conf_dict["FileOutputExperiments"][1]["Single_cell_fasta"] == 0
                  println("WARNING -> ART must need file Fasta...")
                  println("I am save the files")
                  save_Fasta(g_seq, fastaX, tree_red, path_save_file)
            end

            ART_dict = Par_dict["ART"][1]

            if ART_dict["command"] != ""
                  call_ART(ART_dict["command"], path_fasta)
            else
                  profile = ART_dict["profile"]
                  len_read = ART_dict["len_read"]
                  tot_num_reads = ART_dict["tot_num_reads"]
                  outfile_prefix = ART_dict["outfile_prefix"]
                  if SC_noise == 1
                        ef = false
                        sam = true
                  else
                        ef = true
                        sam = false
                  end
                  if Conf_dict["FileOutputExperiments"][1]["Alignment"] == 0
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
                  call_ART(profile,
                           path_fasta,
                           path_save_file,
                           len_read,
                           tot_num_reads,
                           outfile_prefix,
                           paired_end,
                           seed,
                           sam = sam,
                           ef = ef,
                           mate_pair = mate_pair,
                           mean_fragsize = mean_fragsize,
                           std_fragsize = std_fragsize,
                           no_ALN = no_ALN)
            end
            """
      end
      println("J-SPACE: SIMULATION TERMINATED... ($slide_id)")
      println("---------------------------------------------")
end
