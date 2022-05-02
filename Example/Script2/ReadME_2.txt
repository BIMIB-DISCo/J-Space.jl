In this example we show how it is possible to use single function of J-SPACE to performa perform different samplings of the same tumor and generate the sequence of the cells of the sampling using a substitution model based on the signatures. This kind of analysis is useful to study the effect of sampling on the data.

In this example we use in particular 4 functions:

	1) spatial_graph(200, 200, seed, dim = 3, n_cell=1)
	cosa fa
	descrizioni parametri
	cosa genera

	2) simulate_evolution(g_meta,200.0,birth_rate, 0.01, 0.01, 0.00001,  0.2, 0.1,  "contact",seed)
	cosa fa
	descrizioni parametri
	cosa genera
		
	3) sampling_phylogentic_relation(G,
                                                       "Random",
                                                        df,
                                                        100,
                                                        set_mut,
                                                        seed,
                                                        1)
	4)  create_tree(matrix_R, true, 200.0)
	
	
	5) experiment_noISA_sign(
                                  g_sign,
                                  10000, #len genome
                                  "96-SBS", #model
                                  0.001, #mut_rate_avg
                                  0.0, #rate indel
                                  100, #size_indel
                                  seed,
                                  set_mut,
                                  0.5,#lavalette_par
                                  ["SBS6","SBS22"], #used_sign
                                  [0.0], #vector_change_points
                                  [0.5 0.5], #vector_activities
                                  0.8) #ratio_background_signature



The script  work as follow:
- Firstly, we use the command 'using' will load the module or package and make its exported names available for direct use.
In particular, for this example, we use load : 
	*J_Space: for principal function
	*Random, PhyloNetworks, Graph, MetaGraph: to auxialiry function
	*CSV, Tables and DataFrames: to save data

- we Initialize the seed of simulation with function MersenneTwister(Integer Number)
- In this case we want to run the same simulation using different value of the birth rate of the wild type. We set the range of the values into two variables.(e.g. birth_rate_min and birth_rate_max)  
- We choose the number of studied  values in a range  (n_rep)

- Create a range of values with function range by specifying the chosen values(e.g. values_of_birth_rate)

- We initialize all variables that we want to save (e.g., newick files, Newick_finals = [])
- Witha for loop we excute a simulation per point of the range (pay attention to indentation)
- Inside the loop, we create the regular graph with  'g_meta = spatial_graph(200, 200, seed, dim = 3, n_cell=1)'
- Finally, we save all outputfiles that we have obtained from called function.


TARGET: We introduce this example to help perform different samplings on the same configuration of the lattice

Command 'using' will load the module or package and make its exported names available for direct use.
In particular, for this example, we use module : 
	*J_Space: for principal function
	*Random, PhyloNetworks, Graph, MetaGraph: to auxialiry functions
	*CSV, Tables and DataFrames: to save data
	*FASTX: to save file fasta

Initialize the seed of simulation with function MersenneTwister(Integer Number)

Choose the number of sampling repetitions(e.g. n_rep = 10)  


Initialize all variables of values that you want save as the newick file
(initialize them as vector, e.g. Newick_finals = [])

Create lattice that you want to use for your experiment by giving the values 
(e.g. rows = 200, column = 200 and dimension)  'g_meta = spatial_graph(200, 200, seed, dim = 3, n_cell=1)'

We simulate the dynamics of the evolution which returns the configuration final of the graph. 
On this parameter we carried out the samplings

Optionally, you can save the fileout obtain until now.
 
Now, we need a loop for execute the samplings at several times.
"for index in 1:n_rep": the key names are 'for' and 'in', after 'in' put your variable range.(e.g. 1:10).
if you choose n_rep = 10, then you repeat ten times the instruction into loop. 
after 'for' the name of index of range, in which we will refer to into loop.

N.B. attention to indentation

In this case, we don't need a range of values because the sampling done will be influenced by the seed

Foreach sampling, we create tree phylogenetic tree in the Newick format and from the leafs tree 
we save file Fasta created after the simulation of molecular evolution using the signatures.
