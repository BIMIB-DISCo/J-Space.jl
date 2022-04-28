TARGET: We introduce this example to help perform several simulations of the experiments NGS with an range 
        of value for one parameter
			  
Command 'using' will load the module or package and make its exported names available for direct use.
In particular, for this example, we use module : 
	*J_Space: for principal function
	*Random, PhyloNetworks, Graph, MetaGraph: to auxialiry functions
	*CSV, Tables and DataFrames: to save data
	*FASTX: to save file fasta

Initialize the seed of simulation with function MersenneTwister(Integer Number)

Create lattice that you want to use for your experiment by giving the values 
(e.g. rows = 200, column = 200 and dimension)  'g_meta = spatial_graph(200, 200, seed, dim = 3, n_cell=1)'

We simulate the dynamics of the evolution which returns the configuration final of the graph. 

We carried out the Sampling(function 'sampling_phylogentic_relation') after the simulation of dynamics and
we create the phylogenetic tree in the Newickformat (with function 'create_tree').

Optionally, you can save the fileout obtain until now.

We specify the path where read the files Fasta.

We obtain the leafs of the tree with the function in J_Space 'get_leafs'

We call up the function 'experiment_noISA_sign' that simulate a molecular evolution using signature.
This function returns the file Fasta which will be the variables that we will use for the experiment of the NGS

We save the file Fasta obtained previously.

Now you choose the value of the experiment NGS to be changed
(in this example we choose coverage) 
Setup minimum and maximum value into two several variables.(e.g. coverage_min and coverage_max)  
Choose number of repeat between the above values 
Create a range of values with function range by specifying the chosen values(e.g. values_of_covarage)

Thereafter, we need a loop for execute the experiment with several values.
"for single_value in range_values": the key names are 'for' and 'in', after 'in' put your variable range.(e.g. values_of_covarage).
after 'for' the name of single value of range, in which we will refer to into loop.

N.B. attention to indentation

In this case, you must put the path where saved files Fasta are located and the path where you want to save
files output of the experiment.  
When call up the experiment NGS('call_ART'), you must set the parameters that 
should not be change, and into position of parameter that will change put the name of value written after 
'for' (single_value).