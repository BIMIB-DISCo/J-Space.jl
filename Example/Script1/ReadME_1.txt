TARGET: We introduce this example to help recall more times the function to simulate evolution with several values of one parameter 

Command 'using' will load the module or package and make its exported names available for direct use.
In particular, for this example, we use module : 
	*J_Space: for principal function
	*Random, PhyloNetworks, Graph, MetaGraph: to auxialiry function
	*CSV, Tables and DataFrames: to save data

Initialize the seed of simulation with function MersenneTwister(Integer Number)

Now choose the value of simulation to be changed
(in this example we choose birth_rate value) 
Setup minimum and maximum value into two several variables.(e.g. birth_rate_min and birth_rate_max)  
Choose number of repeat between the above values 

Create a range of values with function range by specifying the chosen values(e.g. values_of_birth_rate)

Initialize all variables of values that you want save as the newick file or the final configuration of graph. 
(initialize them as vector, e.g. Newick_finals = [])

Now, we need a loop for execute all the values into created range.
"for single_value in range_values": the key names are 'for' and 'in', after 'in' put your variable range.(e.g. values_of_birth_rate).
after 'for' the name of single value of range, in which we will refer to into loop.

N.B. attention to indentation

Create lattice that you want to use for your experiment by giving the values 
(e.g. rows = 200, column = 200 and dimension)  'g_meta = spatial_graph(200, 200, seed, dim = 3, n_cell=1)'

When call up the simulation function, you have to set the parameters that should not be change, and into
position of parameter that will change put the name of value written after 'for' (single_value)

After that, you can to decide save output or call up other function and save later.

In this example, we insert the output of single run into above declared variables, with function 'push!()'.
We carried out the Sampling(function 'sampling_phylogentic_relation') after the simulation of dynamics and
we create the phylogenetic tree in the Newickformat (with function 'create_tree').

Finally, we save all outputfiles that we have obtained from called function.