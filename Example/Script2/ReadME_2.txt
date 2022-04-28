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