In this example we show how it is possible to use single function of J-SPACE to performa a parameter-scan of a dynamical parameter (e.g., birth rate, death rate, etc.).


In this example we use in particular 2 functions:

	1) spatial_graph(200, 200, seed, dim = 3, n_cell=1)
	cosa fa
	descrizioni parametri
	cosa genera

	2) simulate_evolution(g_meta,200.0,birth_rate, 0.01, 0.01, 0.00001,  0.2, 0.1,  "contact",seed)
	cosa fa
	descrizioni parametri
	cosa genera


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
