# SPACE-SIM
## INTRODUCTION 
SPACE-SIM is a Julia package to simulate the genomic evolution and the spatial growth of a cell population and the procedure of sequencing the genome of the sampled cells. SPACE-SIM operate under a broad set of different conditions, phenomenological rules and experimental setting.
Firstly, the software simulates the spatial dynamics of the cells as a continuous-time multi-type birth-death stochastic process on a graph employing different rules of interaction and an optimised Gillespie algorithm. 
Finally, after mimicking a spatial sampling of the tumour cells, SPACE-SIM  returns the phylogenetic tree of the sample and simulate molecular evolution of the genome under the a infinite-site models or a set of different substitution models  with the possibility to include structural variants as the indels. Finally, employing ART, SPACE-SIM   generate the  synthetic  single-end, paired-end/mate-pair reads of three  next-generation sequencing platforms.
For any of the theoretical details please see xxx
 
## REQUIRED  SOFTWARE AND PACKAGE



## INSTALLATION OF SPACE-SIM

## RUN SPACE-SIM
### RUN A SINGLE SIMULATION
The paramenters and the cofinguration of the simulation are managment by the user by modifing the files Parameters.tml and Config.tml that are detailed below.
### RUN A SINGLE FUNCTION


### RUN THE EXAMPLES

## OUTPUTS OF SPACE-SIM


## THE CONFIGURATION FILE OF SPACE-SIM
In the file Config.tml the user can manage the configuration of SPACE-SIM. This file is useful to choose the path where save the files, the desired plot and outputfile.

The following are the paramenters of such file

- seed, the seed of the simulations.
- generate_reference,  if 0 the user should inser the reference genome in fasta format in the folder path_reference. If 1 SPACE-SIM generate a random sequence.
- path_reference, the path of the reference given by user.
- path_to_save_files, path of the folder where the output files should be saved
- path_to_save_plot,  path of the folder where the output plots should be saved
- Generate_graph, if 0 the user should inser the graph of the dynamics as an adjacency matrix. If 1 SPACE-SIM generate regular lattice, the paramenters of such lattice are specified in the fiele Paramenters.toml.
- Tree_Newick, if 1  the phylogenentic tree of the cells is saved as newick file.
- Final_configuration, if 1 return the metagraph (in format) of the final configuration of the lattice.
- Driver_list, if 1 returns the list of the driver mutations
- Driver_Tree, if 1 returns the tree of the driver mutations
- Single_cell_fasta, if 1 save the fasta of the GT sequences of the sampled cells 
- Single_cell_noise, if 0  the sequencing experiment is not performed, if 1 it returns the FASTQ  files of the reads and SAM files with noise  , if 2 it returns the FASTQ  files of the reads, the SAM files with noise and without noise.
- Alignment, if 1 it returns the GT aligment file in ALN format.
 - VAF_GT, if 1 it returns the VAF of the sampled cells  (working only if isa is used).
- Bulk_noise , if 1 it returns an approximate bulk experiment not using ART (working only if isa is used).
- Time_of_sampling, insert an array of times. SPACE-SIM perform the plot of the state of the lattice in that times.
- Dynamic_Clonal_genotype, if 1 plots the dynamics of the clonal genotypes.
- Graph_configuration, if 1  returns the plot of the state of the lattice.


## THE PARAMETERS FILE OF SPACE-SIM
In the file Parameters.tml the user will find all the paramenters of the dynamics, molecular evolution and experiment. 



### Paramenters of the generation of the lattice

- row, integer number of the rows of the regular lattice, not used if a graph is imported
- col, integer number of the columns of the regular lattice, not used if a graph is imported
- dim, integer number of the rows of the regular lattice, not used if a graph is imported. If =1 the simulation is 2D
- N_starting_cells, integer  the number of starting cells.

- matrix_adjacency,  path to the adjacency matrix ("path/matrix/adjacency") of the graph that should be imported.

### Paramenters of the clonal spatial dynamics

- Model,  for selecting the different  Models of interaction, possible values ["contact", "voter", "hvoter"]
- Max_time, real number it is the maximum time of the simulation 
- rate_birth, real number.  Birth rate per cell per unit of times
- rate_death,  real number . Death rate per cell per unit of times
- rate_migration, real number.  Migration rate per cell per unit of times
- drive_mut_rate, probability of generate a new driver (i.e., subclones) after a division event.
- average_driver_mut_rate,  average birth rate advantage a driver mutation 
- std_driver_mut_rate, standard deviation of the birth rate advantage of a driver mutation 

### Paramenters of the sampling
- Random_sampling. if 1  Random sampling, if 0 cirular/spherical sampling 
- num_cell. Number of sampled cells

#### If Random_sampling = 0
 -  pos_center. Integer number, it is the center of the sampling
- radius_sampling. Integer number, it is size radius of the sampling

### Paramenters of the molecular evolution
- length_genome.  Length of the ancestral genome. Used  if ancestral genome is not given. If ancestral genome is given equal to 0.
- type_isa. Integer number, if 1 SPACE-SIM use the ISA to simulate the molecular evolution, if 0 it is necessary to specify the substituion model below.
#### if type_isa = 1
- neut_mut_rate. Rate of mutation per site and per unit of time
#### if type_isa = 0
- sub_model. A string, variable that select the subistituion model, possible value ->[ "JC69","F81","K80", "HKY85","TN93","K81"]
- indel_size. Real number maximum size of indel 
- Lavelette_param. Real number the parameter of  the Lavelette distribution for the size of indels
- indel_rate. Rate of indel per site and per unit of time
- params. Rates of the substitution models
 if sub_model= "JC69" params=
 if sub_model= "F81" params=
 if sub_models= "K80" params=
 if sub_models= "HKY85" params=
 if sub_models = "TN93" params =
 if sub_models = "K81" params =

### Paramenters of the bulk experiment (working only if ISA approximation is used)
- coverage. Real number, average coverage 
 - FP. Real number,  false positive rate
- FN. Real number, false negative rate

### Paramenters of the sequencing experiment (ART)
- command = "". A string, if the user want to do custom calls of ART I 
#### Otherwise for Illumina  sequencing system is possible to compile the following parmenters
- profile. The name of Illumina sequencing system of the built-in profile used for simulation, e.g., "HS25"
- len_read. The length of reads to be simulated
- tot_num_reads.  Number of reads/read pairs to be generated per sequence
 - outfile_prefix.  The prefix of output filename
- paired_end. Integer number,  0  indicate a paired-end read simulation or to generate reads from both ends of amplicons, if 1 a paried_end simulation is performed
	                 
			    
#### if paired_end == 1,  are required the following 
- mean_fragsize. The mean size of DNA/RNA fragments for paired-end simulations
- std_fragsize. The standard deviation of DNA/RNA fragment size for paired-end simulations
- mate_pair. If  0  mate-pair read simulation (controlla!!!)
   NOTE: art will automatically switch to a mate-pair simulation if the given mean fragment size >= 2000

For all paramenters of ART  please see: https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm

## EXAMPLES

## POSSIBLE ISSUES
 - la glmak usa gpu problema con macchine virtuali, per 



See the file `COPYING` for license information.

