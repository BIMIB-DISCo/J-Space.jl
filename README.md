# SPACE-SIM
## INTRODUCTION 
SPACE-SIM is a Julia package to simulate the genomic evolution and the spatial growth of a cell population and the procedure of sequencing the genome of the sampled cells. SPACE-SIM operate under a broad set of different conditions, phenomenological rules and experimental setting.
Firstly, the software simulates the spatial dynamics of the cells as a continuous-time multi-type birth-death stochastic process on a graph employing different rules of interaction and an optimised Gillespie algorithm. 
Finally, after mimicking a spatial sampling of the tumour cells, SPACE-SIM  returns the phylogenetic tree of the sample and simulate molecular evolution of the genome under the a infinite-site models or a set of different substitution models  with the possibility to include structural variants as the indels. Finally, employing ART, SPACE-SIM   generate the  synthetic  single-end, paired-end/mate-pair reads of three  next-generation sequencing platforms.
For any of the theoretical details please see xxx
 
## REQUIRED  SOFTWARE AND PACKAGE



## INSTALLATION OF SPACE-SIM

## RUN SPACE-SIM

The paramenters and the cofinguration of the simulation are managment by the user by modifing the files Parameters.toml and Config.toml that are detailed below.

## OUTPUTS OF SPACE-SIM


## THE CONFIGURATION FILE OF SPACE-SIM
In the file Config.toml the user can manage the configuration of SPACE-SIM. This file is useful to choose the path where save the files, the desired plot and outputfile.

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

#SONO QUI (FABRIZIO)

## THE PARAMETERS FILE OF SPACE-SIM
In the file Parameters.toml the user will find all the paramenters of the dynamics, molecular evolution and experiment. 



### Paramenters of the generation of the lattice

- row = 21 num row of lattice
- col = 21 num col of lattice
- dim = 1 deep of lattice
- N_starting_cells, integer  the number of starting cells.

- matrix_adjacency = "path/matrix/adjacency" #or 1, if already exist

### Paramenters of the clonal spatial dynamics

- Model, choose the different  Model of interaction ["contact", "voter", "hvoter"]
- Max_time = 200.0 Time of the simulation
- rate_birth = 0.2 Rate of birth
- rate_death = 0.01 Rate of Death
- rate_migration = 0.01 Rate of migration
- drive_mut_rate = 0.01 Rate of mew drive mutation
- average_driver_mut_rate = 0.4 Rate of average of the new driver mutation during proliferation
- std_driver_mut_rate = 0.1 Rate of the standard deviation of average_driver_mut_rate

### Paramenters of the sampling
- Random_sampling = 1 hoose mode sampling: 1 -> Random, 0 -> Neighbourhood
- num_cell = 10 Num of cells of the sampling

# If Random_sampling = 0

 -  pos_center = 0 center of neighborshood, if = 0 -> position random
- radius_sampling = 10 size radius of the sampling

### Paramenters of the molecular evolution
- length_genome = 6000 length of genome reference, if ref not given, else = 0
- type_isa = 0 type of model evolution, if value 0 specificate submodel
*** if type_isa = 1 ***
- neut_mut_rate = 0.000025 rate of neutral mutation
*** if type_isa = 0 ***
- sub_model = "K80" decide substitution matrix for mutation, possible value ->[]
- indel_size = 300
- indel_rate = 0.0005
- branch_length = 10.0
- params = [{"alpha" = 0.5, "beta" = 0.3}]

### Paramenters of the bulk experiment (working only if ISA approximation is used)
- coverage = 10.0
 - FP = 0.0000005 rate false positive
- FN = 0.0000002 rate false negative

### Paramenters of the sequencing experiment (ART)
- command = ""  user want write your call
*** Otherwise write param aters ***
- profile = "HS25" he name of Illumina sequencing system of the built-in profile used for simulation
- len_read = 150 the length of reads to be simulated
- tot_num_reads = 10 number of reads/read pairs to be generated per sequence
 - outfile_prefix = "example" the prefix of output filename

- paired_end = 0  indicate a paired-end read simulation or to generate reads from both ends of amplicons
	                    NOTE: art will automatically switch to a mate-pair simulation if the given mean fragment size >= 2000
*** if paired_end == 1, they are require **
- mean_fragsize = 200 the mean size of DNA/RNA fragments for paired-end simulations
- std_fragsize = 10 the standard deviation of DNA/RNA fragment size for paired-end simulations

- mate_pair = 0 indicate a mate-pair read simulation




## EXAMPLES

## POSSIBLE ISSUES
 - la glmak usa gpu problema con macchine virtuali, per 



See the file `COPYING` for license information.

