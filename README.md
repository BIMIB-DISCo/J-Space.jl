# J-SPACE 
## INTRODUCTION 
J-SPACE is a Julia package to simulate the genomic evolution and the spatial growth of a cell population and sequencing the genome of the sampled cells.
.
Firstly, the software simulates the spatial dynamics of the cells as a continuous-time multi-type birth-death stochastic process on a graph employing different rules of interaction and an optimised Gillespie algorithm. 
After mimicking a spatial sampling of the tumour cells,
J-SPACE  returns the phylogenetic tree of the sample and simulates molecular evolution of the genome under the infinite-site models or a set of different substitution models with the possibility of including structural variants as the indels. Finally, employing ART, J-SPACE generates the  synthetic single-end, paired-end/mate-pair reads of the next-generation sequencing platforms.


 ### Spatial clonal dynamics
 In J-SPACE the dynamics of the spatio-temporal evolution of a tumour is modelled by a stochastic multi-type Birth-Death process over an
arbitrary graph. J-SPACE could generate by itself a 2D or 3D regular lattice. In addition, it is possible to give as input any graph as an adjacency matrix ( an example of the format needed for this matrix is given in path/). 
In this part, the user can tune the birth rate of wild type cells, the death rate of the cells, the rate of migration of cells (not tested), the rule of contact between cells, the probability to develop a driver mutation per division, and the average birth rate advantage of a driver mutation.
Note that every rate inserted in J-SPACE must have the same unit of time both for the spatial dynamics and the molecular evolution.  

 ### Molecular evolution
 
 
J-SPACE simulate the evolution of the sequence of the sample after the simulation of the clonal dynamics. The user can sample the whole population or a subset of it, and the J-SPACE evaluate the phylogenetic tree of the samples. This GT tree is returned as a Newick file. 

The molecular evolution of an ancestral genome  (which can be given by the user as FASTA file or generated randomly) is simulated along the sampled tree via the Doob-Gillespie algorithm.
The user can use an infinite-site model to have fast simulations of situations where the genome is long, the mutational rate is very low (e.g.,<10^-8 substitution for unit of time), and the total simulated time is long.

In the case of finite-site models, J-SPACE takes as input the matrix of instantaneous rates for different substitution models: JC69, F81, K80, HKY85, TN93, and K81.
We suppose that the indels have a size distributed as  a Lavalette distribution.
 Note that using finite-site for long genomes come at the cost of computational performance.
 After this computation, the sequences of the samples (i.e., the leafs of the phylogenetic tree) are returned as FASTA file in the folder xxxx.
 
 

 ### Sequencing experiment

To simulate the reads of a sequencing experiment J-SPACE calls ART (https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm). The user can use a configuration file to specify the error model (for Illumina platforms), the number of reads, the length of the reads, and if the experiment uses single-end paired-end/mate-pair reads.
In addition, in the configuration file, there is the option to insert custom "calls" for ART  with the possibility to use it in any possible configuration.
If the user simulate the experiment, \alg{} returns for each cell, the simulated reads as FASTQ file, the alignment map of the reads over the genome of the sampled cells in SAM and/or ALN format. 
If the infinite-site model is used,  it is possible to obtain the VCF file directly without simulating the reads with ART.
## REQUIRED  SOFTWARE AND PACKAGE



## INSTALLATION OF J-SPACE

## RUN J-SPACE

### RUN A SINGLE SIMULATION
The paramenters and the cofinguration of the simulation are managment by the user by modifing the files Parameters.tml and Config.tml that are detailed below.
### RUN A SINGLE FUNCTION


### RUN THE EXAMPLES

## OUTPUTS OF J-SPACE
 J-SPACE provides the following outputs. 

  -  The state of the lattice at any time of the simulation (as plot and metagraph).
  - The plot of clonal dynamics.
  -  The Ground Truth (GT) genotype of the sampled cells as FASTA files.
  - The GT phylogenetic  tree  of the sampled cells in Newick format.
  -   The mutational tree of the driver mutations (if present).
  -   The list of the driver mutation with their birth rate advantage (if present).
  -   The simulated next-generation sequencing reads as FASTQ files.
  -  The alignment file, which maps the noisy reads on the sequences of the sampled cells both in formats SAM and ALN.
  -   The alignment file, whithout noise in format SAM.



## THE CONFIGURATION FILE OF J-SPACE
In the file Config.tml the user can manage the configuration of J-SPACE. This file is useful to choose the path where save the files, the desired plot and outputfile.

The following are the paramenters of such file

- seed. The seed of the simulations.
- generate_reference.  I 0 the user should inser the reference genome in fasta format in the folder path_reference. If 1 J-SPACE generate a random sequence.
- path_reference.The path of the reference given by user.
- path_to_save_files, path of the folder where the output files should be saved
- path_to_save_plot,  path of the folder where the output plots should be saved
- Generate_graph, if 0 the user should inser the graph of the dynamics as an adjacency matrix. If 1 J-SPACE generate regular lattice, the paramenters of such lattice are specified in the fiele Paramenters.toml.
- Tree_Newick, if 1  the phylogenentic tree of the cells is saved as newick file.
- Final_configuration, if 1 return the metagraph (in format) of the final configuration of the lattice.
- Driver_list, if 1 returns the list of the driver mutations
- Driver_Tree, if 1 returns the tree of the driver mutations
- Single_cell_fasta, if 1 save the fasta of the GT sequences of the sampled cells 
- Single_cell_noise, if 0  the sequencing experiment is not performed, if 1 it returns the FASTQ  files of the reads and SAM files with noise  , if 2 it returns the FASTQ  files of the reads, the SAM files with noise and without noise.
- Alignment, if 1 it returns the GT aligment file in ALN format.
- VAF_GT, if 1 it returns the VAF of the sampled cells  (working only if isa is used).
- Bulk_noise , if 1 it returns an approximate bulk experiment not using ART (working only if isa is used).
- Time_of_sampling, insert an array of times. J-SPACE perform the plot of the state of the lattice in that times.
- Dynamic_Clonal_genotype, if 1 plots the dynamics of the clonal genotypes.
- Graph_configuration, if 1  returns the plot of the state of the lattice.


## THE PARAMETERS FILE OF J-SPACE
In the file Parameters.tml the user will find all the paramenters of the dynamics, molecular evolution and experiment. 



### Paramenters of the generation of the lattice

- row. Integer number of the rows of the regular lattice, not used if a graph is imported
- col. Integer number of the columns of the regular lattice, not used if a graph is imported
- dim. Integer number of the rows of the regular lattice, not used if a graph is imported. If =1 the simulation is 2D
- N_starting_cells. Integer  the number of starting cells.

- matrix_adjacency.  Path to the adjacency matrix ("path/matrix/adjacency") of the graph that should be imported.

### Paramenters of the clonal spatial dynamics

- Model. String, the possible values are ["contact", "voter", "hvoter"] they are the possible different  Models of interaction.
- Max_time. Real number, it is the maximum time of the simulation 
- rate_birth. Real number,  Birth rate per cell per unit of times
- rate_death. Real number,  Death rate per cell per unit of times
- rate_migration. Real number,  Migration rate per cell per unit of times
- drive_mut_rate. Real number, probability of generate a new driver (i.e., subclones) after a division event.
- average_driver_mut_rate. Real number,  average birth rate advantage a driver mutation 
- std_driver_mut_rate. Real number, standard deviation of the birth rate advantage of a driver mutation 

### Paramenters of the sampling
- Random_sampling. Integer, if 1  Random sampling, if 0 cirular/spherical sampling 
- num_cell. Integer, number of sampled cells

#### If Random_sampling = 0
-  pos_center. Integer number, it is the center of the sampling
- radius_sampling. Integer number, it is size radius of the sampling

### Paramenters of the molecular evolution
- length_genome.  Length of the ancestral genome. Used  if ancestral genome is not given. If ancestral genome is given equal to 0.
- type_isa. Integer number, if 1 J-SPACE use the ISA to simulate the molecular evolution, if 0 it is necessary to specify the substituion model below.
#### if type_isa = 1
- neut_mut_rate. Rate of mutation per site and per unit of time
#### if type_isa = 0
- sub_model. A string, variable that select the subistituion model, possible value ->[ "JC69","F81","K80", "HKY85","TN93","K81"]
- indel_size. Real number maximum size of indel 
- Lavelette_param. Real number the parameter of  the Lavelette distribution for the size of indels
- indel_rate. Rate of indel per site and per unit of time. To esclude indel put this parameter to 0.
- params. Rates of the substitution models
 if sub_model= "JC69" params =
 if sub_model= "F81" params =
 if sub_models= "K80" params =
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
## CONTACTS
Please feel free to contact us if you have problems running our tool at .
