### -*- Mode: Julia -*-

### CallART.jl

####################### CALL ART
#function call ART but only tecnology illumina

"Function call_ART
 -ss The name of Illumina sequencing system of the built-in profile used
     for simulation
 -ef Indicate to generate the zero sequencing errors SAM file as well the
     regular one
     NOTE: the reads in the zero-error SAM file have the same alignment positions
 -i  The filename of input DNA/RNA reference
 -p  Indicate a paired-end read simulation or to generate reads from both ends
     of amplicons
	 NOTE: art will automatically switch to a mate-pair simulation if the given
           mean fragment size >= 2000
 -mp Indicate a mate-pair read simulation
 -l  The length of reads to be simulated
 -c  Number of reads/read pairs to be generated per sequence
 -m  The mean size of DNA/RNA fragments for paired-end simulations
 -s  The standard deviation of DNA/RNA fragment size for paired-end simulations
 -o  The prefix of output filename
 -na  Do not output ALN alignment file
 -rs The seed for random number generator (default: system time in second)
     NOTE: using a fixed seed to generate two identical datasets from different
           runs
"
#=input : Profile::String    =>  -ss $Profile
          ef::bool           =>  -ef
          path_ref::String   => -i $file
          paired_end::bool   => -p
          mate_pair::Bool    => -mp
          len_read::Int      => -l $len_read
          tot_num_reads::Int => -c $tot_num_reads
          mean_fragsize::Int => -m $mean_fragsize
          std_fragsize::Int  => -s $std_fragsize
          outfile_prefix::String => -o $outfile_prefix
          no_ALN::Bool           => -na
          seed::MersenneTwister => -rs $seed
=#

function check_input(paired_end, mate_pair, mean_fragsize, std_fragsize)
    control = true
    if paired_end == mate_pair
        control = false
    elseif mean_fragsize == 0 || std_fragsize == 0
        control = false
    end

    #if single-end => file, outputfile prefix, read length, count (-c)
    #sempre vero ?

    #if paired-end => come single-end, mean_fragsize, std_fragsize
    #if paired_end == 1(-> true) non posso avere mate_pair == 1 e viceversa
    return control
end

function call_ART(profile::String,
                  path_ref::String,
                  len_read::Int,
                  tot_num_reads::Int,
                  outfile_prefix::String,
                  paired_end::Bool,
                  seed::MersenneTwister;
                  ef::Bool = true,
                  mate_pair::Bool = false,
                  mean_fragsize::Int = 0,
                  std_fragsize::Int = 0,
                  no_ALN::Bool = 0)

    #=cd("Fasta output\\")        # Cambio directory.
    for file in readdir()       # Scorro tutti i file
        f = hcat(split.(file, ".")...)[1, :]
        if length(f) > 1 && f[2] == "fasta"
            mkpath(f[1])
            cd(f[1])
            command = `art_illumina -sam -i $file -l $len_read -ss HS25
                                                                -f 10 -o $mode`
            #run(command)
            cd("..\\")
        end
    end
    cd("..\\")=#
    if paired_end == true || mate_pair == true
        control = check_input(paired_end,
                              mate_pair,
                              mean_fragsize,
                              std_fragsize)
        if control == false
            return "Error with input file"
        end
    end

    #create command
    mean_std = false
    command = `sudo art_illumina `
    command = "sudo art_illumina "
    if ef
        command = command * "-ef"#`
        ef_c = ["-ef"]
    end
    if paired_end
        command = command * "-p "
        mean_std = true
    end
    if mate_pair
        command = command * "-mp "
        mean_std = true
    end
    if no_ALN
        command = command * "-na "
    end
    if mean_std
        command = command * "-ss $profile -i $file -l $len_read
                             -c $tot_num_reads -m  $mean_fragsize
                             -s $std_fragsize -o $outfile_prefix"
    else
        command = command * "-ss $profile -i $file -l $len_read
                             -c $tot_num_reads -o $outfile_prefix"
    end
    #Cmd(`sudo art_illumina $ef`) #cosi va
    #non funziona, devo trovare il modo
    #cmd = @cmd(command)
    #cmd = @cmd("art_illumina -sam -i file -l len_read -ss HS25 -f 10 -o mode")
end

function call_ART(command::Cmd#= più dove salvare le cose=#)
    #caso in cui ho tanti fasta -> da vedere
    #cd("Fasta output\\")        # Cambio directory in cui voglio salvare.
    for file in readdir()       # Scorro tutti i file
        f = hcat(split.(file, ".")...)[1, :]
        if length(f) > 1 && f[2] == "fasta"
            mkpath(f[1])
            cd(f[1])
            run(command)
            cd("..\\")
        end
    end
    cd("..\\")

end
### end of file -- CallART.jl
