### -*- Mode: Julia -*-

### CallART.jl

####################### CALL ART
#function call ART but only tecnology illumina

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

"""Function call_ART
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
"""
function check_input(paired_end, mate_pair, mean_fragsize, std_fragsize)
    control = true
    if paired_end == mate_pair
        control = false
    end
    if (mean_fragsize == 0 || std_fragsize == 0) && paired_end == 1
        control = false
    end

    #if single-end => file, outputfile prefix, read length, count (-c)

    #if paired-end => come single-end, mean_fragsize, std_fragsize
    #if paired_end == 1(-> true), then mate_pair != 1 and vice versa
    return control
end

function call_ART(command::String, path_fasta::String ,path_fileout::String)
    com = split(command)
    ind_ss = findall(x -> x == "-ss", com)
    next = 0
    if typeof(ind_ss) != Int
        ind_ss = 0
    end
    current_dir = pwd()
    if Sys.islinux()
        for file in readdir(path_fasta)
            f = hcat(split.(file, ".")...)[1, :]
            if length(f) > 1 && f[2] == "fasta"
                path_ref = path_fasta*file
                i = ["-i", path_ref]
                new_com = com[1:ind_ss+1]
                append!(new_com, i)
                append!(new_com, com[ind_ss+2:end])
                cd(path_fileout)
                com_run = `$new_com`
                println("command: ", com_run)
                run(com_run)
                cd(current_dir)
            end
        end
    else
        for file in readdir(path_fasta)
            f = hcat(split.(file, ".")...)[1, :]
            if length(f) > 1 && f[2] == "fasta"
                path_ref = path_fasta*file
                i = ["-i", path_ref]
                new_com = com[1:ind_ss+1]
                append!(new_com, i)
                append!(new_com, com[ind_ss+2:end])
                cd(path_fileout)
                com_run = `$new_com`
                println("command: ", com_run)
                run(com_run)
                cd(current_dir)
            end
        end
    end
end
### end of file -- CallART.jl


function call_ART(profile::String,
                  path_fasta::String,
                  path_fileout::String,
                  len_read::Int,
                  tot_num_reads::Int,
                  outfile_prefix::String, #this is unless
                  paired_end::Bool,
                  seed::MersenneTwister;
                  sam::Bool = false,
                  ef::Bool = true,
                  mate_pair::Bool = false,
                  mean_fragsize::Int = 0,
                  std_fragsize::Int = 0,
                  no_ALN::Bool = false)

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

    ss = ["-ss", profile]
    l = ["-l", len_read]
    c = ["-c", tot_num_reads]
    ef_c = sam_c = p = mp = na = m_s = []

    if ef
        ef_c = ["-ef"]
    elseif sam
        sam_c = ["-sam"]
    end

    if paired_end
        p = ["-p"]
        mean_std = true
    end

    if mate_pair
        mp = ["-mp"]
        mean_std = true
    end

    if no_ALN
        na = ["-na"]
    end

    if mean_std
        m_s = ["-m", mean_fragsize, "-s", std_fragsize]
    end
    current_dir = pwd()
    if Sys.iswindows()
        for file in readdir(path_fasta)       # Scorro tutti i file
            f = hcat(split.(file, ".")...)[1, :]
            if length(f) > 1 && f[2] == "fasta" && f[1] != "reference"
                path_ref = path_fasta*file
                i = ["-i", path_ref]
                outfile_prefix = f[1]*"_"
                o = ["-o", outfile_prefix]
                cd(path_fileout)
                command = `art_illumina $ss $ef_c $p $mp $sam_c $na $i $l $c $m_s $o`
                println("command: ", command)
                run(command)
                cd(current_dir)
            end
        end
    else
        for file in readdir(path_fasta)
            f = hcat(split.(file, ".")...)[1, :]
            if length(f) > 1 && f[2] == "fasta" && f[1] != "reference"
                path_ref = path_fasta*file
                i = ["-i", path_ref]
                outfile_prefix = f[1]*"_"
                o = ["-o", outfile_prefix]
                cd(path_fileout)
                command = `art_illumina $ss $ef_c $p $mp $sam_c $na $i $l $c $m_s $o`
                println("command: ", command)
                run(command)
                cd(current_dir)
            end
        end
    end
end
