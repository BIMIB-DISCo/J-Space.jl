### -*- Mode: Julia -*-

### CallART.jl

####################### CALL ART
#function call ART but only tecnology illumina

"Function call_ART

"
function call_ART(len_read::Int, mode::String; path = "")
    cd("Fasta output\\")        # Cambio directory.
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
    cd("..\\")
end

### end of file -- CallART.jl

