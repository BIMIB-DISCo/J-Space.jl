### -*- Mode: Julia -*-

### Bulk_Experiment.jl

########################## Experiments
#function calculate bulk experiment

"Function create_bulk_groundtruth

"
function create_bulk_groundtruth(G_seq::LongDNASeq,
                                 Fasta_sample::Vector{Any},
                                 Position::Vector{Any})
    # Dataframe -> POS, REF_nucleotide, %A, %C, %G, %T
    #=
    df = DataFrame(Position = Int[], Reference = DNA[],  A = Any[], C = Any[],
                   G = Any[], T = Any[])
    positions = []
    n_sample = length(Fasta_sample)
    for i = 1:len_ROI
        nucl_ref = G_seq[i]
        nucleotides = [f[i] for f in Fasta_sample]
        nucl_counter = counter(nucleotides)
        if nucl_counter[nucl_ref] < n_sample
            pos = i
            k_dict = keys(nucl_counter.map)
            values = []
            for n in names(df)[2:end]
                if n != string(nucl_ref)
                    n_dna = DNA(collect(n)[1])
                    push!(values, (nucl_counter[n_dna]*100)/n_sample)
                else
                    push!(values, "ref")
                end
            end
            push!(df, [pos, nucl_ref, values[1], values[2], values[3], values[4]])
        end
    end
    =#

    df = DataFrame(MUT = Any[], VAF = AbstractFloat[])
    n_sample = length(Fasta_sample)
    Position_sort = sort(Position)
    for pos in Position_sort
        nucl_ref = G_seq[pos]
        nucleotides = [f[pos] for f in Fasta_sample]
        nucl_counter = counter(nucleotides)
        k_dict = keys(nucl_counter.map)
        for k in k_dict
            # println("k :",k,"\t",typeof(k))
            # println("k != nucl_ref -> ",k != nucl_ref)
            if k != nucl_ref
                mut = string(pos) * "_" * string(nucl_ref) * "->" * string(k)
                vaf = nucl_counter[k]/n_sample
                # println("vaf -> ",vaf)
                push!(df, [mut, vaf])
            end
        end
    end
    return df
end


"Function bulk_with_noise

"
function bulk_with_noise(Coverage, VAF, n_sample, FP, FN, Positions, seed)
    df_noise = DataFrame(MUT = Any[], VAF = AbstractFloat[])
    for i in 1:length(Positions)
        n_reads = rand(seed, Poisson(Coverage), 1)[1]
        # println("n_reads: ", n_reads)
        n_read_muts = rand(seed, Binomial(n_reads, VAF[i, :VAF]), 1)[1]
        # println("n_read_muts: ", n_read_muts)
        n_FP = rand(seed, Binomial(n_sample-n_read_muts, FP), 1)[1]
        # println("n_FP: ",n_FP)
        n_FN = rand(seed, Binomial(n_read_muts, FN), 1)[1]
        # println("n_FN: ",n_FN)
        n_VAF_Noise = (n_read_muts + n_FP - n_FN)/n_sample
        # println("n_VAF_Noise: ", n_VAF_Noise)
        push!(df_noise, [VAF[i, :MUT], n_VAF_Noise])
    end
    return df_noise
end


"Function experiment_bulk

"
function experiment_bulk(reference::LongDNASeq,
                         fasta_samples::Vector{Any},
                         position_used::Vector{Any},
                         seed::MersenneTwister;
                         Noise::Int = 0,
                         coverage::AbstractFloat = 0,
                         FP::AbstractFloat = 0,
                         FN::AbstractFloat = 0)
    # calculate bulk experiment if bulk is true
    VAF_GT = create_bulk_groundtruth(reference, fasta_samples, position_used)
    CSV.write("VAF_GT.vcf", VAF_GT, delim = ",")
    if Noise == 1
        if coverage == 0 && FP == 0 && FN == 0
            return "Error: insert covarege, FP and FN"
        end
        n_sample = length(fasta_samples)
        VAF_GT_Noise = bulk_with_noise(coverage,
                                       VAF_GT,
                                       n_sample,
                                       FP,
                                       FN, position_used, seed)
        CSV.write("VAF_GT_Noise.vcf", VAF_GT_Noise, delim = ",")
    end
    return VAF_GT, VAF_GT_Noise #, g_seq
end


### end of file -- Bulk_Experiment.jl
