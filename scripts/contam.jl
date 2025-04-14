using Pkg
Pkg.activate("./")
Pkg.instantiate()
Pkg.precompile()

ENV["MPLBACKEND"] = "Agg"
using PorpidPostproc, NextGenSeqUtils, BioSequences, DataFrames,CSV, DataFramesMeta
# include("../../src/functions.jl")
# include("../../src/contam-filter_functions.jl")

proportion_thresh = snakemake.params["proportion_thresh"]
cluster_thresh = snakemake.params["cluster_thresh"]
dist_thresh = snakemake.params["dist_thresh"]

dir = snakemake.input[1]
run_ID = snakemake.wildcards["dataset"]
#contam_passed_dir = snakemake.output[1] #these not needed for sga sequences as there is only one sequences so all will pass as they can't be different from self consensus
#contam_failed_dir = snakemake.output[2]
#mkpath(contam_passed_dir)
#mkpath(contam_failed_dir)


files = [dir*"/"*f for f in readdir(dir) if f[end-5:end] == ".fasta" && f != "$(run_ID).fasta"]
    if length(files) > 0
    	println("")
        println("Performing contamination check...")
    else
    	println("")
        println("WARNING: no template families found")
        exit()
    end

#files = snakemake.input["files"]
one_run_db = vcat([db_cluster_seqs(f,
                        proportion_thresh = proportion_thresh,
                        cluster_thresh = cluster_thresh) for f in files]...);
contam_db = db_seqs(snakemake.input["panel"])
merged_db = vcat(one_run_db,contam_db);


contam_df = DataFrame(  sample = String[],
                        sequence_name = String[],
                        nearest_nonself_variant = String[],
                        nearest_nonself_distance = Float64[],
                        flagged = Bool[],
                        discarded = Bool[]);
for f in files
    sample = split(basename(f),".fast")[1] #modify this if you want to get the dataset without the file extension.
    flagged_names,flagged_neighbours,flagged_dists,reported,discarded,discarded_bool_inds = contam_check(f,merged_db, thresh=dist_thresh)
    seqnames,seqs = read_fasta_with_descriptors_in_names(f)
    if length(flagged_names) > 0
        for r in 1:length(flagged_names)
            #println((flagged_names,flagged_neighbours,flagged_dists,reported,discarded))
            push!(contam_df,[sample,flagged_names[r],basename(flagged_neighbours[r]),flagged_dists[r],reported[r],discarded[r]])
        end
        #write_fasta(contam_passed_dir*"/"*sample*".fasta",seqs[.!discarded_bool_inds],names = seqnames[.!discarded_bool_inds]) #no need to write these files when using SGA
        #write_fasta(contam_failed_dir*"/"*sample*".fasta",seqs[discarded_bool_inds],names = seqnames[discarded_bool_inds])
    else
        #write_fasta(contam_passed_dir*"/"*sample*".fasta",seqs,names = seqnames)
    end
end
CSV.write(snakemake.output[1],contam_df) #changed to output1 after removing the passed/failed directories
println("")
