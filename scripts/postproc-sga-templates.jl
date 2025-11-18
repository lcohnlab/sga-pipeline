using Pkg
Pkg.activate("./")
Pkg.instantiate()
Pkg.precompile()

ENV["MPLBACKEND"] = "Agg"
using PorpidPostproc, CSV, NextGenSeqUtils, BioSequences, DataFrames, DataFramesMeta


# include("../../src/functions.jl")
# include("../../src/molev_functions.jl")
# include("../../src/postproc_functions.jl")


dir = snakemake.input[1]
fastq_dir = snakemake.input[2]
out = snakemake.output[1]
agreement_thresh = snakemake.params["agreement_thresh"]
max_alignment_reads = snakemake.params["max_alignment_reads"]
dataset = snakemake.wildcards["dataset"]
all_seqs = dir*"/$(dataset).fasta"

passed = out*"/passed_$(agreement_thresh)_agreement/"
failed = out*"/below_$(agreement_thresh)_agreement/"

"""
function that performs minimum agreement filtering of SGA sequences
and separates them into different folders    
  
"""
function sga_template_proc(dir, out; 
	  agreement_thresh=0.7)
    files = [dir*"/"*f for f in readdir(dir) if f[end-5:end] == ".fasta"]
    files = filter!(x -> x != all_seqs, files) #remove file with all sequences from list
    if length(files) > 0
    	println("")
        println("Filtering sequences by $(agreement_thresh) minimum agreement...")
        println("")
    else
        println("WARNING: no sequences found")
        exit()
    end
    
    #filter sequences by minimum agreement  
    for f in files 
    	seqname,seq = read_fasta_with_descriptors_in_names(f) #some update broke this and sequence names became duplicated, see fix below
    	seqname2 = [split(seqname[1], " "; limit=2)[2]] #remove duplication of sequence name  
		filename = replace.(seqname[1], r" .*" => "")
    	min_ag = replace.(seqname[1], r".*min_agreement=" => "")
    	min_ag = parse.(Float64, min_ag) #convert to number
    	if min_ag >= agreement_thresh
    		write_fasta(passed*"$(filename).fasta", seq, names = seqname2)
    	else
    		write_fasta(failed*"$(filename).fasta", seq, names = seqname2)
    	end		
    end
end
	






"""
function that creates alignments for all samples that fail min agreement filtering
alignments are only created for samples with fewer than the max_alignment_reads value to prevent extremely long runtimes for samples with many reads
"""
function alignment_for_failed_filter(fasta_dir, fastq_dir; max_alignment_reads = 1000)
    fasta_files = [fasta_dir*f for f in readdir(fasta_dir) if f[end-5:end] == ".fasta"]
    if length(fasta_files) > 0
        println("Creating alignments from samples below minimum agreement threshold and with fewer than $(max_alignment_reads) reads...")
        println("")
    end
    fastq_files = replace.(fasta_files, r"filtered/below_0.7_agreement" => "demux") #modify file names to read the demux fastq files
    fastq_files = replace.(fastq_files, r"fasta" => "fastq")
    
    t1 = time()
    #create alignments
    for f in fastq_files
    	t3 = time()
    	sample = split(basename(f),".fast")[1] 
    	seqs,phreds,seq_names = read_fastq(f)
    	println("$(sample) - $(length(seqs)) reads")
    	if length(seqs) < max_alignment_reads
    		aligned_seqs = mafft_align(seqs)
    		filename = replace(f, r".fastq" => "_alignment.fasta") 
    		filename = replace(filename, r"demux" => "filtered/below_0.7_agreement/alignments") #change filename back to the filtered directory
    		write_fasta(filename, aligned_seqs, names = seq_names)
    		t4 = time()
    		#println("$(t4 - t3) seconds.")
    	else println("$(sample) has more than the maximum number of allowed reads ($(max_alignment_reads)), no alignment will be created. This value can be adjusted in the snakefile.")
    		filename2 = replace(f, r".fastq" => ".fasta") 
    		filename2 = replace(filename2, r"demux" => "filtered/below_0.7_agreement/alignments") #change filename back to the filtered directory
    		write_fasta(filename2, seqs, names = seq_names)
    	end	
    end
    t2 = time()
    println("")
    println("Creating all alignments took $(t2 - t1) seconds.") 
    println("") 
end


mkpath(out)
mkpath(passed)
mkpath(failed)
mkpath(failed*"/alignments")

#filter sga sequences by minimum agreement
sga_template_proc(dir, out, agreement_thresh=agreement_thresh)

#create alignments for all sequences below min agreement threshold
alignment_for_failed_filter(failed,fastq_dir, max_alignment_reads = max_alignment_reads)
