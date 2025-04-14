using Pkg
Pkg.activate("./")
Pkg.instantiate()
Pkg.precompile()

using Distributed
@everywhere ENV["MPLBACKEND"] = "Agg"
@everywhere using RobustAmpliconDenoising, NextGenSeqUtils

function generateConsensusFromDir(dir, template_name)
    files = [dir*"/"*f for f in readdir(dir) if f[end-5:end] == ".fastq" && f[end-9:end] != "QUAL.fastq" && f[end-10:end] != "issue.fastq"]
    if length(files) == 0
        println("WARNING: no read collections found")
        exit()
    end
    keeps = String[]
    for f in files #select only files with read number >= min_reads
    	sample = split(basename(f),".fast")[1]
    	seqs,phreds,seq_names = read_fastq(f)
    	if length(seqs) >= min_reads
    		push!(keeps, f)
    	else
    		println("No consensus for $(sample): # reads below cutoff of $(min_reads)")	
    	end
    end		
    println("")
    println("Generating consensus for $(length(keeps)) samples")
    println("")
    cons_collection = pmap(ConsensusFromFastq, keeps)
    seq_collection = [i[1] for i in cons_collection]
    seqname_collection = [template_name*i[2] for i in cons_collection]
    return seq_collection, seqname_collection
end

@everywhere function ConsensusFromFastq(file)
    seqs,phreds,seq_names = read_fastq(file)
    sample = split(basename(file),".fast")[1]
    println("$(sample)")
    draft = consensus_seq(seqs)
    draft2 = refine_ref(draft, seqs)
    final_cons = refine_ref(draft2,seqs)
    alignments, maps, matches, matchContent = getReadMatches(final_cons, seqs, 0)
    cons_name = split(basename(file),".fast")[1]*" num_CCS=$(length(seqs)) min_agreement=$(round(minimum(matches); digits = 2))"
    return final_cons, cons_name
end

@everywhere begin
"""
Returns an array of degapped coordinates, such that
coords(ref, read)[i] gives you the position the aligned read/ref
that matches the i'th ungapped position in ref.
"""
function coords(ref, read)
    if length(ref) != length(read)
        error("Aligned strings are meant to be the same length.")
    end
    degappedRef = degap(ref)
    coordMap = zeros(Int64, length(degappedRef))
    count = 1
    for i in 1:length(degappedRef)
        while ref[count] == '-'
            count += 1
        end
        coordMap[i] = count
        count += 1
    end
    return coordMap
end
end

@everywhere begin
"""
Return matches to a candidate reference from a set of reads.
"""
function getReadMatches(candidate_ref, reads, shift; degap_param = true, kmer_align = true)
    alignments = []
    if kmer_align
        alignments = map(i -> kmer_seeded_align(candidate_ref, i), reads)
    else
        alignments = map(i -> nw_align(candidate_ref, i), reads)
    end

    maps = [coords(i...) for i in alignments]

    if (degap_param)
        matchContent = [[degap(alignments[i][2][maps[i][k]:maps[i][k+shift]]) for i in 1:length(maps)] for k in 1:length(candidate_ref)-shift]
        matches = [freq(matchContent[k], degap(candidate_ref[k:k+shift])) for k in 1:length(matchContent)]
    else
       matchContent = [[(alignments[i][2][maps[i][k]:maps[i][k+shift]]) for i in 1:length(maps)] for k in 1:length(candidate_ref)-shift]
       matches = [freq(matchContent[k], candidate_ref[k:k+shift]) for k in 1:length(matchContent)]
    end
    return alignments, maps, matches, matchContent
end
end

function write_ind_files(cons, name)
	seq = String[] #had to create a vector for the string or write_fasta returns length error
	push!(seq, cons)
	seq_name = String[]
	push!(seq_name, name)
	samplename = split(basename(name),' ')[1]
	write_fasta(snakemake.output[1]*"/$(samplename).fasta",
	seq,
	names = seq_name)
end

#Calculate consensus sequences for each fastq file.
t1 = time()
dataset = snakemake.wildcards["dataset"]
base_dir = snakemake.input[1]
all_seqs = snakemake.output[1]*"/$(dataset).fasta"
min_reads = snakemake.params["min_reads"]
println("The min reads is: $(min_reads)")
cons_collection, name_collection = generateConsensusFromDir(base_dir, "")

#write out individual sequence files and a collection of all sequences
mkdir(snakemake.output[1])
map(write_ind_files, cons_collection, name_collection)
write_fasta(all_seqs,
    cons_collection;
    names = name_collection)
t2 = time()
println("")
println("Consensus generation for all samples took $(t2-t1) seconds.")
println("")

