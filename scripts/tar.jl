using Pkg
Pkg.activate("./")
Pkg.instantiate()
Pkg.precompile()

using PorpidPostproc, NextGenSeqUtils, FASTX

# zip Dataset directory for easy download

function my_read_fasta_records(filename)
    stream = open(FASTA.Reader, filename)
    records = FASTA.Record[]
    for entry in stream
        push!(records, entry)
    end
    return records
end

function read_fasta_with_names_and_descriptions(filename; seqtype=String)
    records = my_read_fasta_records(filename)
    return [FASTA.identifier(r) for r in records],
            [FASTA.description(r) for r in records],
             seqtype[FASTA.sequence(seqtype, r) for r in records]
end

function write_fasta_with_names_and_descripts(filename::String, seqs; names = String[], descripts = String[])
    if length(names) > 0 && (length(names) != length(seqs) || length(names) != length(descripts) )
        error("number of sequences does not match number of names")
    end
    if length(names) == 0
        names = ["seq_$i" for i in 1:length(seqs)]
        descripts = ["" for i in 1:length(seqs)]
    end
    stream = open(FASTA.Writer, filename)
    for (name, seq, descript) in zip(names, seqs, descripts)
        write(stream, FASTA.Record(name, descript, seq))
    end
    close(stream)
end

dataset = snakemake.wildcards["dataset"]
dataset_dir = "dataset/$(dataset)"
samples = snakemake.params["samples"]


# now rename dataset_dir, zip and rename back
run(`mv $(dataset_dir) $(dataset_dir)--sga-index`)
run(`tar -C dataset -czf dataset/$(dataset)--sga-index.tar.gz $(dataset)--sga-index`)
run(`mv $(dataset_dir)--sga-index $(dataset_dir)`)
println("dataset directory archived and zipped ...")

