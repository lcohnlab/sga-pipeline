using Pkg
Pkg.activate("./")
Pkg.instantiate()
Pkg.precompile()

ENV["MPLBACKEND"] = "Agg"
using PorpidPostproc, NextGenSeqUtils, StatsBase
using BioSequences, FASTX
using DataFrames, CSV
using CodecZlib: GzipDecompressorStream

println("using Julia version: $(VERSION)")

t1 = time()

SAMPLE_CONFIGS = snakemake.params["config"]
samples = snakemake.params["samples"]
index_type = snakemake.params["index_type"]
mkdir(snakemake.output[1])


f_kwargs = [
    :demux_dir => snakemake.output[1],
    :samples => SAMPLE_CONFIGS,
    :verbose => false,
    :index_type => snakemake.params["index_type"],
    :error_rate => snakemake.params["error_rate"],
    :min_length => snakemake.params["min_length"],
    :max_length => snakemake.params["max_length"],
    :min_reads => snakemake.params["min_reads"],
    :label_prefix => "seq",
    :error_out => true
    ]

println("performing chunked quality filtering and demux on $(snakemake.input[1])")
chunk_size = snakemake.params["chunk_size"]
min_reads = snakemake.params["min_reads"]
println("")
println("Index type: $(index_type)")
println("Minimum read length $(snakemake.params["min_length"])")
println("Maximum read length $(snakemake.params["max_length"])")

reads = chunked_filter_apply(snakemake.input[1], chunked_quality_demux; chunk_size=chunk_size, f_kwargs, verbose = false)

# now do some accounting
total_reads = reads[1]
quality_reads = reads[2]
demuxed_reads = reads[3]

# report on total sequences for each sample
println()
println("total reads => $(total_reads)")
println("quality reads => $(quality_reads)")
println("-------------------------------------")
println("demuxed reads:")
no_assigned = 0
no_rejected = 0
filepaths = readdir(snakemake.output[1],join=true)
df_demux = DataFrame(Sample = [], Count = [])
for s in samples #add samples to df
	push!(df_demux,[s,0])
end
df_reject = DataFrame(Sample = [], Count = [])
for path in filepaths
    stream = open(FASTQ.Reader, path)
    count = 0
    for record in stream
        count += 1
        if occursin("REJECT",path)
            global no_rejected += 1
        else
            global no_assigned += 1
        end
    end
    sample_name = replace(basename(path), ".fastq" => "")
    #println(sample_name," => ",count)
    if occursin("REJECT",path)
        push!(df_reject,[sample_name,count])
    else
        push!(df_demux,[sample_name,count]) #this adds duplicate lines to df for each sample
    end
end

#select only one line of df_demux for each sample with greatest read value. This allows samples with no reads to show
df_demux = combine(df_demux -> df_demux[argmax(df_demux.Count), :], groupby(df_demux, :Sample))

println("$(df_demux)")
println("-------------------------------------")
println("total demuxed => $(no_assigned)")
println("total rejected => $(no_rejected)")
println("")
println("If read counts are lower than expected check that min and max read lengths match amplicon size")
println("")

# down sample reads to specified max
df_demux_sampled = DataFrame(Sample = [], Count = [], Sampled = [])
max_reads = snakemake.params["max_reads"]
if max_reads < 1
    max_reads = 10000000
end
println("downsampling to about $(max_reads) reads")
global no_lost=0
global no_retained=0
for path in filepaths
    out_path = path[1:end-9]*"_sampled.fastq.gz"
    if ! occursin("REJECT",path)
        if endswith(path, ".gz")
            stream = FASTQ.Reader(GzipDecompressorStream(open(path)))
            out_stream = FASTQ.Writer(GzipCompressorStream(open(out_path,"w")))
        else
            stream = FASTQ.Reader(open(path))
            out_stream = FASTQ.Writer(open(out_path,"w"))
        end
        sample_name = replace( replace( basename(path), ".gz" => "" ), ".fastq" => "" )
        count = 0
        sampled=0
        ac = df_demux[df_demux.Sample .== sample_name, :Count][1]
        sp = 1.0
        ac > 0 ? sp = min(1.0, max_reads / ac) : sp = 1.0
        if sp < 1.0
            for record in stream
                count += 1
                if rand() < sp
                    write(out_stream, record)
                    sampled += 1
                    global no_retained += 1
                else
                    global no_lost += 1
                end
            end
            close(stream)
            close(out_stream)
            println("$(sample_name), sampled $(sampled) reads")
            rm(path)
            mv(out_path,path)
            push!(df_demux_sampled,[sample_name,count,sampled])
        else
            close(stream)
            close(out_stream)
            println("$(sample_name), retained all $(ac) reads")
            global no_retained += ac
            rm(out_path)
            push!(df_demux_sampled,[sample_name,ac,ac])
        end
    end
end

CSV.write("$(snakemake.output[3])", df_demux)
CSV.write("$(snakemake.output[4])", df_reject)

# save a quality_filter report
# no_of_fails = no_of_reads - no_assigned # fix this, get chunked_filter_apply to return no_of_fails
df_qual = DataFrame(Description = [], Count = [])
push!(df_qual,["Initial number of raw reads",total_reads])
push!(df_qual,["Number of quality reads",quality_reads])
push!(df_qual,["Number assigned by demux",no_assigned])
CSV.write("$(snakemake.output[2])", df_qual)

t2 = time()
println("")
println("Quality filtering and demultiplexing took $(t2 - t1) seconds.")
