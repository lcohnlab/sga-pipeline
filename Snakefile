import subprocess, sys
configfile: "config.yaml"
DATASETS = [d for d in config for s in config[d]]
SAMPLES = [s for d in config for s in config[d]]
VERSION = "1.10.5"
COMMIT = subprocess.check_output(['git', 'rev-parse', '--verify', 'HEAD']).strip().decode()

sys.stderr.write("Running PorpidPostproc\n")
sys.stderr.write("Version: {0}\n".format(VERSION))
sys.stderr.write("Commit ID: {0}\n".format(COMMIT))

# sga-index-consensus parameters
# demux
chunk_size = 100000  		# default 100000
index_type = "Index_primer" # default "Index_Primer", also accepts "Nextera_primer"
error_rate = 0.01    		# default 0.01
min_length = 2000    		# default 2100
max_length = 6000    		# default 4000
max_reads = 1000       		# default 1000 reads per sample,
                         	# use something large for no downsampling
# consensus
min_reads = 5				# default 5
# contam
cluster_thresh = 0.015   	# default 0.015
proportion_thresh = 0.2  	# default 0.2
dist_thresh = 0.015      	# default 0.015
# postproc
agreement_thresh = 0.7   	# default 0.7
max_alignment_reads = 1100  # default 1100
                            # be sure value is ~10% larger than max_reads to avoid conflicts

rule all:
    input:
        expand("dataset/{dataset}--sga-index.tar.gz",
            dataset = DATASETS)

rule demux:
    input:
        "raw-reads/{dataset}.fastq.gz"
    output:
        directory("dataset/{dataset}/demux"),
        "dataset/{dataset}/quality_report.csv",
        "dataset/{dataset}/demux_report.csv",
        "dataset/{dataset}/reject_report.csv"
    params:
        samples = SAMPLES,
        chunk_size = chunk_size,
        error_rate = error_rate,
        min_length = min_length,
        max_length = max_length,
        max_reads = max_reads,
        min_reads = min_reads,
        index_type = index_type,
        config = lambda wc: config[wc.dataset]
    script:
        "scripts/demux.jl"

rule consensus:
    input:
    	"dataset/{dataset}/demux"
    output:
        directory("dataset/{dataset}/consensus")
    params:
        min_reads = min_reads,
        config = lambda wc: config[wc.dataset]
    script:
        "scripts/consensus.jl"

rule contam:
    input:
        "dataset/{dataset}/consensus",
        panel = "panels/contam_panel.fasta"
    output:
        "dataset/{dataset}/contam_report.csv"
    params:
        proportion_thresh = proportion_thresh,
        cluster_thresh = cluster_thresh,
        dist_thresh = dist_thresh,
        config = lambda wc: config[wc.dataset]
    script:
        "scripts/contam.jl"

rule postproc:
    input:
        "dataset/{dataset}/consensus",
        "dataset/{dataset}/demux",
        "dataset/{dataset}/contam_report.csv"
    output:
        directory("dataset/{dataset}/filtered")
    params:
        agreement_thresh = agreement_thresh,
        max_alignment_reads = max_alignment_reads
    script:
        "scripts/postproc-sga-templates.jl"
                          
rule tar:
    input:
        "dataset/{dataset}/filtered"
    output:
        "dataset/{dataset}--sga-index.tar.gz"
    params:
        datasets = DATASETS,
        samples = SAMPLES
    script:
        "scripts/tar.jl"

