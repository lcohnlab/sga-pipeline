## branch sga-pipeline

This branch has been adapted by Dylan Westfall and Anna Farrell-Sherman
to enable the pipeline to demultiplex SGA samples labeled with Index
primers rather than unique Sample IDs present in the cDNA primer. This
primarily revolved around replacing the demultiplexing code in
demux_functions.jl with code which Alec Pankow had written previously
for an older version of PORPIDpipeline that could accept Indexed samples
and removing code which handled UMIs and postprocessing reports.

Based upon code from PORPIDpipeline by Alec Pankow and Ben Murrell, now
maintained by Hugh Murrell.
(https://github.com/MurrellGroup/PORPIDpipeline.git)

now upgraded to Julia version 1.10.5

## Quick start


### Dependencies (on an ubuntu machine)

- first update all apps
   - `apt update`
   - `apt upgrade`
- Snakemake
   - `apt-get install -y snakemake`
- python3 packages
  - `apt-get install python3-pandas`
  - `apt-get install python3-seaborn`


### Julia version 1.10.5

We recommend you use the `juliaup` version manager to install julia.
from a terminal you can do this as follows:

```bash
curl -fsSL https://install.julialang.org | sh
```

This should install the Julia version manager, `juliaup` as well as
the latest version of Julia. To find out how to use the version manager 
to makesure you have version 1.10.5 as your default, go here:

[https://github.com/JuliaLang/juliaup]

Once Julia is installed, make sure you can enter the julia REPL from 
the command line and check the version number by logging out and in again 

```bash
exit
ssh root@.....
```

and then from your new terminal session:

```bash
juliaup status
julia --version
```

If the version number is not 1.10.5 then you need to use `juliaup` to install
that version and make it the default. 

```bash
juliaup add 1.10.5
juliaup default 1.10.5
```

for further details concerning `juliaup` go here:

[https://github.com/JuliaLang/juliaup?tab=readme-ov-file#using-juliaup]

### cloning the PORPIDpipeline repository

Now that the dependencies are setup we clone the repository

```bash
cd ~
git clone https://github.com/lcohnlab/sga-pipeline.git
```

### setting up the Julia package environment

then navigate to the `PORPIDpipeline` project folder and start the Julia REPL. 
Enter the package manager using `]` and then enter

```julia
activate .
instantiate
precompile
```

This will activate, install, and precompile the `julia` environment specified by the 
`Project.toml` and `Manifest.toml` files. The `precompile` command
above is not strictly needed but is useful if there are issues with installing
the `julia` packages listed in `Project.toml`

### Configuration

To configure the workflow, use the format below to define
your library construction. 
It should follow the same format as shown below:

```yaml
Dataset1:
  Sample1:
    fwd_index: Index_F21
    index_type: Index_primer
    rev_index: Index_R02
    rev_primer: CAAGGCAAGCTTTATTGAGGCTTAAS
    sec_str_primer: GATTGTGTGGCARGTAGACAGRATG
  Sample2:
    fwd_index: Index_F22
    index_type: Index_primer
    rev_index: Index_R02
    rev_primer: CAAGGCAAGCTTTATTGAGGCTTAAS
    sec_str_primer: GATTGTGTGGCARGTAGACAGRATG
  Sample3:
    fwd_index: Index_F23
    index_type: Index_primer
    rev_index: Index_R02
    rev_primer: CAAGGCAAGCTTTATTGAGGCTTAAS
    sec_str_primer: GATTGTGTGGCARGTAGACAGRATG
```

The primer sequences provided will be used for demultiplexing and will be trimmed
from the final sequences. **sec_str_primer** and **rev_primer** are the PCR primers 
used in the 2nd rd of PCR. **fwd_index** and **rev_index** are the index primers used to 
tag each sample. Index primer set is defined in the file CohnNexteraIndexes.csv.

Raw CCS .fastq files should be placed in the `raw-reads/` subdirectory and named 
according to the the dataset name used in the `config.yaml` file, ie, `Dataset1.fastq`
for the example dataset.

The contamination check will compare all sequences to one another, but if it is desired for
another sequence (for example a lab strain) to be included any sequence can be added to the
`contam_panel.fasta` file in the panel folder. The sequence must be trimmed to the same length
as the target sequence it is to be compared to. Note that the contamination check will flag
sequences that are within the distance threshold `dist_thresh` specified in the snakefile
(default 1.5%).


### Preview and Execution

Preview jobs with Snakemake and run with {n} cores.

```bash
#preview jobs
snakemake -np

#run
snakemake -j{n}
```

For more info on Snakemake, see https://snakemake.readthedocs.io/en/stable/



