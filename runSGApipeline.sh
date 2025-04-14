#!/bin/bash
#SBATCH --job-name==sga-index-consensus
#SBATCH --time=1:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=7

source /app/lmod/lmod/init/profile # to load modules

export TMPDIR='/fh/scratch/delete90/cohn_l/tree'

#snakemake also loads Python
module load snakemake/7.18.2-foss-2021b
module load MAFFT/7.490-GCC-11.2.0-with-extensions
module load FastTree/2.1.11-GCCcore-11.2.0
# module load Julia/1.7.3-linux-x86_64
module load Julia/1.10.5-linux-x86_64
module load matplotlib/3.5.2-foss-2021b
module load Seaborn/0.11.2-foss-2021b

if [ "$#" -lt 1 ]; then
    echo "please supply a config file name as first parameter"
    exit
fi
echo "config file is $1"

echo "${SLURM_JOB_NAME} job submited using ${SLURM_NTASKS} cores"

echo "Script executed from: ${PWD}"
BASEDIR=$(dirname $0)
echo "Script location: ${BASEDIR}"

# create a symbolic link for the snakemake config file to point to the config for the current study (Not necessary for FH due to storage permission)
#rm -f /tools/PORPIDpipeline/porpidpipeline/config.yaml
#ln -s /RAW/PORPID/CONFIG/$1.yaml /tools/PORPIDpipeline/porpidpipeline/config.yaml

# tell slurm where anaconda is and conda activate the PORPIDpipeline environment
#source /tools/PORPIDpipeline/anaconda3/etc/profile.d/conda.sh
#conda activate PORPIDpipeline

# navigate to the porpidpipeline directory and run snakemake
# add -F to to the snakemake command to force re-run of all rules
snakemake --rerun-incomplete -j${SLURM_NTASKS}  

echo "DONE"
