#!/bin/bash
#SBATCH --job-name="11_contig_processing"
#SBATCH --time=4:00:00  # CHANGE walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=10   # CHANGE for processor core(s) per node
#SBATCH --mail-user="rsschwartz@uri.edu" #CHANGE user email address
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_%x_%j.out
#SBATCH --error=slurm_%x_%j.err
#SBATCH -p uri-cpu

cd $SLURM_SUBMIT_DIR

module purge

#CHANGE IF NOT ON UNITY
module load uri/main
module load BCFtools/1.17-GCC-12.2.0
module load MAFFT/7.505-GCC-11.3.0-with-extensions

module load conda/latest
conda activate sisrs

#Optionally specify coverage threshold: '-c', '--cov' (default set to 3) and heterozygosity threshold: '-z', '--hz'(default set to 0.01)

T=4 #CHANGE THIS to the threshold number
P=$SLURM_JOB_CPUS_PER_NODE

OUTFOLDER=../../SISRS_Small_test2 #CHANGE to analysis directory

../scripts/sisrs contig_processing -t $T -p $P -d $OUTFOLDER
