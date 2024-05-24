#!/bin/sh
#SBATCH --job-name="14_filter2"
#SBATCH --time=45:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # processor core(s) per node
#SBATCH --mail-user="rsschwartz@uri.edu" #CHANGE user email address
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_%x_%j.out
#SBATCH --error=slurm_%x_%j.err

cd $SLURM_SUBMIT_DIR

module purge

module load Biopython/1.75-foss-2019b-Python-3.7.4
module load BLAST+/2.8.1-foss-2018b

OUTFOLDER=../../SISRS_Small_test/ #CHANGE THIS to the analysis directory
MAX_CONTIGS=ALL #CHANGE to reduce output

../scripts/sisrs filter2 -d $OUTFOLDER -m $MAX_CONTIGS
