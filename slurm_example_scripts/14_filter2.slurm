#!/bin/bash
#SBATCH --job-name="14_filter2"
#SBATCH --time=5:00:00  # CHANGE - walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # processor core(s) per node
#SBATCH --mail-user="rsschwartz@uri.edu" #CHANGE user email address
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_%x_%j.out
#SBATCH --error=slurm_%x_%j.err
#SBATCH -p uri-cpu

cd $SLURM_SUBMIT_DIR

module purge

#CHANGE IF NOT ON UNITY
module load uri/main bio/Biopython/1.79-intel-2022a
module load BLAST+/2.15.0-gompi-2023a

OUTFOLDER=../../SISRS_Small_test/ #CHANGE THIS to the analysis directory
MAX_CONTIGS=ALL #CHANGE to reduce output

../scripts/sisrs filter2 -d $OUTFOLDER -m $MAX_CONTIGS
