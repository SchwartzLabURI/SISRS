#!/bin/bash
#SBATCH --job-name="10_output"
#SBATCH --time=5:00:00  # CHANGE as needed - walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # processor core(s) per node - just use 1 - it doesn't parallel
#SBATCH --mail-user="rsschwartz@uri.edu" #CHANGE to your user email address
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm_%x_%j.out
#SBATCH --error=slurm_%x_%j.err
#SBATCH -p uri-cpu

cd $SLURM_SUBMIT_DIR

module purge

#CHANGE IF NOT ON UNITY
module load uri/main bio/Biopython/1.79-intel-2022a

OUTFOLDER=../../SISRS_Small_test #CHANGE to analysis directory

../scripts/sisrs output -m 1 2 3 -d $OUTFOLDER #CHANGE the numbers as needed
