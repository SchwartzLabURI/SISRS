#!/bin/bash
#SBATCH --job-name="10_output"
#SBATCH --time=5:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # processor core(s) per node - just use 1 - it doesn't parallel
#SBATCH --mail-user="rsschwartz@uri.edu" #CHANGE to your user email address
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm_%x_%j.out
#SBATCH --error=slurm_%x_%j.err

cd $SLURM_SUBMIT_DIR

module purge

#CHANGE if not on a URI system
module load Biopython/1.83-foss-2023b

#CHANGE following example
#python3 sisrs_07_output_sisrs.py -m <list of number of species allowed to be missing> -d <main output folder>
OUTFOLDER=../../SISRS_Small_test #CHANGE to analysis directory

../scripts/sisrs output -m 1 2 3 -d $OUTFOLDER
