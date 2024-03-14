#!/bin/bash
#SBATCH --job-name="3_step"
#SBATCH --time=1:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # CHANGE THIS to the number of processors
#SBATCH --mail-user="rsschwartz@uri.edu" #CHANGE THIS to your email address
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --output=slurm_%x_%j.out
#SBATCH --error=slurm_%x_%j.err

cd $SLURM_SUBMIT_DIR

module purge

#CHANGE THIS if not on a URI system
module load Biopython/1.83-foss-2023b
module load BBMap/39.01-GCC-11.3.0

GENOME=1000000000 #CHANGE THIS to the approximate size of your genome
DIR=../../SISRS_Small_test #CHANGE this to your analysis folder

python3 sisrs_03_read_subsetter.py -gs $GENOME -d $DIR
