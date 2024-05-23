#!/bin/bash
#SBATCH --job-name="4_step"
#SBATCH --time=100:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # CHANGE THIS to the number of processor core(s) per node
#SBATCH --mail-user="rsschwartz@uri.edu" #CHANGE THIS to your email address
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_%x_%j.out
#SBATCH --error=slurm_%x_%j.err
#SBATCH --exclusive
#SBATCH --mem=500GB

cd $SLURM_SUBMIT_DIR

module purge

#CHANGE THIS if not on a URI system
module load Biopython/1.83-foss-2023b
module load Ray/2.3.1-foss-2019b

P=36 #CHANGE THIS to the number of processors
DIR=../../SISRS_Small_test #CHANGE THIS to the analysis directory
./sisrs composite -p $P -d $DIR
