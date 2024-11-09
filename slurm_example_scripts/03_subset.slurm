#!/bin/bash
#SBATCH --job-name="03_subset"
#SBATCH --time=1:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=10   # CHANGE THIS to the number of processors
#SBATCH -p uri-cpu
#SBATCH --mail-user="rsschwartz@uri.edu" #CHANGE THIS to your email address
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_%x_%j.out
#SBATCH --error=slurm_%x_%j.err

cd $SLURM_SUBMIT_DIR

module purge

#CHANGE THIS IF NOT ON UNITY
module load uri/main bio/Biopython/1.79-intel-2022a
module load FastQC/0.11.9-Java-11
module load bbmap/39.01

GENOME=1000000000 #CHANGE THIS to the approximate size of your genome
DIR=../../SISRS_Small_test #CHANGE this to your analysis folder

../scripts/sisrs subset -gs $GENOME -d $DIR
