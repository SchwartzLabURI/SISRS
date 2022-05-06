#!/bin/sh
#SBATCH --job-name="7_b"
#SBATCH --time=45:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=20   # CHANGE processor core(s) per node
#SBATCH --mail-user="user@example.edu" #CHANGE user email address
#SBATCH --mail-type=END,FAIL
#SBATCH --exclusive

cd $SLURM_SUBMIT_DIR

module purge

module load MAFFT/7.475-gompi-2020b-with-extensions

#python sisrs_07_b_contigs_alignment.py -p <number of processors to use> -d <path to the output directory>

P=20 #CHANGE THIS to the number of processors

OUTFOLDER=~/SISRS_Small_test/ #CHANGE THIS to the analysis directory

python3 sisrs_07_b_contigs_alignment.py -p $P -d $OUTFOLDER
