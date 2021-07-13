#!/bin/bash
#SBATCH --job-name="RAxML"
#SBATCH --time=96:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=20   # processor core(s) per node 
#SBATCH --mail-user="user@example.edu"
#SBATCH --mail-type=END,FAIL

cd $SLURM_SUBMIT_DIR

module purge

module load RAxML/8.2.12-intel-2019b-hybrid-avx2

raxmlHPC -s <path to alignment file> -n <name of output file> -m ASC_GTRGAMMA --asc-corr="lewis" -T 20 -f a -p $RANDOM -N 100 -x $RANDOM
