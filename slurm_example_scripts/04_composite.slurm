#!/bin/bash
#SBATCH --job-name="04_composite"
#SBATCH --time=10:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=20   # CHANGE THIS to the number of processor core(s) per node
#SBATCH -p uri-cpu
#SBATCH --mail-user="rsschwartz@uri.edu" #CHANGE THIS to your email address
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_%x_%j.out
#SBATCH --error=slurm_%x_%j.err
#SBATCH --mem=50GB #CHANGE THIS for larger genomes

cd $SLURM_SUBMIT_DIR

module purge

#CHANGE THIS IF NOT ON UNITY
module load uri/main 
#module load bio/Biopython/1.79-intel-2022a
module load Ray/2.3.1-iimpi-2022a
module load conda/latest
conda activate sisrs

P=$SLURM_JOB_CPUS_PER_NODE

DIR=../../SISRS_Small_test2 #CHANGE THIS to the analysis directory

../scripts/sisrs composite -p $P -d $DIR
