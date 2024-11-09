#!/bin/bash
#SBATCH --job-name="05_setup2"
#SBATCH --time=10:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=20   # CHANGE THIS to processor core(s) per node
#SBATCH --mail-user="rsschwartz@uri.edu" #CHANGE THIS to your user email address
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_%x_%j.out
#SBATCH --error=slurm_%x_%j.err
#SBATCH -p uri-cpu

cd $SLURM_SUBMIT_DIR

module purge

#CHANGE THIS IF NOT ON UNITY
module load uri/main bio/Biopython/1.79-intel-2022a
module load bbmap/39.01
module load bio/SAMtools/1.14-GCC-11.2.0
module load bio/Bowtie2/2.4.5-GCC-11.3.0

P=$SLURM_JOB_CPUS_PER_NODE

DIR=../../SISRS_Small_test #CHANGE THIS to the analysis directory

../scripts/sisrs setup2 -p $P -d $DIR
