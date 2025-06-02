#!/bin/bash
#SBATCH --job-name="02_trim"
#SBATCH --time=3:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH -p uri-cpu
#SBATCH --ntasks-per-node=10   # CHANGE THIS to the number of processors on your node
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

P=$SLURM_JOB_CPUS_PER_NODE

D=../SISRS_Small #CHANGE THIS to the directory with input data OR the output directory for a prior run when combined with --link
DIR=../../SISRS_Small_test #CHANGE THIS to the analysis (output) directory

#Note - the only taxa used will be in TaxonList.txt so edit that for fewer taxa

../scripts/sisrs trimmer -p $P -d $D -dir $DIR
