#!/bin/bash
#SBATCH --job-name="2_step"
#SBATCH --time=100:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # CHANGE THIS to the number of processors on your node
#SBATCH --mail-user="rsschwartz@uri.edu" #CHANGE THIS to your email address
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --output=slurm_%x_%j.out
#SBATCH --error=slurm_%x_%j.err

cd $SLURM_SUBMIT_DIR

module purge

#CHANGE THIS IF NOT ON A URI SYSTEM
module load Biopython/1.83-foss-2023b
module load FastQC/0.11.9-Java-11
module load BBMap/39.01-GCC-11.3.0

P=36 #CHANGE THIS to the number of processors

D=../../SISRS_Small #CHANGE THIS to the directory with input data
DIR=../../SISRS_Small_test #CHANGE THIS to the analysis (output) directory

#Note - the only taxa used will be in TaxonList.txt so edit that for fewer taxa

python3 sisrs_02_read_trimmer.py -p $P -d $D -dir $DIR
