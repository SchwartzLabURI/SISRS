#!/bin/bash
#SBATCH --job-name="5_step"
#SBATCH --time=10:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # CHANGE THIS to processor core(s) per node
#SBATCH --mail-user="rsschwartz@uri.edu" #CHANGE THIS to your user email address
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_%x_%j.out
#SBATCH --error=slurm_%x_%j.err
#SBATCH --exclusive

cd $SLURM_SUBMIT_DIR

module purge

#CHANGE THIS IF NOT on a URI system
module load Biopython/1.83-foss-2023b
module load Bowtie2/2.4.5-GCC-11.3.0
module load BBMap/39.01-GCC-11.3.0
module load SAMtools/1.16.1-GCC-11.3.0

P=36 #CHANGE THIS to the number of processors
DIR=../../SISRS_Small_test #CHANGE THIS to the analysis directory
python3 sisrs_05_setup_sisrs.py -p $P -d $DIR
