#!/bin/sh
#SBATCH --job-name="7_a"
#SBATCH --time=4:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=10   # processor core(s) per node
#SBATCH --mail-user="rsschwartz@uri.edu" #user email address
#SBATCH --mail-type=END,FAIL
#SBATCH --exclusive
cd $SLURM_SUBMIT_DIR

module purge

module load BCFtools/1.15.1-GCC-11.3.0
module load Biopython/1.83-foss-2023b

#Optionally specify coverage threshold: '-c', '--cov' (default set to 3) and heterozygosity threshold: '-z', '--hz'(default set to 0.01)

T=1 #CHANGE THIS to the threshold number

OUTFOLDER=../../SISRS_Small_test #CHANGE to analysis directory

python3 sisrs_07_a_contigs_processing.py -t $T -d $OUTFOLDER
