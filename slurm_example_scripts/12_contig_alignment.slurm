#!/bin/sh
#SBATCH --job-name="12_contig_align"
#SBATCH --time=4:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # CHANGE processor core(s) per node
#SBATCH --mail-user="rsschwartz@uri.edu" #CHANGE user email address
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_%x_%j.out
#SBATCH --error=slurm_%x_%j.err
#SBATCH --exclusive

cd $SLURM_SUBMIT_DIR

module purge

module load Biopython/1.83-foss-2023b
module load MAFFT/7.505-GCC-11.3.0-with-extensions

#python sisrs_07_b_contigs_alignment.py -p <number of processors to use> -d <path to the output directory>

P=36 #CHANGE THIS to the number of processors

OUTFOLDER=../../SISRS_Small_test/ #CHANGE THIS to the analysis directory

../scripts/sisrs contig_alignment -p $P -d $OUTFOLDER
