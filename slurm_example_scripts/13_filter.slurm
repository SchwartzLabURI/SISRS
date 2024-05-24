#!/bin/sh
#SBATCH --job-name="13_filter"
#SBATCH --time=45:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # CHANGE processor core(s) per node
#SBATCH --mail-user="rsschwartz@uri.edu" #CHANGE user email address
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_%x_%j.out
#SBATCH --error=slurm_%x_%j.err
#SBATCH --exclusive

cd $SLURM_SUBMIT_DIR

module purge

module load Biopython/1.78-foss-2020b
module load MAFFT/7.475-gompi-2020b-with-extensions


P=36 #CHANGE THIS to the number of processors

OUTFOLDER=../../SISRS_Small_test/ #CHANGE THIS to the analysis directory
ALIGNMENT_INFO=alignment_pi_locs_m2.txt #CHANGE
THRESHOLD=3 #CHANGE
LOC_LENGTH=160
MAX_DIST=0.08 #CHANGE - 0.08 = throw out loci more than 8% dif from composite
NUM_MISS=0 #CHANGE number missing in locus alignment

../scripts/sisrs filter -d $OUTFOLDER -a $ALIGNMENT_INFO -t $THRESHOLD -p $P -l $LOC_LENGTH -m $MAX_DIST -n $NUM_MISS
