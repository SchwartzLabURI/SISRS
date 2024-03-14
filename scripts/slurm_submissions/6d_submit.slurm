#!/bin/bash
#SBATCH --job-name="6d_step"
#SBATCH --time=120:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # processor core(s) per node - only use 1 for this!
#SBATCH --mail-user="user@example.edu" #user email address
#SBATCH --mail-type=ALL
#SBATCH -o slurm_%x_%A_%a.out
#SBATCH -e slurm_%x_%A_%a.err

cd $SLURM_SUBMIT_DIR

module purge

#CHANGE if not on a URI system
module load Biopython/1.83-foss-2023b
module load SAMtools/1.16.1-GCC-11.3.0

OUTFOLDER=~/SISRS_Small_test #CHANGE to analysis directory
MINREAD=3 #CHANGE to set the threshold for the number of reads to call a site
THRESHOLD=1 #CHANGE to set the threshold [<=1] for the proportion of sites required to be the same to call a site

SPP=($(cat ${OUTFOLDER}/TaxonList.txt))
ARRLEN=${#SPP[@]}

for (( i = 0; i < $ARRLEN; i++ ))
do
    python3 sisrs_06d_pileup2.py -d $OUTFOLDER -m $MINREAD -t $THRESHOLD -s ${SPP[i]}
done
