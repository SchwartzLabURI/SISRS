#!/bin/bash
#SBATCH --job-name="6d_step"
#SBATCH --time=12:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # processor core(s) per node - only use 1 for this!
#SBATCH --mail-user="rsschwartz@uri.edu" #user email address
#SBATCH --mail-type=ALL
#SBATCH --array=[0-11]%12 #CHANGE as 6c - note each job only takes 1 core so you can run them all simultaneously
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
cd $SLURM_SUBMIT_DIR

#for advice on array jobs see https://github.com/nreid/using_array_jobs

module purge

#CHANGE if not on a URI system
module load Biopython/1.83-foss-2023b
module load SAMtools/1.16.1-GCC-11.3.0

OUTFOLDER=../../SISRS_Small_test #CHANGE to analysis directory
MINREAD=3 #CHANGE to set the threshold for the number of reads to call a site
THRESHOLD=1 #CHANGE to set the threshold [<=1] for the proportion of sites required to be the same to call a site

SPP=($(cat ${OUTFOLDER}/TaxonList.txt))
echo ${SPP[@]}
echo ${SPP[$SLURM_ARRAY_TASK_ID]}

./sisrs pileup2 -d $OUTFOLDER -m $MINREAD -t $THRESHOLD -s ${SPP[$SLURM_ARRAY_TASK_ID]}

