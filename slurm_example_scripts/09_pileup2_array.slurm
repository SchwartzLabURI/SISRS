#!/bin/bash
#SBATCH --job-name="09_pileup"
#SBATCH --time=12:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # processor core(s) per node - only use 1 for this!
#SBATCH --mail-user="rsschwartz@uri.edu" #user email address
#SBATCH --mail-type=ALL
#SBATCH --array=[0-11]%12 #CHANGE to match the number of taxa - note each job only takes 1 core so you can run them all simultaneously
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH -p uri-cpu

cd $SLURM_SUBMIT_DIR

#for advice on array jobs see https://github.com/nreid/using_array_jobs

module purge

#CHANGE IF NOT ON UNITY
module load uri/main bio/Biopython/1.79-intel-2022a
module load samtools/1.14

OUTFOLDER=../../SISRS_Small_test #CHANGE to analysis directory
MINREAD=3 #CHANGE to set the threshold for the number of reads to call a site
THRESHOLD=1 #CHANGE to set the threshold [<=1] for the proportion of sites required to be the same to call a site

SPP=($(cat ${OUTFOLDER}/TaxonList.txt))
echo ${SPP[@]}
echo ${SPP[$SLURM_ARRAY_TASK_ID]}

../scripts/sisrs pileup2 -d $OUTFOLDER -m $MINREAD -t $THRESHOLD -s ${SPP[$SLURM_ARRAY_TASK_ID]}

