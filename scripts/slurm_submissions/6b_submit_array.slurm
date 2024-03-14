#!/bin/bash
#SBATCH --job-name="6b_step"
#SBATCH --time=12:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # processor core(s) per node - only use 1 for this!
#SBATCH --mail-user="rsschwartz@uri.edu" #CHANGE to user email address
#SBATCH --mail-type=ALL
#SBATCH --array=[0-11]%4 #bracketed numbers indicate number of total jobs you need(taxa-1) - 0 based, inclusive; following number = simultaneous
#SBATCH -o slurm_%x_%A_%a.out
#SBATCH -e slurm_%x_%A_%a.err
cd $SLURM_SUBMIT_DIR

#for advice on array jobs see https://github.com/nreid/using_array_jobs

module purge

#CHANGE IF not on a URI system
module load Biopython/1.83-foss-2023b
module load SAMtools/1.16.1-GCC-11.3.0

OUTFOLDER=../../SISRS_Small_test #CHANGE to analysis directory

SPP=($(cat ${OUTFOLDER}/TaxonList.txt))
echo ${SPP[@]}
echo ${SPP[$SLURM_ARRAY_TASK_ID]}

python3 sisrs_06b_pileup.py -d $OUTFOLDER -s ${SPP[$SLURM_ARRAY_TASK_ID]}
