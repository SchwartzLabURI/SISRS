#!/bin/bash
#SBATCH --job-name="07_pileup"
#SBATCH --time=12:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # processor core(s) per node - only use 1 for this!
#SBATCH --mail-user="rsschwartz@uri.edu" #CHANGE to user email address
#SBATCH --mail-type=ALL
#SBATCH --array=[0-11]%12 #bracketed numbers indicate number of total jobs you need(taxa-1) - 0 based, inclusive; following number = simultaneous
#SBATCH -o slurm_%x_%A_%a.out
#SBATCH -e slurm_%x_%A_%a.err
#SBATCH -p uri-cpu

cd $SLURM_SUBMIT_DIR

#for advice on array jobs see https://github.com/nreid/using_array_jobs

module purge

#CHANGE IF NOT ON UNITY
module load uri/main bio/Biopython/1.79-intel-2022a
module load samtools/1.14

OUTFOLDER=../../SISRS_Small_test #CHANGE to analysis directory

SPP=($(cat ${OUTFOLDER}/TaxonList.txt))
echo ${SPP[@]}
echo ${SPP[$SLURM_ARRAY_TASK_ID]}

../scripts/sisrs pileup -d $OUTFOLDER -s ${SPP[$SLURM_ARRAY_TASK_ID]}
