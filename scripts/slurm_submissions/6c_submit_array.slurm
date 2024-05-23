#!/bin/bash
#SBATCH --job-name="6c_step"
#SBATCH --time=12:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # CHANGE processor core(s) per node
#SBATCH --mail-user="rsschwartz@uri.edu" #CHANGE user email address
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --array=[0-11]%4 #CHANGE like 6b - bracketed numbers indicate number of total jobs you need(taxa-1) - 0 based, inclusive; following number = simultaneous
#SBATCH -o slurm_%x_%A_%a.out
#SBATCH -e slurm_%x_%A_%a.err

cd $SLURM_SUBMIT_DIR

#for advice on array jobs see https://github.com/nreid/using_array_jobs

module purge

#CHANGE if not on a URI system
module load Biopython/1.83-foss-2023b
module load Bowtie2/2.4.5-GCC-11.3.0
module load SAMtools/1.16.1-GCC-11.3.0

PROCESSORS=36 #CHANGE to the number of processors
OUTFOLDER=../../SISRS_Small_test #CHANGE to the analysis directory

SPP=($(cat ${OUTFOLDER}/TaxonList.txt))
echo ${SPP[@]}
echo ${SPP[$SLURM_ARRAY_TASK_ID]}

./sisrs align2 -d $OUTFOLDER -p $PROCESSORS -f ${SPP[$SLURM_ARRAY_TASK_ID]}
