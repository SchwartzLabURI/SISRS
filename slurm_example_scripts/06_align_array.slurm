#!/bin/bash
#SBATCH --job-name="06_align"
#SBATCH --time=12:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # CHANGE THIS to the number of processors
#SBATCH --mail-user="rsschwartz@uri.edu" #CHANGE THIS to your user email address
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --array=[0-11]%4 # CHANGE this second bracketed number to the total jobs you need(taxa-1); CHANGE after % to number of  simultaneous jobs
#SBATCH -o slurm_%x_%A_%a.out
#SBATCH -e slurm_%x_%A_%a.err

cd $SLURM_SUBMIT_DIR

#for advice on array jobs see https://github.com/nreid/using_array_jobs

module purge

#CHANGE IF NOT on a URI system
module load Biopython/1.83-foss-2023b
module load Bowtie2/2.4.5-GCC-11.3.0
module load SAMtools/1.16.1-GCC-11.3.0

# create an array variable containing the folders names
PROCESSORS=36 #CHANGE to number of processors
OUTFOLDER=../../SISRS_Small_test #CHANGE to analysis directory

SPP=($(cat ${OUTFOLDER}/TaxonList.txt))
echo ${SPP[@]}
echo ${SPP[$SLURM_ARRAY_TASK_ID]}

../scripts/sisrs align -d $OUTFOLDER -p $PROCESSORS -f ${SPP[$SLURM_ARRAY_TASK_ID]}
