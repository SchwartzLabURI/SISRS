#!/bin/bash
#SBATCH --job-name="6d_step"
#SBATCH --time=120:00:00  # walltime limit (HH:MM:SS)
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
module load Python/3.7.4-GCCcore-8.3.0
module load SciPy-bundle/2019.10-foss-2019b-Python-3.7.4
module load Bowtie2/2.3.5.1-GCC-8.3.0
module load FastQC/0.11.8-Java-1.8
module load BBMap/38.81-foss-2019b-Java-1.8
module load Biopython/1.75-foss-2019b-Python-3.7.4
module load Ray/2.3.1-foss-2019b
module load SAMtools/1.10-GCC-8.3.0
module load BEDTools/2.29.2-GCC-8.3.0

OUTFOLDER=~/SISRS_Small_test #CHANGE to analysis directory
MINREAD=3 #CHANGE to set the threshold for the number of reads to call a site
THRESHOLD=1 #CHANGE to set the threshold [<=1] for the proportion of sites required to be the same to call a site

SPP=($(cat ${OUTFOLDER}/TaxonList.txt))
echo ${SPP[@]}
echo ${SPP[$SLURM_ARRAY_TASK_ID]}

python3 sisrs_06d_pileup2.py -d $OUTFOLDER -m $MINREAD -t $THRESHOLD -s ${SPP[$SLURM_ARRAY_TASK_ID]}
