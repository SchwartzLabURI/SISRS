#!/bin/bash
#SBATCH --job-name="08_align"
#SBATCH --time=12:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=20   # CHANGE processor core(s) per node
#SBATCH --mail-user="rsschwartz@uri.edu" #CHANGE user email address
#SBATCH --mail-type=ALL
#SBATCH -p uri-cpu
#SBATCH --array=[0-11]%4 #CHANGE like 6b - bracketed numbers indicate number of total jobs you need(taxa-1) - 0 based, inclusive; following number = simultaneous
#SBATCH -o slurm_%x_%A_%a.out
#SBATCH -e slurm_%x_%A_%a.err

cd $SLURM_SUBMIT_DIR

#for advice on array jobs see https://github.com/nreid/using_array_jobs

module purge

#CHANGE THIS IF NOT ON UNITY
module load uri/main bio/Biopython/1.79-intel-2022a
module load bio/SAMtools/1.14-GCC-11.2.0
module load bio/Bowtie2/2.4.5-GCC-11.3.0

PROCESSORS=$SLURM_JOB_CPUS_PER_NODE

OUTFOLDER=../../SISRS_Small_test #CHANGE to the analysis directory

SPP=($(cat ${OUTFOLDER}/TaxonList.txt))
echo ${SPP[@]}
echo ${SPP[$SLURM_ARRAY_TASK_ID]}

../scripts/sisrs align2 -d $OUTFOLDER -p $PROCESSORS -f ${SPP[$SLURM_ARRAY_TASK_ID]}
