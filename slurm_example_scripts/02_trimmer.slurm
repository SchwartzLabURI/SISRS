#!/bin/bash
#SBATCH --job-name="02_trim"
#SBATCH --time=3:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH -p uri-cpu
#SBATCH --ntasks-per-node=10   # CHANGE THIS to the number of processors on your node
#SBATCH --mail-user="rsschwartz@uri.edu" #CHANGE THIS to your email address
#SBATCH --mail-type=ALL
#SBATCH --array=[0-11]%4 # CHANGE this second bracketed number to the total jobs you need(taxa-1); CHANGE after % to number of  simultaneous jobs
#SBATCH -o slurm_%x_%A_%a.out
#SBATCH -e slurm_%x_%A_%a.err

cd $SLURM_SUBMIT_DIR

module purge

#CHANGE THIS IF NOT ON UNITY
module load uri/main 
#module load bio/Biopython/1.83-foss-2023a
module load FastQC/0.11.9-Java-11
module load bbmap/39.01
module load conda/latest
conda activate sisrs

P=$SLURM_JOB_CPUS_PER_NODE

D=~/SISRS_Small #CHANGE THIS to the directory with input data OR the output directory for a prior run when combined with --link
OUTFOLDER=~/SISRS_Small_test2 #CHANGE THIS to the analysis (output) directory

# create an array variable containing the folders names
SPP=($(cat ${OUTFOLDER}/TaxonList.txt))
echo ${SPP[@]}
echo ${SPP[$SLURM_ARRAY_TASK_ID]}

#Note - the only taxa used will be in TaxonList.txt so edit that for fewer taxa

../scripts/sisrs trimmer -p $P -d $D -dir $OUTFOLDER -sp ${SPP[$SLURM_ARRAY_TASK_ID]}
