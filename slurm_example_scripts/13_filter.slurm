#!/bin/bash
#SBATCH --job-name="13_filter"
#SBATCH --time=5:00:00  # CHANGE as needed - walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=20   # CHANGE processor core(s) per node
#SBATCH --mail-user="rsschwartz@uri.edu" #CHANGE user email address
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_%x_%j.out
#SBATCH --error=slurm_%x_%j.err
#SBATCH -p uri-cpu

cd $SLURM_SUBMIT_DIR

module purge

#CHANGE IF NOT ON UNITY
module load uri/main 
module load conda/latest
if ! conda env list | grep -q sisrs; then
  conda config --add channels bioconda
  conda create -n sisrs pyvcf3 Biopython numpy pandas psutil -y
  pip install cyvcf2
fi
conda activate sisrs

module load MAFFT/7.505-GCC-11.3.0-with-extensions

P=$SLURM_JOB_CPUS_PER_NODE

OUTFOLDER=../../SISRS_Small_test2/ #CHANGE THIS to the analysis directory
ALIGNMENT_INFO=alignment_pi_locs_m2.txt #CHANGE
THRESHOLD=3 #CHANGE
LOC_LENGTH=160
MAX_DIST=0.08 #CHANGE - 0.08 = throw out loci more than 8% dif from composite
NUM_MISS=0 #CHANGE number missing in locus alignment

../scripts/sisrs filter -d $OUTFOLDER -a $ALIGNMENT_INFO -t $THRESHOLD -p $P -l $LOC_LENGTH -m $MAX_DIST -n $NUM_MISS
