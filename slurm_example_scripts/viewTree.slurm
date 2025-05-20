#!/bin/bash
#SBATCH --job-name="tree"
#SBATCH --time=1:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1  # number of nodes
#SBATCH --ntasks-per-node=2 # R likes 2?
#SBATCH --mail-user="rsschwartz@uri.edu" #CHANGE to your email address
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_%x_%j.out
#SBATCH --error=slurm_%x_%j.err
#SBATCH -p uri-cpu

cd $SLURM_SUBMIT_DIR

#CHANGE IF NOT ON UNITY
module load uri/main R/4.3.2-gfbf-2023a

#local space for R packages and it won't ask about install location
mkdir -p ~/R-packages
export R_LIBS=~/R-packages

TREE=alignment_pi_m1_nogap.cf.tree #CHANGE to your filename
OUTGROUP=AotNan #CHANGE to your outgroup

Rscript viewTree.R $TREE $OUTGROUP
