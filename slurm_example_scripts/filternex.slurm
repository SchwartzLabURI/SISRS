#!/bin/bash
#SBATCH --job-name="filter_nexus"
#SBATCH --time=45:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # processor core(s) per node - just use 1 - it doesn't parallel
#SBATCH --mail-user="rsschwartz@uri.edu" #CHANGE to your user email address
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_%x_%j.out
#SBATCH --error=slurm_%x_%j.err

cd $SLURM_SUBMIT_DIR

module purge

#CHANGE if not on a URI system
module load Biopython/1.81-foss-2022b
module load MAFFT/7.475-gompi-2020b-with-extensions

OUTFOLDER=../../../Burmeistera_markers #CHANGE to analysis directory

../scripts/sisrs filter_nexus -m 4 -f ${OUTFOLDER}/SISRS_Run/alignment_pi.nex
