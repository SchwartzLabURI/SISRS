#!/bin/bash
#SBATCH --job-name="12_contig_align"
#SBATCH --time=4:00:00  # CHANGE as needed - walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=20   # CHANGE processor core(s) per node
#SBATCH --mail-user="rsschwartz@uri.edu" #CHANGE user email address
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_%x_%j.out
#SBATCH --error=slurm_%x_%j.err
#SBATCH -p uri-cpu
#SBATCH --constraint=avx512

cd $SLURM_SUBMIT_DIR

module purge

#CHANGE IF NOT ON UNITY
module load uri/main bio/Biopython/1.79-intel-2022a
module load MAFFT/7.505-GCC-11.3.0-with-extensions

P=$SLURM_JOB_CPUS_PER_NODE

OUTFOLDER=../../SISRS_Small_test/ #CHANGE THIS to the analysis directory

../scripts/sisrs contig_alignment -p $P -d $OUTFOLDER
