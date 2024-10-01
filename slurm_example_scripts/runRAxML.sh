#!/bin/bash
#SBATCH --job-name="RAxML"
#SBATCH --time=6:00:00  # CHANGE - walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=20   # CHANGE processor core(s) per node 
#SBATCH --mail-user="rsschwartz@uri.edu" #CHANGE to your email
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_%x_%j.out
#SBATCH --error=slurm_%x_%j.err
#SBATCH -p uri-cpu

cd $SLURM_SUBMIT_DIR

module purge

#CHANGE IF NOT ON UNITY
module load intel-oneapi-mpi/2021.6.0 apptainer/latest shpc/0.1.26 raxml/8.2.12

ALIGNMENT=../../SISRS_Small_test/SISRS_Run/alignment_pi_m1_nogap.phylip-relaxed  #CHANGE - Put the path to your alignment file here

OUTPUT=$(basename ${ALIGNMENT%.*})  #output filename relevant to the alignment
P=$SLURM_JOB_CPUS_PER_NODE

raxmlHPC -s $ALIGNMENT -n $OUTPUT -m ASC_GTRGAMMA --asc-corr="lewis" -T $P -f a -p $RANDOM -N 100 -x $RANDOM
