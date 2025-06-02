#!/bin/bash
#SBATCH --job-name="scfRun"
#SBATCH --time=4:00:00  # CHANGE walltime limit (HH:MM:SS)
#SBATCH --nodes=1  # number of nodes
#SBATCH --ntasks-per-node=20   # CHANGE to the number of processor core(s) per node
#SBATCH --mail-user="rsschwartz@uri.edu" #CHANGE to your email address
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_%x_%j.out
#SBATCH --error=slurm_%x_%j.err
#SBATCH -p uri-cpu


cd $SLURM_SUBMIT_DIR

#CHANGE IF NOT ON UNITY
module load mpich/4.2.1
module load uri/main
module load iq-tree/2.3.1

TREE=RAxML_bestTree.alignment_pi_m1_nogap #CHANGE to your filename (include the path as needed)
ALIGNMENT=../../SISRS_Small_test/SISRS_Run/alignment_pi_m1_nogap.phylip-relaxed #CHANGE to the file used to generate the tree
OUTPUT=alignment_pi_m1_nogap #CHANGE to a name that will indicate the output information relevant to this tree
T=$SLURM_JOB_CPUS_PER_NODE

iqtree2-mpi -te $TREE -s $ALIGNMENT --scf 100 --prefix $OUTPUT


