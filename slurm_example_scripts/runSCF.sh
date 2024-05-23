#!/bin/bash
#SBATCH --job-name="scfRun"
#SBATCH --time=48:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1  # number of nodes
#SBATCH --ntasks-per-node=36   # CHANGE to the number of processor core(s) per node
#SBATCH --mail-user="user@example.edu" #CHANGE to your email address
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_%x_%j.out
#SBATCH --error=slurm_%x_%j.err

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

cd $SLURM_SUBMIT_DIR

module load IQ-TREE/2.1.2-foss-2020a

TREE=RAxML_bestTree.alignment_pi_m3_nogap #CHANGE to your filename (include the path as needed)
ALIGNMENT=../../SISRS_Small_test/SISRS_Run/alignment_pi_m3_nogap.phylip-relaxed #CHANGE to the file used to generate the tree
OUTPUT=alignment_pi_m3_nogap #CHANGE to a name that will indicate the output information relevant to this tree
T=36 #CHANGE to the number of processors (20 or 36)

iqtree2 -t $TREE -s $ALIGNMENT --scf 100 --prefix $OUTPUT -nt $T
