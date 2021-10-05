#!/bin/bash
#SBATCH --job-name="scfRun"
#SBATCH --time=48:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1  # number of nodes
#SBATCH --ntasks-per-node=20   # CHANGE to the number of processor core(s) per node
#SBATCH --mail-user="jewelvoyer@uri.edu" #CHANGE to your email address
#SBATCH --mail-type=END,FAIL
#SBATCH --output="out_scfRun"
#SBATCH --error="err_scfRun"

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

cd $SLURM_SUBMIT_DIR

module load IQ-TREE/1.7-beta16-foss-2018b

TREE=RAxML_bipartitionsBranchLabels.M2_nogap_clean #CHANGE to your filename (include the path as needed)
ALIGNMENT=alignment_pi_m2_nogap_clean.phylip-relaxed #CHANGE to the file used to generate the tree
OUTPUT=ALGN1_SCF #CHANGE to a name that will indicate the output information relevant to this tree
T=20 #CHANGE to the number of processors (20 or 36)

iqtree -t $TREE -s $ALIGNMENT --scf 100 --prefix ALGN1_SCF -nt $T

