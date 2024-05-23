#!/bin/bash
#SBATCH --job-name="svdq"
#SBATCH --time=150:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # CHANGE processor core(s) per node
#SBATCH --mail-user="user@example.edu" #CHANGE to your email
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_%x_%j.out
#SBATCH --error=slurm_%x_%j.err

cd $SLURM_SUBMIT_DIR
module load PAUP/4.0a168-centos64

alignmentFile="../../SISRS_Small_test/SISRS_Run/alignment_pi_m3_nogap.phylip-relaxed" #Change to the path to your phylip file - don't make the path too long
logFile="${alignmentFile}.svdq.log"
treeFile="${alignmentFile}.svdq.tre"

nBoot=1000
nCpus=36 #CHANGE to match ntasks line if needed

date

echo "
log file=${logFile};
tonexus format=RelPHYLIP fromfile=${alignmentFile} tofile=${alignmentFile}.svdq.nex interleaved=yes replace=yes;
exe ${alignmentFile}.svdq.nex;
svdq bootstrap nreps=${nBoot} nthreads=${nCpus};
saveTrees file=${treeFile};
quit;
" | paup

date
