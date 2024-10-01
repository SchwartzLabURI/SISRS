#!/bin/bash
#SBATCH --job-name="svdq"
#SBATCH --time=10:00:00  # CHANGE - walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=20   # CHANGE processor core(s) per node
#SBATCH --mail-user="rsschwartz@uri.edu" #CHANGE to your email
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_%x_%j.out
#SBATCH --error=slurm_%x_%j.err
#SBATCH -p uri-cpu

cd $SLURM_SUBMIT_DIR

#CHANGE IF NOT ON UNITY
module load uri/main PAUP/4.0a168-centos64


alignmentFile="../../SISRS_Small_test/SISRS_Run/alignment_pi_m3_nogap.phylip-relaxed" #Change to the path to your phylip file - don't make the path too long
logFile="${alignmentFile}.svdq.log"
treeFile="${alignmentFile}.svdq.tre"

nBoot=1000
nCpus=$SLURM_JOB_CPUS_PER_NODE

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
