#!/bin/bash
#SBATCH --job-name="2_step"
#SBATCH --time=100:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=20   # CHANGE THIS to the number of processors on your node
#SBATCH --mail-user="rsschwartz@uri.edu" #CHANGE THIS to your email address
#SBATCH --mail-type=END,FAIL
#SBATCH --exclusive

cd $SLURM_SUBMIT_DIR

module purge

#CHANGE THIS IF NOT ON A URI SYSTEM
module load Python/3.7.4-GCCcore-8.3.0
module load SciPy-bundle/2019.10-foss-2019b-Python-3.7.4
module load Bowtie2/2.3.5.1-GCC-8.3.0
module load FastQC/0.11.8-Java-1.8
module load BBMap/38.81-foss-2019b-Java-1.8
module load Biopython/1.75-foss-2019b-Python-3.7.4
module load Ray/2.3.1-foss-2019b
module load SAMtools/1.10-GCC-8.3.0
module load BEDTools/2.29.2-GCC-8.3.0

P=20 #CHANGE THIS to the number of processors 
D=../SISRS_Small #CHANGE THIS to the directory with input data
DIR=~/SISRS_Small_test #CHANGE THIS to the analysis (output) directory

#Note - the only taxa used will be in TaxonList.txt so edit that for fewer taxa

python3 sisrs_02_read_trimmer.py -p $P -d $D -dir $DIR
