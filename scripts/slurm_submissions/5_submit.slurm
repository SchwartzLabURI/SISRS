#!/bin/bash
#SBATCH --job-name="5_step"
#SBATCH --time=10:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=10   # processor core(s) per node
#SBATCH --mail-user="user@email.edu" #user email address
#SBATCH --mail-type=END,FAIL
#SBATCH --exclusive
cd $SLURM_SUBMIT_DIR

module purge

module load Python/3.7.4-GCCcore-8.3.0
module load SciPy-bundle/2019.10-foss-2019b-Python-3.7.4
module load Bowtie2/2.3.5.1-GCC-8.3.0
module load FastQC/0.11.8-Java-1.8
module load BBMap/38.81-foss-2019b-Java-1.8
module load Biopython/1.75-foss-2019b-Python-3.7.4
module load Ray/2.3.1-foss-2019b
module load SAMtools/1.10-GCC-8.3.0
module load BEDTools/2.29.2-GCC-8.3.0

#python3 sisrs_05_setup_sisrs.py -dir <path to the output directory>
