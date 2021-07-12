#!/bin/sh
#SBATCH --job-name="7_a"
#SBATCH --time=45:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=10   # processor core(s) per node
#SBATCH --mail-user="user@email.edu" #user email address
#SBATCH --mail-type=END,FAIL
#SBATCH --exclusive
cd $SLURM_SUBMIT_DIR

module purge

module load SAMtools/1.12-GCC-10.2.0
module load BCFtools/1.12-GCC-10.2.0
module load BBMap/38.87-foss-2020b
module load BEDTools/2.30.0-GCC-10.2.0
module load Biopython/1.78-foss-2020b



#python sisrs_07_a_contigs_processing.py -trh <number of threads to use> -dir <path to the output directory>
