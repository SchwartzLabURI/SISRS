#!/bin/bash
#SBATCH --job-name="2c_step_chloroplast2"
#SBATCH --time=24:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # processor core(s) per node
#SBATCH --mail-user="" #CHANGE TO user email address
#SBATCH --mail-type=END,FAIL
#SBATCH --array=[0-1]%2 #CHANGE like 6b - bracketed numbers indicate number of total jobs you need(taxa-1) - 0 based, inclusive; following number = simultaneous
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err

cd $SLURM_SUBMIT_DIR

module purge

#CHANGE THESE IF NOT ON A URI SYSTEM

module load SPAdes/3.15.2-GCC-10.3.0
module load SciPy-bundle/2021.10-foss-2021b
module load BLAST+/2.11.0-gompi-2021a 
module load Bowtie2/2.3.5.1-GCC-8.3.0
module load BBMap/38.96-foss-2020b
module load sympy/1.8-foss-2021a

#    fo = sys.argv[1]
#    taxon = sys.argv[2]
#    proc = sys.argv[3]

PROC=36 #CHANGE to number of processors (see ntasks above)
DIR=~/SISRS_plant_test/ #CHANGE THIS to the path for the output of your analysis

SPP=($(cat ${DIR}/TaxonList.txt))
echo ${SPP[$SLURM_ARRAY_TASK_ID]}


python3 sisrs_02b_filter_chloroplasts.py $DIR ${SPP[$SLURM_ARRAY_TASK_ID]} $PROC
