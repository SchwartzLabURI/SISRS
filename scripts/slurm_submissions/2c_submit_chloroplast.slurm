#!/bin/bash
#SBATCH --job-name="2c_step_chloroplast2"
#SBATCH --time=24:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # processor core(s) per node
#SBATCH --mail-user="rsschwartz@uri.edu" #CHANGE TO user email address
#SBATCH --mail-type=END,FAIL
#SBATCH --array=[0-11]%4 #CHANGE like 6b - bracketed numbers indicate number of total jobs you need(taxa-1) - 0 based, inclusive; following number = simultaneous
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err

cd $SLURM_SUBMIT_DIR

module purge

#CHANGE THESE IF NOT ON A URI SYSTEM

module load Python/3.7.4-GCCcore-8.3.0
module load SciPy-bundle/2019.10-foss-2019b-Python-3.7.4
module load BBMap/38.81-foss-2019b-Java-1.8
module load Biopython/1.75-foss-2019b-Python-3.7.4
module load SAMtools/1.10-GCC-8.3.0
module load sympy/1.7.1-foss-2020b sympy/1.8-foss-2021a
module load scikit-learn/0.23.2-fosscuda-2020b
module load SPAdes/3.14.0-GCC-8.3.0-Python-3.7.4
module load Bowtie2/2.3.5.1-GCC-8.3.0

module load BLAST+/2.10.1-gompi-2020a

#    fo = sys.argv[1]
#    taxon = sys.argv[2]
#    proc = sys.argv[3]

PROC=36 #CHANGE to number of processors (see ntasks above)
DIR=~/SISRS_Small_test #CHANGE THIS to the path for the output of your analysis

SPP=($(cat ${DIR}/TaxonList.txt))
echo ${SPP[$SLURM_ARRAY_TASK_ID]}


python3 sisrs_02b_filter_chloroplasts.py $DIR ${SPP[$SLURM_ARRAY_TASK_ID]} $PROC 
