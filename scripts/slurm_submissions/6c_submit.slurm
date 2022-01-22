#!/bin/bash
#SBATCH --job-name="6c_step"
#SBATCH --time=120:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # CHANGE processor core(s) per node
#SBATCH --mail-user="yana_hrytsenko@uri.edu" #CHANGE user email address
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err

cd $SLURM_SUBMIT_DIR


module purge

#CHANGE if not on a URI system
module load Python/3.7.4-GCCcore-8.3.0
module load SciPy-bundle/2019.10-foss-2019b-Python-3.7.4
module load Bowtie2/2.3.5.1-GCC-8.3.0
module load FastQC/0.11.8-Java-1.8
module load BBMap/38.81-foss-2019b-Java-1.8
module load Biopython/1.75-foss-2019b-Python-3.7.4
module load Ray/2.3.1-foss-2019b
module load SAMtools/1.10-GCC-8.3.0
module load BEDTools/2.29.2-GCC-8.3.0

#housekeeping
#if no 6_out_files dir mkdir 6_out_files
#if no 6_err_files dir mkdir 6_err_files

DIRECTORY_OUT_6="6_out_files"

if [ ! -d "$DIRECTORY_OUT_6" ]; then
  mkdir $DIRECTORY_OUT_6
fi

DIRECTORY_ERR_6="6_err_files"

if [ ! -d "$DIRECTORY_ERR_6" ]; then
  mkdir $DIRECTORY_ERR_6
fi

PROCESSORS=36 #CHANGE to the number of processors
OUTFOLDER=~/SISRS_Small_test #CHANGE to the analysis directory

SPP=($(cat ${OUTFOLDER}/TaxonList.txt))

ARRLEN=${#SPP[@]}

for (( i = 0; i < $ARRLEN; i++ ))
do
    python3 sisrs_06c_align2.py -d $OUTFOLDER -p $PROCESSORS -f ${SPP[i]}
done

#move all .out and .err files to their correct dirs so they dont clutter the slurm submission directory

mv 6*.out $DIRECTORY_OUT_6

mv 6*.err $DIRECTORY_ERR_6
