#!/bin/bash
#SBATCH --job-name="6b_step"
#SBATCH --time=12:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # processor core(s) per node - only use 1 for this!
#SBATCH --mail-user="user@example.edu" #CHANGE to user email address
#SBATCH --mail-type=ALL
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
cd $SLURM_SUBMIT_DIR


module purge

#CHANGE IF not on a URI system
module load Biopython/1.83-foss-2023b
module load SAMtools/1.16.1-GCC-11.3.0

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


OUTFOLDER=~/SISRS_Small_test #CHANGE to analysis directory

SPP=($(cat ${OUTFOLDER}/TaxonList.txt))

ARRLEN=${#SPP[@]}

for (( i = 0; i < $ARRLEN; i++ ))
do
    python3 sisrs_06b_pileup.py -d $OUTFOLDER -s ${SPP[i]}
done


#move all .out and .err files to their correct dirs so they dont clutter the slurm submission directory

mv 6*.out $DIRECTORY_OUT_6

mv 6*.err $DIRECTORY_ERR_6
