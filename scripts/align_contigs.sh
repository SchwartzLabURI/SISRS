#!/bin/bash
#SBATCH --job-name="7b"
#SBATCH --time=2:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=20   # processor core(s) per node
#SBATCH --mail-user="yana_hrytsenko@uri.edu"
#SBATCH --mail-type=END,FAIL

cd $SLURM_SUBMIT_DIR

#module load MAFFT/7.475-gompi-2020b-with-extensions

PTH=$1 #input path with all the contig files

PTH_OUT=$2 #output path to store alignments

THR=$3 #num threads to use #1 node, 20 processors per node -p 20;

#get all the files in the dir
FILELIST=( $( find $PTH -maxdepth 1 -type f ) )

#get nuber of contigs files to process
ARRLEN=${#FILELIST[@]}

for (( i = 0; i < $ARRLEN; i++ ))
do
  out_file="$(basename ${FILELIST[i]})" # contig file
  mafft --auto --thread $THR ${FILELIST[i]} > $PTH_OUT${out_file}
done
