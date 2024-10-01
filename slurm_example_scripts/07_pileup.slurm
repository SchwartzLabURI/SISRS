#!/bin/bash
#SBATCH --job-name="07_pileup"
#SBATCH --time=12:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # processor core(s) per node - only use 1 for this!
#SBATCH --mail-user="user@example.edu" #CHANGE to user email address
#SBATCH --mail-type=ALL
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH -p uri-cpu

cd $SLURM_SUBMIT_DIR

module purge

#CHANGE IF NOT ON UNITY
module load uri/main bio/Biopython/1.79-intel-2022a
module load samtools/1.14

OUTFOLDER=../../SISRS_Small_test #CHANGE to analysis directory

SPP=($(cat ${OUTFOLDER}/TaxonList.txt))

ARRLEN=${#SPP[@]}

for (( i = 0; i < $ARRLEN; i++ ))
do
    ../scripts/sisrs pileup -d $OUTFOLDER -s ${SPP[i]}
done
