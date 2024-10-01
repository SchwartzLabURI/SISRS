#!/bin/bash
#SBATCH --job-name="09_pileup"
#SBATCH --time=12:00:00  # CHANGE as needed - walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # processor core(s) per node - only use 1 for this!
#SBATCH --mail-user="user@example.edu" #user email address
#SBATCH --mail-type=ALL
#SBATCH -o slurm_%x_%A_%a.out
#SBATCH -e slurm_%x_%A_%a.err
#SBATCH -p uri-cpu

cd $SLURM_SUBMIT_DIR

module purge

#CHANGE IF NOT ON UNITY
module load uri/main bio/Biopython/1.79-intel-2022a
module load samtools/1.14

OUTFOLDER=../../SISRS_Small_test #CHANGE to analysis directory
MINREAD=3 #CHANGE to set the threshold for the number of reads to call a site
THRESHOLD=1 #CHANGE to set the threshold [<=1] for the proportion of sites required to be the same to call a site

SPP=($(cat ${OUTFOLDER}/TaxonList.txt))
ARRLEN=${#SPP[@]}

for (( i = 0; i < $ARRLEN; i++ ))
do
    ../scripts/sisrs pileup2 -d $OUTFOLDER -m $MINREAD -t $THRESHOLD -s ${SPP[i]}
done
