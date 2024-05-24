#!/bin/bash
#SBATCH --job-name="06_align"
#SBATCH --time=120:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=20   # CHANGE THIS to the number of processors
#SBATCH --mail-user="user@example.edu" #CHANGE THIS to your user email address
#SBATCH --mail-type=ALL
#SBATCH --exclusive

cd $SLURM_SUBMIT_DIR

module purge

#CHANGE IF NOT on a URI system
module load Biopython/1.83-foss-2023b
module load Bowtie2/2.4.5-GCC-11.3.0
module load SAMtools/1.16.1-GCC-11.3.0

PROCESSORS=36 #CHANGE to number of processors
OUTFOLDER=../../SISRS_Small_test #CHANGE to analysis directory

SPP=($(cat ${OUTFOLDER}/TaxonList.txt))

ARRLEN=${#SPP[@]}

for (( i = 0; i < $ARRLEN; i++ ))
do
    ../scripts/sisrs align -d $OUTFOLDER -p $PROCESSORS -f ${SPP[i]}
done

