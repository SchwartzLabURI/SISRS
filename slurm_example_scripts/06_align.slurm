#!/bin/bash
#SBATCH --job-name="06_align"
#SBATCH --time=12:00:00  # CHANGE as needed - walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=20   # CHANGE THIS to the number of processors
#SBATCH --mail-user="user@example.edu" #CHANGE THIS to your user email address
#SBATCH --mail-type=ALL
#SBATCH -p uri-cpu

cd $SLURM_SUBMIT_DIR

module purge

#CHANGE THIS IF NOT ON UNITY
module load uri/main bio/Biopython/1.79-intel-2022a
module load bio/SAMtools/1.14-GCC-11.2.0
module load bio/Bowtie2/2.4.5-GCC-11.3.0

PROCESSORS=$SLURM_JOB_CPUS_PER_NODE

OUTFOLDER=../../SISRS_Small_test #CHANGE to analysis directory

SPP=($(cat ${OUTFOLDER}/TaxonList.txt))

ARRLEN=${#SPP[@]}

for (( i = 0; i < $ARRLEN; i++ ))
do
    ../scripts/sisrs align -d $OUTFOLDER -p $PROCESSORS -f ${SPP[i]}
done

