#!/bin/bash
#SBATCH --job-name="08_align"
#SBATCH --time=12:00:00  # CHANGE as needed - walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=20   # CHANGE processor core(s) per node
#SBATCH --mail-user="rsschwartz@uri.edu" #CHANGE user email address
#SBATCH --mail-type=ALL
#SBATCH -p uri-cpu
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err

cd $SLURM_SUBMIT_DIR

module purge

#CHANGE THIS IF NOT ON UNITY
module load uri/main #bio/Biopython/1.79-intel-2022a
module load bio/SAMtools/1.14-GCC-11.2.0
module load bio/Bowtie2/2.4.5-GCC-11.3.0
module load conda/latest
conda activate sisrs

PROCESSORS=$SLURM_JOB_CPUS_PER_NODE

OUTFOLDER=../../SISRS_Small_test #CHANGE to the analysis directory

SPP=($(cat ${OUTFOLDER}/TaxonList.txt))

ARRLEN=${#SPP[@]}

for (( i = 0; i < $ARRLEN; i++ ))
do
    ../scripts/sisrs align2 -d $OUTFOLDER -p $PROCESSORS -f ${SPP[i]}
done
