#!/bin/bash
#SBATCH --job-name="consensus"
#SBATCH --time=2:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=20   # processor core(s) per node
#SBATCH --mail-user="yana_hrytsenko@uri.edu"
#SBATCH --mail-type=END,FAIL

cd $SLURM_SUBMIT_DIR

#module purge
#module load BCFtools/1.12-GCC-10.2.0

PTH=$1

#get all the directories in the path except Composite_Genome and SISRS_Run and contigs_outputs
FILELIST=( $( find $PTH -maxdepth 1 -type d -not -name Composite_Genome -not -name contigs_outputs -not -name SISRS_Run) )

#len of the list is 2 AotNan, CalJac
ARRLEN=${#FILELIST[@]}

for (( i = 0; i < $ARRLEN; i++ ))
do
  taxon_name="$(basename ${FILELIST[i]})" #get directory basename

  bam_file="${FILELIST[i]}/${taxon_name}.bam" #input

  vcf_zipped="${FILELIST[i]}/${taxon_name}.vcf.gz" #output

  #do pileup using no reference and report allelic depth
  #then call variants using all sites (M) and zip the output
  bcftools mpileup -Ou --no-reference -a FORMAT/AD ${bam_file} | bcftools call -Oz -mM -o ${vcf_zipped}

done
