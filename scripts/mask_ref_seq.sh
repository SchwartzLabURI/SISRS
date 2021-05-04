#!/bin/bash
#SBATCH --job-name="masking"
#SBATCH --time=2:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=20   # processor core(s) per node
#SBATCH --mail-user="yana_hrytsenko@uri.edu"
#SBATCH --mail-type=END,FAIL

cd $SLURM_SUBMIT_DIR

#module purge
#module load BEDTools/2.30.0-GCC-10.2.0

PTH=$1

#get all the directories in the path except Composite_Genome, SISRS_Run, contigs_outputs
FILELIST=( $( find $PTH -maxdepth 1 -type d -not -name Composite_Genome -not -name contigs_outputs -not -name SISRS_Run) )


#len of the list is 2 AotNan, CalJac
ARRLEN=${#FILELIST[@]}

for (( i = 0; i < $ARRLEN; i++ ))
do
  taxon_name="$(basename ${FILELIST[i]})" #get directory basename


  bam_file="${FILELIST[i]}/${taxon_name}.bam" #extract the taxon name


  ref_fa="${FILELIST[i]}/contigs.fa"


  masked_ref="${FILELIST[i]}/masked_contigs.fa"


  bed_file="${FILELIST[i]}/to_mask.bed"


  #step 0 Mask the reference sequence
  bedtools genomecov -ibam ${bam_file} -bga | awk '{if ($4 == 0) print $0}' > ${bed_file}
  bedtools maskfasta  -fi ${ref_fa} -bed ${bed_file} -fo ${masked_ref}
done
