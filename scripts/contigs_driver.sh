#!/bin/bash
#SBATCH --job-name="consensus"
#SBATCH --time=2:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=20   # processor core(s) per node
#SBATCH --mail-user="yana_hrytsenko@uri.edu"
#SBATCH --mail-type=END,FAIL

cd $SLURM_SUBMIT_DIR

#module purge
#module load SAMtools/1.12-GCC-10.2.0
#module load BCFtools/1.12-GCC-10.2.0
#module load BBMap/38.87-foss-2020b


PTH=$1

#get all the directories in the path except Composite_Genome and SISRS_Run and contigs_outputs
FILELIST=( $( find $PTH -maxdepth 1 -type d -not -name Composite_Genome -not -name contigs_outputs -not -name SISRS_Run) )

#len of the list is 2 AotNan, CalJac
ARRLEN=${#FILELIST[@]}

for (( i = 0; i < $ARRLEN; i++ ))
do
  taxon_name="$(basename ${FILELIST[i]})" #get directory basename

  masked_ref="${FILELIST[i]}/masked_contigs.fa"

  indexed_masked_ref="${FILELIST[i]}/masked_contigs.fai"

  bam_file="${FILELIST[i]}/${taxon_name}.bam" #extract the taxon name

  vcf_file="${FILELIST[i]}/${taxon_name}.vcf"

  vcf_zipped="${FILELIST[i]}/${taxon_name}.vcf.gz"

  zipped_masked_fa="${FILELIST[i]}/masked_ref_seq_contigs.fa.gz"

  consensus_seq="${FILELIST[i]}/${taxon_name}_consensus.fa"

  consensus_seq_formatted="${FILELIST[i]}/${taxon_name}_single_line_format_consensus.fa"

  #step 1
  bcftools mpileup -Ou -f ${masked_ref} ${bam_file} | bcftools call -Ou -mv -o ${vcf_file}


  #step 2
  bgzip -c ${vcf_file} > ${vcf_zipped}


  #step 3
  bgzip -c ${masked_ref} > ${zipped_masked_fa}


  #step 4
  tabix -fp vcf ${vcf_zipped}


  #step 5
  samtools faidx ${masked_ref} -o ${indexed_masked_ref}


  #step 6
  zcat ${zipped_masked_fa} | bcftools consensus ${vcf_zipped} > ${consensus_seq}


  #step 7
  reformat.sh in=${consensus_seq} out=${consensus_seq_formatted} fastawrap=0

  #uncomment the line before to remove unnecessary files
  #rm ${zipped_masked_fa} ${masked_ref} ${indexed_masked_ref} ${vcf_zipped} ${vcf_file} ${consensus_seq}

done
