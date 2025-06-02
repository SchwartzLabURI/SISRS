#!/bin/sh
#SBATCH --job-name="optimize"
#SBATCH --time=100:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=10   # processor core(s) per node
#SBATCH --mail-user="rsschwartz@uri.edu" #CHANGE user email address
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_%x_%j.out
#SBATCH --error=slurm_%x_%j.err

cd $SLURM_SUBMIT_DIR

module purge

module load uri/main bio/Biopython/1.79-intel-2022a
module load MAFFT/7.505-GCC-11.3.0-with-extensions
module load BLAST+/2.15.0-gompi-2023a

P=10 #CHANGE THIS to the number of processors
OUTFOLDER=../../../Centropogon_markers/ #CHANGE THIS to the analysis directory
MAX_CONTIGS=ALL #CHANGE to reduce output
LOC_LENGTH=1000
MAX_DIST=0.08 #0.08 = throw out loci more than 8% dif from composite

echo NUM_MISS ALIGNMENT_INFO THRESHOLD SISRS_sites num_contigs num_contigs2 contig_length num_pi_sites > ${OUTFOLDER}/optimize_markers2.tsv #collect results

../scripts/sisrs filter_nexus -m 20 24 26 30 -f ${OUTFOLDER}/SISRS_Run/alignment_pi.nex

for THRESHOLD in {1..3}; do #min num variable sites in contig
  for N in {0..30..4}; do #num spp missing in site alignment
    for NUM_MISS in {0..30..4}; do  # number spp missing in locus alignment
      ALIGNMENT_INFO=alignment_pi_locs_m${N}.txt

      ../scripts/sisrs filter -d $OUTFOLDER -a $ALIGNMENT_INFO -t $THRESHOLD -p $P -l $LOC_LENGTH -m $MAX_DIST -n "$NUM_MISS"

      ../scripts/sisrs filter2 -d $OUTFOLDER -m $MAX_CONTIGS

      SISRS_sites=$(head -1 ${OUTFOLDER}/SISRS_Run/alignment_pi_m${N}.phylip-relaxed | cut -d ' ' -f 3)

      num_contigs=$(grep -c '^>' ${OUTFOLDER}/SISRS_Run/contigs_for_probes.fa)

      num_contigs2=$(grep -c '^>' ${OUTFOLDER}/SISRS_Run/contigs_for_probes.fa.fasta)

      contig_length=$(grep -v '^>' ${OUTFOLDER}/SISRS_Run/contigs_for_probes.fa.fasta | wc -c)

      num_snps=$(python3 ../scripts/count_snps.py -f ${OUTFOLDER}) # number of variable pi sites across all locus alignments

      # information about splits for each site

      echo $NUM_MISS $ALIGNMENT_INFO ${THRESHOLD} ${SISRS_sites} $num_contigs $num_contigs2 $contig_length $num_snps >> ${OUTFOLDER}/optimize_markers2.tsv #collect results

      mv ${OUTFOLDER}/SISRS_Run/aligned_contigs2 ${OUTFOLDER}/SISRS_Run/aligned_contigs2_${NUM_MISS}_${ALIGNMENT_INFO}_${THRESHOLD}
      mv ${OUTFOLDER}/SISRS_Run/contigs_for_probes.fa.fasta ${OUTFOLDER}/SISRS_Run/contigs_for_probes_${NUM_MISS}_${ALIGNMENT_INFO}_${THRESHOLD}.fa.fasta
      mv ${OUTFOLDER}/SISRS_Run/contigs_for_probes.fa ${OUTFOLDER}/SISRS_Run/contigs_for_probes_${NUM_MISS}_${ALIGNMENT_INFO}_${THRESHOLD}.fa
      mv ${OUTFOLDER}/SISRS_Run/variation_in_contigs.txt ${OUTFOLDER}/SISRS_Run/variation_in_contigs_${NUM_MISS}_${ALIGNMENT_INFO}_${THRESHOLD}.txt

    done
  done
done



