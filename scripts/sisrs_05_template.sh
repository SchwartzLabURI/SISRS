#!/bin/sh
bowtie2 -p PROCESSORS -N 1 --local -x BOWTIE2-INDEX -U READS | samtools view -Su -@ PROCESSORS -F 4 - | samtools sort -@ PROCESSORS - -o SISRS_DIR/TAXA/TAXA_Temp.bam

samtools view -@ PROCESSORS -H SISRS_DIR/TAXA/TAXA_Temp.bam > SISRS_DIR/TAXA/TAXA_Header.sam

samtools view -@ PROCESSORS SISRS_DIR/TAXA/TAXA_Temp.bam | grep -v "XS:" | cat SISRS_DIR/TAXA/TAXA_Header.sam - | samtools view -@ PROCESSORS -b - > SISRS_DIR/TAXA/TAXA.bam

rm SISRS_DIR/TAXA/TAXA_Temp.bam
rm SISRS_DIR/TAXA/TAXA_Header.sam

samtools mpileup -f COMPOSITE_GENOME SISRS_DIR/TAXA/TAXA.bam > SISRS_DIR/TAXA/TAXA.pileups

python3 SCRIPT_DIR/specific_genome.py SISRS_DIR/TAXA COMPOSITE_GENOME

samtools faidx SISRS_DIR/TAXA/contigs.fa
bowtie2-build SISRS_DIR/TAXA/contigs.fa SISRS_DIR/TAXA/contigs -p PROCESSORS

bowtie2 -p PROCESSORS -N 1 --local -x SISRS_DIR/TAXA/contigs -U READS | samtools view -Su -@ PROCESSORS -F 4 - | samtools sort -@ PROCESSORS - -o SISRS_DIR/TAXA/TAXA_Temp.bam

samtools view -@ PROCESSORS -H SISRS_DIR/TAXA/TAXA_Temp.bam > SISRS_DIR/TAXA/TAXA_Header.sam
samtools view -@ PROCESSORS SISRS_DIR/TAXA/TAXA_Temp.bam | grep -v "XS:" | cat SISRS_DIR/TAXA/TAXA_Header.sam - | samtools view -@ PROCESSORS -b - > SISRS_DIR/TAXA/TAXA.bam

rm SISRS_DIR/TAXA/TAXA_Temp.bam
rm SISRS_DIR/TAXA/TAXA_Header.sam

samtools index SISRS_DIR/TAXA/TAXA.bam

samtools mpileup -f COMPOSITE_GENOME SISRS_DIR/TAXA/TAXA.bam > SISRS_DIR/TAXA/TAXA.pileups

python3 SCRIPT_DIR/get_pruned_dict.py SISRS_DIR/TAXA COMPOSITE_DIR MINREAD THRESHOLD
