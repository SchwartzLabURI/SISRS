import sys
import os

#This script acts as a wrapper for get_organelle_from_reads.py. It takes three arguments: the path to the forward read, the path to
#the reverse read, and the path to output to. Input paths can be gzipped or not, and the output path is created upon calling.

#After chloroplasts have been identified, bbmap is called to filter the plastome from the reference genome
#

forward_end, reverse_end, output = sys.argv[1], sys.argv[2], sys.argv[3]

#These modules are required
#sympy/1.7.1-foss-2020b sympy/1.8-foss-2021a
#scikit-learn/0.23.2-fosscuda-2020b
#SPAdes/3.15.2-GCC-10.2.0
#Bowtie2/2.4.4-GCC-10.2.0
#BLAST+/2.8.1-foss-2018b

os.mkdir(f"{output}")

find_chloroplasts = f"python3 GetOrganelle/get_organelle_from_reads.py -1 {forward_end} -2 {reverse_end} -o {output}/GetChloroplast/ -F embplant_pt"
os.system(find_chloroplasts)


name="embplant_pt.K115.complete.graph1.1.path_sequence.fasta"

filter_chloroplasts = f"bbmap.sh in={forward_end} in2={reverse_end} ref={output}/GetChloroplast/{name} ambiguous=all threads=36 outm={output}/removed_chloro/mappedChloroplast1.fastq outm2={output}/removed_chloro/mappedChloroplast2.fastq outu={output}/removed_chloro/SRR7716091__Nuclear_Trim_1.fastq.gz outu2={output}/removed_chloro/SRR7716091__Nuclear_Trim_2.fastq.gz"
os.system(filter_chloroplasts)
