#!/usr/bin/env python3

import sys
import os
import glob
from sisrs_02_read_trimmer import trimHelper

#This script acts as a wrapper for get_organelle_from_reads.py. It takes three arguments: the path to the forward read, the path to
#the reverse read, and the path to output to. Input paths can be gzipped or not, and the output path is created upon calling.

#After chloroplasts have been identified, bbmap is called to filter the plastome from the reference genome
#

if __name__ == '__main__':

    fo = sys.argv[1]
    taxon = sys.argv[2]
    proc = sys.argv[3]

    folder = fo+"/Reads/TrimReads/"+taxon
    output = fo+"/Reads/TrimReads_nocl/"+taxon
    #forward_end, reverse_end, output = sys.argv[1], sys.argv[2], sys.argv[3]
	
    os.makedirs(output) #make dir and enclosing folders

    #get paired reads
    tax_dir = folder #has fastqs
    trim_read_dir = folder
    newData = [] 
    out_trim_dir,left_pairs,right_pairs,single_end = trimHelper(tax_dir,trim_read_dir,newData)

    for left in left_pairs:
        forward_end = left +'_Trim_1.fastq.gz'
        reverse_end = left +'_Trim_2.fastq.gz'
        #at some point fix this so we don't call python from python
        find_chloroplasts = f"python3 GetOrganelle/get_organelle_from_reads.py -1 {forward_end} -2 {reverse_end} -o {output}/ -F embplant_pt"
        os.system(find_chloroplasts)

        cl_genome = glob.glob(output+'/*.fasta') #find fasta in output folder  #name="embplant_pt.K115.complete.graph1.1.path_sequence.fasta"

        filter_chloroplasts = f"bbmap.sh in={forward_end} in2={reverse_end} ref={output}/{cl_genome[1]} ambiguous=all threads={proc} outm={output}/mappedChloroplast1.fastq outm2={output}/mappedChloroplast2.fastq outu={output}/{left}_Nuclear_Trim_1.fastq.gz outu2={output}/{left}_Nuclear_Trim_2.fastq.gz"
        os.system(filter_chloroplasts)

    for singlef in single_end:
        reads = singlef + '_Trim.fastq.gz'
        find_chloroplasts = f"python3 GetOrganelle/get_organelle_from_reads.py -u {reads} -o {output}/ -F embplant_pt"
        os.system(find_chloroplasts)

	cl_genome = glob.glob(output+'/*.fasta')

	filter_chloroplasts = f"bbmap.sh in={reads} ref={output}/{cl_genome[1]} ambiguous=all threads={proc} outm={output}/mappedChloroplast.fastq outu={output}/{singlef}_Nuclear_Trim.fastq.gz
