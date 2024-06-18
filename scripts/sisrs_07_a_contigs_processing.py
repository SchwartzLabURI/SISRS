#!/usr/bin/env python3

import os
import sys
import gzip
import subprocess
import argparse
from subprocess import call
from Bio.Seq import Seq
from Bio import SeqIO
from os import listdir
from os.path import isfile, join
import statistics
from sisrs_08_filter_contigs import get_taxon_list

def format_consensus_output(output_path, taxa_threshold):
    '''
    This function formats the consensus files outputs.

    Arguments: path to the output directory, taxa threshold.

    Returns: none.
    '''

    #where to put contigs
    p = output_path + 'SISRS_Run/contigs_outputs/'
    if not os.path.exists(p):
        os.makedirs(p)

    # list of taxa
    taxon_list = get_taxon_list(output_path+'/TaxonList.txt')

    # open all taxa handles at once
    consensus_handles = {}
    for dir in taxon_list:
        taxon_consensus_file = output_path + 'SISRS_Run/' + dir + '/' + dir + '_single_line_format_consensus.fa'
        consensus_handles[dir] = open(taxon_consensus_file, 'r')

    # iterate over any (the first taxon) consensus, other files assumed to have exactly same order and contents
    for line in consensus_handles[taxon_list[0]]:
        # if the line starts with ">" then it's a contig name, record it
        if line[0] == ">":
            contigName = line[1:].strip()
            # just a sanity check - get same line from other taxa to verify that contig is same
            for other_taxon in taxon_list[1:]:
                if contigName != consensus_handles[other_taxon].readline()[1:].strip():
                    print ("contig order must be different between consensus files")
                    sys.exit()
            # if no check is performed then other handles should be iterated for 1 step, perhaps by invoking .next
            # to match the first handle which we iterate in the for loop
        # if the line does not start with ">" then we assume it's a sequence following the contigName - save for all taxa
        else:
            # create a dictionary to temporary hold the sequences for this locus
            contigs_dict = {}

            # save the first taxon (over which we iterate in the for loop) - we already have this line in RAM
            # remove N and save only if not empty
            taxon0seq = line.replace("N", "").strip()
            if len(taxon0seq) > 0:
                contigs_dict[taxon_list[0]] = taxon0seq

            # read 1 line from all other taxa and save to the same dict
            for dir in taxon_list[1:]:
                taxonseq = consensus_handles[dir].readline().replace("N", "").strip()
                # only save if after removing all N there's something left
                if len(taxonseq) > 0:
                    contigs_dict[dir] = taxonseq

            # check how many taxa passed
            if len(contigs_dict) >= taxa_threshold:
                # if passing the threshold, create a file
                contigs_file = output_path + 'SISRS_Run/contigs_outputs/' + contigName + '.fasta'
                out_file = open(contigs_file,'w')

                # populate the file
                for dir in sorted(contigs_dict):
                    out_file.write(">" + dir + '\n')
                    out_file.write(contigs_dict[dir] + '\n')

                # close the output handle
                out_file.close()

    # close all input handles
    for dir in taxon_list:
        consensus_handles[dir].close()

# get vcf-based consensus
def vcf_consensus(output_path, coverage_threshold, hz_threshold):
    '''
    This function produces the consensus sequences from VCF data.
    The function assumes all sites (including non variable) were
    output.
    First, get contig names from the header, this is needed since
    even with all sites output, loci without any reads are completely
    missing from the VCF body.
    Then iterate over VCF body to get sites, keep track of current
    contig in the body vs in the header, fill in empty contigs (loci),
    for non empty, filter sites with low allelic coverage and
    positions marked with '.', as well as filter out (replace with N)
    sequences with high heterozygosity.

    Arguments: path to the output directory, coverage threshold, heterozygosity threshold.

    Returns: none.
    '''

    # list of taxa
    taxon_list = get_taxon_list(output_path + '/TaxonList.txt')

    # open all taxa handles at once
    consensus_handles = {}
    for directory in taxon_list:
        taxon_vcf_file = output_path + 'SISRS_Run/' + directory + \
            '/' + directory + '.vcf.gz'
        taxon_consensus_file = output_path + 'SISRS_Run/' + directory + \
            '/' + directory + '_single_line_format_consensus.fa'
        infile = gzip.open(taxon_vcf_file, 'rt')
        outfile = open(taxon_consensus_file, "w")
        outputseq = ""
        locname = ""
        loclen = 0
        hznum = 0
        #data on heterozygosity - for debugging purposes
        hzhandle = open(output_path + 'SISRS_Run/' + directory + \
            '/' + directory +'_hztable.csv', 'w')
        sisrs_contig_list = []
        sisrs_contig_counter = 0
        #iterate over VCF
        for line in infile:
            if line[0] == "#":
                if line[0:13] == "##contig=<ID=":
                    #get info on all contigs (loci) to detect completely
                    #missing loci in the VCF
                    sisrs_contig_list.append(line[13:].split(',')[0])
                else:
                    #ignore most of the header
                    continue
            else:
                splitline = line.strip().split("\t")
                if locname != splitline[0] and locname != "":
                    #next locus in VCF - output the previous
                    while sisrs_contig_list[sisrs_contig_counter] != locname:
                        #if current locus in the header is diffent,
                        #produce empty (N) output for loci that were up to this
                        #point missing in the VCF body
                        print (">"+sisrs_contig_list[sisrs_contig_counter], file=outfile)
                        print ("N", file=outfile)
                        print (sisrs_contig_list[sisrs_contig_counter]+",NA", file=hzhandle)
                        sisrs_contig_counter += 1
                    # filter by heterozygosity
                    if loclen > 0:
                        hzval = hznum/loclen
                        print (locname+","+str(hznum/loclen), file=hzhandle)
                        if hzval > hz_threshold:
                            outputseq = "N"
                    else:
                        print (locname+",NA", file=hzhandle)
                        outputseq = "N"
                    # write output and reset vars
                    print (">"+locname, file=outfile)
                    print (outputseq, file=outfile)
                    sisrs_contig_counter += 1
                    outputseq = ""
                    hznum = 0
                    loclen = 0
                #actual VCF line parsing and saving
                locname = splitline[0]
                basecall = [splitline[3]] + splitline[4].split(",")
                info = splitline[8].split(":")
                ad_pos = info.index("AD")
                coverage_list = [int(x) for x in splitline[9].split(":")[ad_pos].split(",")]
                max_value = max(coverage_list)
                #skip sites with low cov or '.'
                if max_value >= coverage_threshold and basecall[1] != ".":
                    max_index = coverage_list.index(max_value)
                    outputseq += basecall[max_index]
                    if len(basecall) > 2:
                        hznum += 1
                    loclen += 1
        #produce the output for the final locus
        while sisrs_contig_list[sisrs_contig_counter] != locname:
            print(">"+sisrs_contig_list[sisrs_contig_counter], file=outfile)
            print("N", file=outfile)
            print(sisrs_contig_list[sisrs_contig_counter]+",NA", file=hzhandle)
            sisrs_contig_counter += 1
        if loclen > 0:
            hzval = hznum/loclen
            print (locname+","+str(hznum/loclen), file=hzhandle)
            if hzval > hz_threshold:
                outputseq = "N"
        else:
            print(locname+",NA", file=hzhandle)
            outputseq = "N"
        print(">"+locname, file=outfile)
        print(outputseq, file=outfile)
        sisrs_contig_counter += 1
        #given the contig list in the header, fill in the remaining missing loci
        #in the VCF body
        while sisrs_contig_counter < len(sisrs_contig_list):
            print(">"+sisrs_contig_list[sisrs_contig_counter], file=outfile)
            print("N", file=outfile)
            print(sisrs_contig_list[sisrs_contig_counter]+",NA", file=hzhandle)
            sisrs_contig_counter += 1
        #close handles
        infile.close()
        outfile.close()
        hzhandle.close()

def contig_driver(output_path):
    '''
    This funcction is designed to call mpileup command and perform pileups on each of the bam files.

    Arguments: path to the output directory.

    Returns: none.

    '''

    # list of taxa
    taxon_list = get_taxon_list(output_path + '/TaxonList.txt')

    for tax in taxon_list:
        bam_file = output_path+'SISRS_Run/'+tax+'/'+tax+'.bam'
        vcf_zipped = output_path+'SISRS_Run/'+tax+'/'+tax+'.vcf.gz'

        #do pileup using no reference and report allelic depth
        #then call variants using all sites (M) and zip the output
        bcf_command = f"bcftools mpileup -Ou --no-reference -a FORMAT/AD {bam_file} | bcftools call -Oz -mM -o {vcf_zipped}"
        os.system(bcf_command)

def run7a(output_path, taxa_threshold, cov_thresh, hz_thresh):
    contig_driver(output_path)

    #call consensus and filter by allelic coverage and ratio of heterozygous sites
    vcf_consensus(output_path, coverage_threshold=cov_thresh, hz_threshold=hz_thresh)

    #recompile sequence of each taxon by locus
    format_consensus_output(output_path, taxa_threshold)

if __name__ == '__main__':

    # python3 sisrs_07_a_contigs_processing.py -t $T -dir $OUTFOLDER
    # Get arguments
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-t','--threshold',action='store',nargs="?", type = int)
    my_parser.add_argument('-d', '--dir', action='store', nargs="?")
    my_parser.add_argument('-c', '--cov', action='store', default=3,nargs="?", type=int)
    my_parser.add_argument('-z', '--hz', action='store', default=0.01, nargs="?", type=float)

    args = my_parser.parse_args()

    run7a(args.dir +'/', args.threshold, args.cov, args.hz)
