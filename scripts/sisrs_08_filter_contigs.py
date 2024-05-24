#!/usr/bin/env python3

import csv
import argparse
import sys
from collections import Counter
from Bio import SeqIO
from os import path, mkdir, listdir, system
from math import ceil


def alignContigs(new_contig_folder, newnew_contig_folder, processors):
    '''
    Align contigs and output to folder
    Args:
        new_contig_folder (string): outPath + "/SISRS_Run/contigs_outputs2/"
        newnew_contig_folder (string): outPath + "/SISRS_Run/aligned_contigs2/"
        processors: num processors

    Returns:
        none
    '''
    contig_list = listdir(new_contig_folder)  # list of contigs
    for con in contig_list:
        con_p_in = path.join(new_contig_folder, con)
        con_p_out = path.join(newnew_contig_folder, con)
        mafft_command = f"mafft --auto --thread {processors} {con_p_in} > {con_p_out}"
        system(mafft_command)

def write_alignment_plus_composite(high_count_contigs, outPath, num_sp, composite, new_contig_folder):
    """

    Args:
        high_count_contigs (list): contigs that have more than the min_threshold (number of variable sites) and are longer than length_of_locus
        outPath (str): where we are outputting data e.g. '../../SISRS_Small_test'
        num_sp (int): number of species required to be present given number of missing allowed (default is half)
        composite (dict):  SeqIO.to_dict(SeqIO.parse(outPath + "/SISRS_Run/Composite_Genome/contigs.fa", "fasta"))
        new_contig_folder (string): outPath + "/SISRS_Run/contigs_outputs2/"

    Returns:
        none

    """

    for k in high_count_contigs:
        a_file = outPath + '/SISRS_Run/contigs_outputs/' + k + '.fasta'
        if path.exists(a_file):
            alignment = list(SeqIO.parse(a_file, "fasta"))
            if len(alignment) == num_sp: #check we have all spp
                composite_seq = composite[k] #get composite seq for this contig
                alignment.append(composite_seq)
                SeqIO.write(alignment, new_contig_folder + k + '.fasta', "fasta")

def print_probe_info(outPath, good_contigs, composite):
    '''
    write the composite sequences to contigs_for_probes.fa
    print out probe info

    Args:
        outPath (str): where we are outputting data e.g. '../../SISRS_Small_test'
        good_contigs (dict): contigs where all seqs are less than the max dist from the composite
        composite (dict):  SeqIO.to_dict(SeqIO.parse(outPath + "/SISRS_Run/Composite_Genome/contigs.fa", "fasta"))

    Returns:
        none
    '''

    total_length_probes = 0
    total_tiles_80 = 0
    total_tiles_60 = 0
    with open(outPath + "/SISRS_Run/contigs_for_probes.fa", 'w') as f:
        for k, v in sorted(good_contigs.items(), key=lambda item: item[1]):  # line in good_contigs.items():
            print(v)
            f.write('>' + k + "\n")
            f.write(str(composite[k].seq) + "\n")
            total_length_probes += len(str(composite[k].seq))
            total_tiles_80 += ceil(len(str(composite[k].seq)) / 80)
            total_tiles_60 += ceil(len(str(composite[k].seq)) / 60)

    print(total_length_probes)
    print(total_tiles_60)
    print(total_tiles_80)

def filter_contigs_distance(high_count_contigs, newnew_contig_folder, taxon_list, max_dist, contigcounts):
    '''

    Args:
        high_count_contigs (list): contigs that have more than the min_threshold (number of variable sites) and are longer than length_of_locus:
        newnew_contig_folder (string): outPath + "/SISRS_Run/aligned_contigs2/":
        taxon_list (list): list of taxa
        max_dist (float): maximum % difference allowed for sample seqs from composite seq
        contigcounts (dict): contig: number of variable sites - comes from the number of times a contig is seen in the specified alignment file

    Returns:
        dict: contigs where all sample sequences are within the max_dist of the composite seq
    '''

    good_contigs = {}
    for k in high_count_contigs:
        distances = {}
        a_file = newnew_contig_folder + k + '.fasta'
        if path.exists(a_file):
            alignment = SeqIO.to_dict(SeqIO.parse(a_file, "fasta"))
            composite_seq = str(alignment[k].seq).upper()

            for sp, seq in alignment.items():
                if sp in taxon_list: #not composite
                    seq2 = str(seq.seq).upper()
                    #print(seq2)
                    #print(str(composite_seq.seq))
                    distances[sp] = 0
                    if len(seq2) != len(composite_seq):
                        print(composite_seq)
                        print(seq2)
                        sys.exit("Seqs don't match")
                    for i in range(len(seq2)): #go through seq site by site to get distance
                        if seq2[i] != '-':
                            if composite_seq[i] != seq2[i]:
                                distances[sp]+=1
                    distances[sp] = distances[sp] / len(seq2.replace('-', ''))
                    #print(distances[sp])
            if max(distances.values()) < max_dist:
                good_contigs[k] = contigcounts[k]
            else:
                print(sorted(distances.values()))
    return(good_contigs)

def all_of_step8(outPath, alignment_file, min_threshold, processors, length_of_locus, max_dist, num_miss):

    path_to_seq_lengths = outPath + '/SISRS_Run/Composite_Genome/contigs_SeqLength.tsv'
    path_to_sites = outPath + '/SISRS_Run/' + alignment_file #name of file eg alignment_pi_locs_m2.txt - use something minimally missing
    path_to_output = outPath + '/SISRS_Run/'

    #read in informative site names
    with open(path_to_sites) as file:
        sites = dict(csv.reader(file, delimiter='/')) #SISRS_contig-2/44

    #read in contig names and their overall lengths
    with open(path_to_seq_lengths, newline = '') as file:
        contigs = dict(csv.reader(file, delimiter='\t')) #SISRS_contig-2    351

    contigcounts = dict(Counter(sites.keys())) #number of variable sites
    high_count_contigs = [k for k, v in contigcounts.items() if int(v) > min_threshold and int(contigs[k]) > length_of_locus]

    composite = SeqIO.to_dict(SeqIO.parse(outPath + "/SISRS_Run/Composite_Genome/contigs.fa", "fasta"))     #>SISRS_contig-2000001 186 nucleotides

    with open(outPath+'/TaxonList.txt') as f:
        taxon_list = f.readlines()
        taxon_list = sorted([x.rstrip() for x in taxon_list])

    num_sp = len(taxon_list)
    if num_miss == -1:
        num_miss = int(0.5*num_sp)
    num_sp = num_sp - num_miss


    #see https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.5260

    #goals:
    #1. at least 160 bp = 2 probes
    #2. each sp no more than 7% different from the composite
    #3. 60k probes total = 4m bp
    #all sp in contig alignment
    #less than X% missing per site - info for each spp

    new_contig_folder = outPath + "/SISRS_Run/contigs_outputs2/"
    if not path.exists(new_contig_folder):
        mkdir(new_contig_folder)

    newnew_contig_folder = outPath + "/SISRS_Run/aligned_contigs2/"
    if not path.exists(newnew_contig_folder):
        mkdir(newnew_contig_folder)

    write_alignment_plus_composite(high_count_contigs, outPath, num_sp, composite, new_contig_folder)

    alignContigs(new_contig_folder, newnew_contig_folder, processors)

    #good ideas here: http://www.mossmatters.com/blog/SequenceClusters.html

    good_contigs = filter_contigs_distance(high_count_contigs, newnew_contig_folder, taxon_list, max_dist, contigcounts)

    print_probe_info(outPath, good_contigs, composite)

if __name__ == '__main__':

    # Get arguments
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-d', '--outPath', action='store', nargs="?")
    my_parser.add_argument('-a', '--alignment', action='store', nargs="?")
    my_parser.add_argument('-t', '--min_threshold', type=int, action='store', default=1, nargs="?")
    my_parser.add_argument('-p','--processors', type=int, action='store', default=1,nargs="?")
    my_parser.add_argument('-l', '--length_of_locus', type=int, action='store', default=160, nargs="?")
    my_parser.add_argument('-m', '--max_dist', type=float, action='store', default=0.08, nargs="?")
    my_parser.add_argument('-n', '--num_miss', type=int, action='store', default=-1, nargs="?") #will later convert a default of -1 to 50% of samples
    args = my_parser.parse_args()

    outPath = args.outPath  #outPath = '../../SISRS_Small_test'
    alignment = args.alignment
    min_threshold = args.min_threshold #number of variable sites
    processors = args.processors
    length_of_locus = args.length_of_locus #160 = 2 probes
    max_dist = args.max_dist #probably 0.08 or less than 8% difference from composite
    num_miss = args.num_miss #number missing in locus alignment

    all_of_step8(outPath, alignment, min_threshold, processors, length_of_locus, max_dist, num_miss)