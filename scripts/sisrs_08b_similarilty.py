#!/usr/bin/env python3

import csv
import argparse
import sys
from Bio import SeqIO
from os import system
from math import ceil

def blast_and_filter(outpath):
    '''
    get the list of loci to drop from the contigs based on blast overlap of 90% of the smaller contig

    Args:
        outpath (string): output folder for work

    Returns:
        string: fasta file contigs_for_probes
        list: loci to drop from dataset
    '''

    fasta = outpath + "/SISRS_Run/contigs_for_probes.fa"
    matrix_file = outpath + "/SISRS_Run/contig_similarity.txt"

    blast_command = f"blastn -query {fasta} -subject {fasta} -task blastn -outfmt \"6 qseqid sseqid qlen slen length pident nident\" > {matrix_file}"
    # print(blast_command)
    system(blast_command)

    loci_to_drop = []

    print('The following are blastn results for the locus pairs where one should be dropped')
    print('Query_ID Subject_ID Query_len Subject_len Overlap_length %ident #ident')

    with open(matrix_file, newline='') as file:  # qseqid sseqid qlen slen length pident nident
        contigs = csv.reader(file, delimiter='\t')  # iterator
        for line in contigs:
            if line[0] != line[1]:  # not the same contig
                if float(line[5]) > 90:  # 90% identity for region that overlaps
                    if float(line[4]) > 0.9 * min(float(line[2]), float(
                            line[3])):  # length of overlap is greater than 90% of the smaller contig
                        print(line)
                        if int(line[2]) > int(line[3]):
                            loci_to_drop.append(line[1])
                        else:
                            loci_to_drop.append(line[0])

    loci_to_drop = list(set(loci_to_drop))
    print('The following loci are dropped due to blast overlap', loci_to_drop)

    return(fasta, loci_to_drop)


def heterozyg_filter(outpath, loci_to_drop):
    '''
    finds loci with high heterozygosity based on prior output file (_hztable.csv) and adds them to the list - these might be paralogs

    Args:
        outpath (string): output folder for work
        loci_to_drop (list): loci to drop from dataset

    Returns:
    list: updated loci to drop from dataset

    '''

    # list of taxa
    with open(outpath+'/TaxonList.txt') as f:
        taxonlist = f.readlines()
        taxonlist = sorted([x.rstrip() for x in taxonlist])

    hz_dict = {} #keep track of which contigs have heterozygosity (= paralogs?)
    for taxon in taxonlist:
        f = outpath + 'SISRS_Run/' + taxon + '/' + taxon + '_hztable.csv' #het table output by step 7
        with open(f, newline = '') as fi:
            hz_table = csv.reader(fi, delimiter=',') #iterator: SISRS_contig-5000011,0.09424083769633508
            for line in hz_table:
                #print(line)
                if line[1] != 'NA': #seq not found in sp
                    if float(line[1]) > 0.01: #high het
                        if line[0] in hz_dict:
                            hz_dict[line[0]].append(taxon)
                        else:
                            hz_dict[line[0]] = [taxon]

    for k, v in hz_dict.items():
        if len(v) > 0.5 * len(taxonlist): #half taxon list may have paralogs - this is arbitrary
            loci_to_drop.append(k)

    loci_to_drop = list(set(loci_to_drop))
    print('The following loci are dropped due to blast overlap or excess heterozygosity (potential paralogy)', loci_to_drop)

    return(loci_to_drop)

def filter_probes(fasta, loci_to_drop, final_total):
    '''

    Args:
        fasta (string): fasta file contigs_for_probes
        loci_to_drop (list): loci to drop from dataset
        final_total: number of loci to keep for probes - may be ALL or an int

    Returns:
        none

    '''

    contig_seqs = []
    total_length_probes = 0
    total_tiles_80 = 0
    total_tiles_60 = 0
    with open(fasta) as handle: #contigs_for_probes
        for record in SeqIO.parse(handle, "fasta"):
            if record.id not in loci_to_drop:   #keep contig if not in list to drop
                contig_seqs.append(record)
                total_length_probes += len(str(record.seq))
                total_tiles_80 += ceil(len(str(record.seq)) / 80)
                total_tiles_60 += ceil(len(str(record.seq)) / 60)

    print('The total length of probes is', total_length_probes)
    print('The expected number of tiles of size 60 is', total_tiles_60)
    print('The expected number of tiles of size 80 is', total_tiles_80)

    contig_seqs.sort(key=lambda r: -len(r))

    if final_total == 'ALL':
        SeqIO.write(contig_seqs, fasta + '.fasta', "fasta")
    else:
        final_total = int(final_total)
        SeqIO.write(contig_seqs[0:final_total], fasta + '.fasta', "fasta")

def all_of_8b(outpath, final_total):

    fasta, loci_to_drop = blast_and_filter(outpath)

    loci_to_drop = heterozyg_filter(outpath, loci_to_drop)

    filter_probes(fasta, loci_to_drop, final_total)

if __name__ == '__main__':
    # Get arguments
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-d', '--outpath', action='store', nargs="?")
    my_parser.add_argument('-m', '--final_total', action='store', default='ALL')
    args = my_parser.parse_args()

    all_of_8b(args.outpath, args.final_total)