#!/usr/bin/env python3

import argparse
from collections import Counter
from Bio import SeqIO, AlignIO

def count_snps(fa, aligned_contigs2, outf):
    numsnps = 0
    num_pi_sites = 0
    num_pi_sites_dict = {}

    with open(fa) as contigs:
        for record in SeqIO.parse(contigs, "fasta"):
            contig_num_pi_sites = 0
            alignment = AlignIO.read(aligned_contigs2+record.id+'.fasta', "fasta")
            for i in range(alignment.get_alignment_length()):
                site = list(alignment[:, i]) #get site for all spp
                site = [s for s in site if s in ['a', 'c', 'g', 't']] #remove gaps
                if len(site) > 3: #present in 4 spp
                    if len(list(set(site))) > 1: #variation (although singletons allowed)
                        numsnps += 1
                        if Counter(site).most_common(2)[1][1] > 1: # the number of times the second-most common item is found
                            num_pi_sites += 1 #only keeping non-singletons
                            contig_num_pi_sites += 1
                            #print("".join(site))
            num_pi_sites_dict[record.id] = contig_num_pi_sites
    print(num_pi_sites)

    outfile = open(outf, 'w')
    for item in sorted(num_pi_sites_dict.items(), key=lambda item: item[1]):
        outfile.write(item[0] + ' ' + str(item[1]) + '\n') #sort dict by number of sites and write to a file
    outfile.close()

    return(num_pi_sites)

def main_count(folder):
    fa = folder + '/SISRS_Run/contigs_for_probes.fa.fasta'
    aligned_contigs2 = folder+'/SISRS_Run/aligned_contigs2/'
    outf = folder + '/SISRS_Run/variation_in_contigs.txt'
    count_snps(fa, aligned_contigs2, outf)

if __name__ == '__main__':

    # Get arguments
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-f', '--folder', action='store', nargs="?")
    args = my_parser.parse_args()

    main_count(args.folder)