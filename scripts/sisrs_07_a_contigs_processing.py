#!/usr/bin/env python3


'''
Last edit: Yana Hrytsenko & Alexander Knyshov May 6th, 2021
'''
import os
import sys
import subprocess
from subprocess import call
from cmdCheck import *
from Bio.Seq import Seq
from Bio import SeqIO
from os import listdir
from os.path import isfile, join


# format the file outputs
def format_consensus_output(output_path, taxa_threshold):
    ''' This function formats the consensus files outputs.  '''

    path_to_contigs_LocList_file = output_path + 'SISRS_Run/Composite_Genome/'
    path_to_taxon_dirs = output_path + "Reads/RawReads/"


    # get a list of taxon dirs except Composite_Genome
    taxon_list = []
    taxon_list = os.listdir(path_to_taxon_dirs)
    taxon_list.remove('fastqcOutput') # remove this directory from the list for easy traversal
    taxon_list.remove('trimOutput') # remove this directory from the list for easy traversal

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
            # save only if not full of N - with additional strip
            if len(line.replace("N", "").strip()) > 0:
                contigs_dict[taxon_list[0]] = line.strip()

            # read 1 line from all other taxa and save to the same dict but only if not complete N
            for dir in taxon_list[1:]:
                taxonseq = consensus_handles[dir].readline().strip()
                # only save if after removing all N there's something left - have to do additional strip
                if len(taxonseq.replace("N", "").strip()) > 0:
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

if __name__ == '__main__':

    # Store the command line in a seperate argument
    cmd = sys.argv

    if len(cmd) < 4:
        print("THIS SCRIPT REQUIERS 4 ARGUMENTS (-trh <taxon threshold> -dir <path to output data directory>)")
        exit()

    out_dir = ""
    sisrs = ""

    if '-dir' in cmd or '--directory' in cmd:
        out_dir = isFound('-dir','--directory',cmd)
        sisrs = out_dir
        if (sisrs[-1] != '/'):
            output_path = sisrs + '/'
        else:
            output_path = sisrs

        contigs_from_consensus = output_path + "SISRS_Run/" + "contigs_outputs"

        if not os.path.exists(contigs_from_consensus):
            os.mkdir(contigs_from_consensus)
    else:
        print("SPECIFY THE OUTPUT PATH (-dir, --directory). PROGRAM EXITING.")
        exit()

    if '-trh' in cmd or '--threshold' in cmd:
        num_taxa = isFound('-trh','--threshold',cmd)
        taxa_threshold = int(num_taxa)
    else:
        print("SPECIFY THE TAXA THRESHOLD (-trh, --threshold). PROGRAM EXITING.")
        exit()


    #first mask the ref sequence
    os.chmod('mask_ref_seq.sh', 0o755)
    subprocess.call("./mask_ref_seq.sh " + output_path + 'SISRS_Run/', shell=True)

    #generate consensus sequence
    os.chmod('contigs_driver.sh', 0o755)
    subprocess.call("./contigs_driver.sh " + output_path + 'SISRS_Run/', shell=True)

    format_consensus_output(output_path, taxa_threshold)
