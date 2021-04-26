#!/usr/bin/env python3
import os
import sys
import subprocess
from subprocess import call
from cmdCheck import *
from Bio.Seq import Seq
from Bio import SeqIO
from os import listdir
from os.path import isfile, join


#format the file outputs
def format_consensus_output(output_path):
    path_to_contigs_LocList_file = output_path + 'Composite_Genome/'
    path_to_taxon_dirs = output_path
    #get list of taxon dirs except Composite_Genome
    taxon_list = []
    taxon_list = os.listdir(path_to_taxon_dirs)
    taxon_list.remove('Composite_Genome') #remove this directory from the list for easy traversal

    contigs_LocList_file = path_to_contigs_LocList_file + 'contigs_LocList'

    #get names of contigs for naming of the files and fasta headers
    contigs_dict = {} #store contigs in a dictionary
    with open(contigs_LocList_file, 'r') as in_file:
        for line in in_file:
            split1_line = line.strip().split('-')
            split2_line = split1_line[1].split('/')

            key = split2_line[0]
            val = split2_line[1]
            contigs_dict[key] = val

    in_file.close()

    counter = 1
    for key in sorted(contigs_dict.keys()):
        #store each sequence for each taxon by contig name
        contigs_file = path_to_contigs_LocList_file + '/' + 'SISRS_contig-' + key + '.fasta'
        out_file = open(contigs_file,'a+')

        for dir in taxon_list:
            out_file.write(">" + dir + '\n')

            #read from the consensus sequence file by taxon
            taxon_consensus_file = path_to_taxon_dirs + '/' + dir + '/' + dir + '_single_line_format_consensus.fa'

            in_file = open(taxon_consensus_file, 'r')
            all_lines = in_file.readlines()

            #"interlace" sequences by taxons into a file by contig position
            out_file.write(all_lines[counter] + '\n')
        if(counter + 2 <= len(all_lines) - 1):
            counter += 2
        else:
            break;

    in_file.close()
    out_file.close()

#remove sequences that contain only N's
def remove_Ns_contigs(path_to_contigs, taxa_threshold):

    if (path_to_contigs[-1] != '/'):
        path_to_contigs = path_to_contigs + '/'

    list_of_contig_files = [f for f in listdir(path_to_contigs) if isfile(join(path_to_contigs, f))]
    #list_of_contig_files.remove('contigs_LocList')

    list_of_contig_files.sort() #sort list to make sure to start from the first contig


    for i in range(len(list_of_contig_files)):
        new_seqs=[]
        for seq_record in SeqIO.parse(path_to_contigs + list_of_contig_files[i], "fasta"):
            if seq_record.seq.count("N") != len(seq_record):
                new_seqs.append(seq_record)

        file_out = path_to_contigs + list_of_contig_files[i] #will overwrite the old file
        if(len(new_seqs) >= taxa_threshold): #if the number of sequences is below or equals given threshold overwrite the file
            SeqIO.write(new_seqs, file_out, "fasta")
        else:
            os.remove(file_out) #otherwise, delete that file from the path



if __name__ == '__main__':

    # Store the command line in a seperate argument
    cmd = sys.argv

    if len(cmd) < 2:
        print("THIS SCRIPT REQUIERS 2 ARGUMENTS (-dir <path to output data directory>)")
        exit()

    out_dir = ""
    sisrs = ""

    if '-dir' in cmd or '--directory' in cmd:
        out_dir = isFound('-dir','--directory',cmd)
        sisrs = out_dir
        if (sisrs[-1] == '/'):
            output_path = sisrs + "SISRS_Run/"
        else:
            output_path = sisrs + "/SISRS_Run/"

        contigs_from_consensus = output_path + contigs_outputs
        print(contigs_from_consensus)
        os.mkdir(contigs_from_consensus)#creates a directory "contigs_outputs"
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
    subprocess.call("./mask_ref_seq.sh " + output_path, shell=True)

    #generate consensus sequence
    os.chmod('contigs_driver.sh', 0o755)
    subprocess.call("./contigs_driver.sh " + output_path, shell=True)


    format_consensus_output(output_path)
    remove_Ns_contigs(output_path, taxa_threshold)
