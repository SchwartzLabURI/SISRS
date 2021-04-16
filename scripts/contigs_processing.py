#!/usr/bin/env python3
import os
import sys
import subprocess
from subprocess import call
#from cmdCheck import *

'''
# Store the command line in a seperate argument
cmd = sys.argv

if len(cmd) < 2:
    print("THIS SCRIPT REQUIERS 2 ARGUMENTS (-dir <path to output data directory>)")
    exit()

out_dir = ""
sisrs = ""

if '-dir' in cmd or '--directory' in cmd:
    out_dir = isFound('-dir','--directory',cmd)
    sisrs = os.path.dirname(out_dir)
else:
    print("SPECIFY THE OUTPUT PATH (-dir, --directory).PROGRAM EXITING.")
    exit()

#output path
sisrs = sys.argv[1]
output_path = sisrs + "/SISRS_Run/"

#first mask the ref sequence
os.chmod('mask_ref_seq.sh', 0o755)
subprocess.call("./mask_ref_seq.sh" + output_path, shell=True)

#generate consensus sequence
os.chmod('contigs_driver.sh', 0o755)
subprocess.call("./contigs_driver.sh " + output_path, shell=True)
'''

sisrs = sys.argv[1]

output_path = sisrs + "/SISRS_Run/"
#format the file outputs
path_to_contigs_LocList_file = output_path + 'Composite_Genome/'
path_to_taxon_dirs = output_path #/data3/schwartzlab/yana/SISRS_reRun/SISRS_Walkthrough_test_BOBs_scripts/SISRS_Run/
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



for key in sorted(contigs_dict.keys()):

    contigs_file = path_to_contigs_LocList_file + '/' + 'SISRS_contig-' + key + '.fasta'
    out_file = open(contigs_file,'a+')

    for dir in taxon_list:
        out_file.write(">" + dir + '\n') #>AotNan
        taxon_consensus_file = path_to_taxon_dirs + '/' + dir + '/' + dir + '_single_line_format_consensus.fa'

        in_file = open(taxon_consensus_file, 'r')
        header_line_f = in_file.readline() #>SISRS_contig-0 discard this line
        sequence_line_f = in_file.readline() #NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
        out_file.write(sequence_line_f + '\n')

    in_file.close()
    out_file.close()
