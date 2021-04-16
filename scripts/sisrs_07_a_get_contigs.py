import os
import sys
import itertools


path_to_contigs_LocList_file = "/Users/yanahrytsenko/Documents/GitHub/SISRS/SISRS_Run_yana/Composite_Genome" #sys.argv[1] #/data3/schwartzlab/yana/SISRS_reRun/SISRS_Walkthrough_test_BOBs_scripts/SISRS_Run/Composite_Genome/
path_to_taxon_dirs = "/Users/yanahrytsenko/Documents/GitHub/SISRS/SISRS_Run_yana" #sys.argv[2] #/data3/schwartzlab/yana/SISRS_reRun/SISRS_Walkthrough_test_BOBs_scripts/SISRS_Run

#get list of taxon dirs except Composite_Genome
taxon_list = []
taxon_list = os.listdir(path_to_taxon_dirs)


taxon_list.remove('Composite_Genome') #remove this directory from the list for easy traversal


#put contigs into dict
contigs_LocList_file = path_to_contigs_LocList_file + '/' + 'contigs_LocList'

contigs_dict = {} #store contigs in a dictionary
with open(contigs_LocList_file, 'r') as in_file:
    for line in in_file:
        split1_line = line.strip().split('-')
        split2_line = split1_line[1].split('/')

        key = split2_line[0]
        val = split2_line[1]
        contigs_dict[key] = val

in_file.close()

i = 0

for key in sorted(contigs_dict.keys()):

    contigs_file = path_to_contigs_LocList_file + '/' + 'SISRS_contig-' + key + '.fasta'

    out_file = open(contigs_file,'a+')

    for dir in taxon_list:

        out_file.write(">" + dir + '\n') #>AotNan

        taxon_LocList_file = path_to_taxon_dirs + '/' + dir + '/' + dir + '_LocList'

        in_file2 = open(taxon_LocList_file, 'r')

        mylist = in_file2.read().splitlines()

        contig_string = ""

        lines = mylist[i : i + int(contigs_dict[key])] #get a slice of the file

        contig_string = contig_string.join(lines)

        out_file.write(contig_string + '\n')


    i += int(contigs_dict[key])

in_file2.close()
out_file.close()
