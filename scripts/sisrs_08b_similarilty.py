#!/usr/bin/env python3

import csv
import sys
from os import system

outpath = sys.argv[1]
contigs = outpath + "/SISRS_Run/contigs_for_probes.fa"
matrix_file = outpath + "/SISRS_Run/contig_similarity.txt"

blast_command = f"blastn -query {contigs} -subject {contigs} -task blastn -outfmt \"6 qseqid sseqid qlen slen length pident nident\" > {matrix_file}"
#print(blast_command)
system(blast_command)

print('Query_ID Subject_ID Query_len Subject_len Overlap_length %ident #ident')

with open(matrix_file, newline = '') as file: #qseqid sseqid qlen slen length pident nident
    contigs = csv.reader(file, delimiter='\t') #iterator
    for line in contigs:
        if line[0] != line[1]: #not the same contig
            if float(line[5]) > 90: #90% identity for region that overlaps
                if float(line[4]) > 0.9 * min(float(line[2]), float(line[3])): #length of overlap is greater than 90% of the smaller contigi
                    print(line)

print('done comparison')
