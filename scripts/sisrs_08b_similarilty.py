#!/usr/bin/env python3

import csv
import sys
from Bio import SeqIO
from os import system

outpath = sys.argv[1]
fasta = outpath + "/SISRS_Run/contigs_for_probes.fa"
matrix_file = outpath + "/SISRS_Run/contig_similarity.txt"

blast_command = f"blastn -query {fasta} -subject {fasta} -task blastn -outfmt \"6 qseqid sseqid qlen slen length pident nident\" > {matrix_file}"
#print(blast_command)
system(blast_command)

loci_to_drop = []

print('Query_ID Subject_ID Query_len Subject_len Overlap_length %ident #ident')

with open(matrix_file, newline = '') as file: #qseqid sseqid qlen slen length pident nident
    contigs = csv.reader(file, delimiter='\t') #iterator
    for line in contigs:
        if line[0] != line[1]: #not the same contig
            if float(line[5]) > 90: #90% identity for region that overlaps
                if float(line[4]) > 0.9 * min(float(line[2]), float(line[3])): #length of overlap is greater than 90% of the smaller contigi
                    print(line)
                    if int(line[2]) > int(line[3]):
                        loci_to_drop.append(line[1])
                    else:
                        loci_to_drop.append(line[0])

loci_to_drop = list(set(loci_to_drop))
print(loci_to_drop)

#read fasta and remove loci
contig_seqs = []
with open(fasta) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        if record.id not in loci_to_drop:
            contig_seqs.append(record)
contig_seqs.sort(key=lambda r: -len(r))
SeqIO.write(contig_seqs, fasta + '.fasta', "fasta")

print('done comparison')
