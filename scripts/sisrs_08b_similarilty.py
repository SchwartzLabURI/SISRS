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
print(loci_to_drop)

#read fasta and remove loci
contig_seqs = []
with open(fasta) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        if record.id not in loci_to_drop:
            contig_seqs.append(record)
contig_seqs.sort(key=lambda r: -len(r))

if sys.argv[2] == 'ALL':
    SeqIO.write(contig_seqs, fasta + '.fasta', "fasta")
else:
    final_total = int(sys.argv[2])
    SeqIO.write(contig_seqs[0:final_total], fasta + '.fasta', "fasta")

print('done comparison')
