#!/usr/bin/env python3

#!/usr/bin/env python3

import os
from os import path
import sys
from glob import glob
import pandas as pd

from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment,AlignInfo
from Bio.Align.Applications import MafftCommandline
from Bio.Alphabet import generic_dna,IUPAC, Gapped
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

script_dir = sys.path[0]
processors = str(int(sys.argv[1]))
ref_species = sys.argv[2]
contig_count = str(int(sys.argv[3]))

site_id = 'SISRS_Biallelic_NoMissing'

post_processing_dir = path.dirname(path.abspath(script_dir))+"/Post_SISRS_Processing"
sisrs_site_dir = path.dirname(path.abspath(script_dir))+"/Post_SISRS_Processing/SISRS_Sites"

if(not path.isdir(post_processing_dir+"/Chronos_"+contig_count)):
    os.mkdir(post_processing_dir+"/Chronos_"+contig_count)
chronos_dir = post_processing_dir+"/Chronos_"+contig_count

sisrs_dir = path.dirname(path.abspath(script_dir))+"/SISRS_Run"
sisrs_tax_dirs = sorted(glob(sisrs_dir+"/*/"))
sisrs_tax_dirs = [x for x in sisrs_tax_dirs if not x.endswith('Composite_Genome/')]

all_sisrs=  pd.read_csv(sisrs_site_dir+'/'+ref_species+'_'+site_id+'_NoGaps_LocList.txt',header=None)[0].str.split("/", n = 1, expand = True)
all_sisrs.columns=['Contig','Site']

top_sisrs = pd.DataFrame(all_sisrs.groupby(['Contig']).size()).sort_values(0,ascending=False).head(int(contig_count)).index.tolist()

with open(chronos_dir+'/SISRS_'+contig_count+'_ContigList.txt','w') as outfile:
    for contig in top_sisrs:
        outfile.write(contig+"\n")

for tax_dir in sisrs_tax_dirs:
    taxa = path.basename(tax_dir[:-1])
    if(not path.isdir(chronos_dir+'/'+taxa)):
        os.mkdir(chronos_dir+'/'+taxa)
    chronos_tax_dir = chronos_dir+'/'+taxa
    sub_command = ['filterbyname.sh',
                   'in={}'.format(tax_dir+'contigs.fa'),
                   'out={}'.format(chronos_tax_dir+'/ref_genes.fa'),
                   'names={}'.format(chronos_dir+'/SISRS_'+contig_count+'_ContigList.txt'),
                   'include=t',
                   'ow=t']
    os.system(' '.join(sub_command))

fafiles = sorted(list(glob(chronos_dir+'/*/ref_genes.fa')))

if(not path.isdir(chronos_dir+'/unaligned')):
    os.mkdir(chronos_dir+'/unaligned')
un_outfolder = chronos_dir+'/unaligned'

if(not path.isdir(chronos_dir+'/aligned')):
    os.mkdir(chronos_dir+'/aligned')
al_outfolder = chronos_dir+'/aligned'

speciesList = list()
for f in fafiles:
    speciesList.append(os.path.basename(os.path.dirname(f)))

#Add species FASTAs to list, changing record IDs to species ID
fastaList = []

#Append all renamed records to new list per species
for i in range(0,len(fafiles)):
    spFasta = SeqIO.parse(fafiles[i],"fasta")
    tempList = []
    for record in spFasta:
        record.id = speciesList[i]
        record.name = ''
        record.description = ''
        tempList.append(record)
    fastaList.append(tempList)

#Create multifasta for each individual contig
alignList = []
for i in range(0,len(top_sisrs)):
    tempList = []
    for sp in fastaList:
        tempList.append(sp[i])
    alignList.append(tempList)

for i in range(0,len(top_sisrs)):
    SeqIO.write(alignList[i], un_outfolder+'/'+top_sisrs[i]+'_Unaligned.fa', "fasta")

#Generate MAFFT command list
unfiles = sorted(list(glob(un_outfolder+'/*Unaligned.fa')))
alfiles = []
for un in unfiles:
    alfiles.append(str.replace(un,'Unaligned','Aligned').replace('unaligned','aligned').replace('.fa','.phylip-relaxed'))

with open(un_outfolder+'/Mafft_Commands',"w") as outfile:
    for i in range(0,len(unfiles)):
        outfile.write('mafft --auto '+unfiles[i]+' | python '+script_dir+'/pipe_fasta_to_phyliprelaxed.py ' + alfiles[i]+"\n")

mafft_command = [
    'parallel',
    '--jobs',
    '{}'.format(processors),
    '<',
    '{}'.format(un_outfolder+'/Mafft_Commands')]
os.system(' '.join(mafft_command))

#Generate partition file while concatenating alignments
partitionList = []
i=0
t=1

cat_data = {s:list() for s in speciesList}

for f in alfiles:
    data = SeqIO.to_dict(SeqIO.parse(f.rstrip(), "phylip-relaxed"))
    alignLength = len(data[list(data.keys())[0]].seq)
    partitionList.append('DNA, p'+str(i+1)+'='+str(t)+'-'+str(((t-1)+alignLength)))
    i+=1
    t+=alignLength
    for k in speciesList:
        if k in data:
            cat_data[k].extend(list(data[k].seq))
        else:
            cat_data[k].extend(['N'] * len(data[list(data.keys())[0]].seq))

with open(chronos_dir+"/Concat_"+contig_count+"_Partition", "w") as outfile:
    for s in partitionList:
        outfile.write("%s\n" % s)

datalist=[]
for k,v in cat_data.items():
    seq = SeqRecord(Seq(''.join(v)), id=k)
    datalist.append(seq)

SeqIO.write(datalist,chronos_dir+"/Concat_"+contig_count+".phylip-relaxed", "phylip-relaxed")
