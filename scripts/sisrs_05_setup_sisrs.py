#!/usr/bin/env python3

'''

This script prepares data for a SISRS run by setting the data up and creating
mapping scripts
Contigs are renamed and moved to the SISRS_Run/Composite_Genome directory
The composite genome is indexed by Bowtie2 and Samtools
SISRS scripts are generated from a template and saved to the SISRS_Run/TAXA folder
'''

import os
from os import path
import sys
from glob import glob
from cmdCheck import *
import pandas as pd
from Bio import SeqIO
import argparse
import re

'''
This function is designed to obtain all of the directories needed to finish the
setup of sisrs. This function only needs the working sisrs directory as a
parameter.
'''
def obtainDir(outPath):
    ''' This function is designed to obtain all of the directories needed to finish the setup of SISRS. '''

    #Set TrimRead + SISRS directories based off of script folder location
    trim_read_dir =  outPath+"/Reads/TrimReads"

    #taxa come from taxonlist file
    with open(outPath+'/TaxonList.txt') as f:
        taxa = f.readlines()
    taxa = sorted([x.rstrip() for x in taxa]) 

    trim_read_tax_dirs = [trim_read_dir+'/'+x for x in taxa]
    ray_dir = outPath+"/Ray_Composite_Genome"
    sisrs_dir = outPath+"/SISRS_Run"

    trim_read_tax_dirs = [x for x in trim_read_tax_dirs if not x.endswith('trimOutput/')]
    trim_read_tax_dirs = [x for x in trim_read_tax_dirs if not x.endswith('fastqcOutput/')]
    trim_read_tax_dirs = sorted([x for x in trim_read_tax_dirs if not x.endswith('subsetOutput/')])

    #Create composite genome folder
    if(not path.isdir(sisrs_dir+"/Composite_Genome")):
        os.mkdir(sisrs_dir+"/Composite_Genome")
    composite_dir =sisrs_dir+"/Composite_Genome"

    return trim_read_tax_dirs,ray_dir,sisrs_dir,composite_dir

'''
This function is designed to rename the contig files and make copies of the
length files that are a product of the output from Ray. It requieres as input
the directory were ray output all of its dated and the composite directory.
'''
def fileChanges(ray_dir,composite_dir):
    ''' This function renames the contig files and make copies of the length files from Ray output. '''

    #Create renamed contig file in SISRS_Run/Composite_Genome
    rename_command = [
        'rename.sh',
        'in={}'.format(ray_dir+"/Contigs.fasta"),
        'out={}'.format(composite_dir+'/contigs.fa'),
        'prefix=SISRS',
        'addprefix=t',
        'ow=t']
    os.system(" ".join(rename_command))

    #Copy sequence length file from Ray
    contig_length = pd.read_csv(ray_dir+'/ContigLengths.txt',sep="\t",header=None)
    contig_length.columns = ['Contig','Length']
    contig_length.Contig = ("SISRS_" + contig_length.Contig)
    contig_length.to_csv(composite_dir+"/contigs_SeqLength.tsv",sep="\t",header=None,index=False)

'''
This function is designed to call the bowtie command and the samtools command
in order to index the composite genomes correctly. This function requiers two
parameters the composite genome directory and the number of threads that it can
run with
'''
def indexCompGenome(composite_dir,threads):
    ''' This function calls Bowtie and Samtools for indexing the composite genome. '''

    #Index composite genome using Bowtie2 and Samtools
    bowtie_command = [
        'bowtie2-build',
        '{}'.format(composite_dir+'/contigs.fa'),
        '{}'.format(composite_dir+'/contigs'),
        '-p',
        '{}'.format(str(threads))]
    os.system(" ".join(bowtie_command))

    samtools_command = [
        'samtools',
        'faidx',
        '{}'.format(composite_dir+'/contigs.fa')]
    os.system(" ".join(bowtie_command))

'''
This function is designed to start the file setups required for sisrs to run
properly. It requires as parameters the composite_dir
'''
def beginSetUp(composite_dir):
    ''' This function starts the file setups required for sisrs to run properly. '''

    #Create file with an entry for each individual site in the alignment
    siteCount=0
    locListFile = open(composite_dir+'/contigs_LocList','w') #originally it was an append mode 'a+'
    with open(composite_dir +"/contigs_SeqLength.tsv","r") as filein:
        for line in iter(filein):
            splitline=line.split()
            for x in range(1,(int(splitline[1])+1)):
                locListFile.write((splitline[0] +'/'+str(x)+'\n'))
                siteCount+=1
    locListFile.close()

    print("==== Site list created: " + str(siteCount) + " total sites ==== \n",flush=True)

if __name__ == '__main__':

    # Get arguments
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-d','--directory',action='store',nargs="?")
    my_parser.add_argument('-p','--processors',action='store',default=1,nargs="?")
    args = my_parser.parse_args()

    sis = args.directory
    proc = args.processors

    trim_read_tax_dirs,ray_dir,sisrs,composite_dir = obtainDir(sis)
    fileChanges(ray_dir,composite_dir)
    indexCompGenome(composite_dir,proc)
    beginSetUp(composite_dir)
