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
import pandas as pd
import argparse
from sisrs_08_filter_contigs import get_taxon_list

def obtainDir(outPath):
    '''
    This function is designed to obtain all of the directories needed to finish the
    setup of SISRS. This function only needs the working sisrs directory as a
    parameter.

    Arguments: 
    outPath (string): path to the output directory.

    Returns: 
    list: paths to the trimmed reads directories
    string: path to ray composite directory
    string: path to SISRS Run folder 
    string: path to the composite genome directory.
    '''

    #Set TrimRead + SISRS directories based off of script folder location

    if path.isdir(outPath + "/Reads/TrimReads_nocl"):
        trim_read_dir = outPath + "/Reads/TrimReads_nocl"
    else:
        trim_read_dir = outPath + "/Reads/TrimReads"

    #taxa come from taxonlist file
    taxa = get_taxon_list(outPath+'/TaxonList.txt')

    trim_read_tax_dirs = [trim_read_dir+'/'+x for x in taxa]
    ray_dir = outPath+"/Ray_Composite_Genome"
    sisrs_dir = outPath+"/SISRS_Run"

    #Create composite genome folder
    if(not path.isdir(sisrs_dir+"/Composite_Genome")):
        os.mkdir(sisrs_dir+"/Composite_Genome")
    composite_dir = sisrs_dir+"/Composite_Genome"

    return trim_read_tax_dirs,ray_dir,sisrs_dir,composite_dir


def fileChanges(ray_dir,composite_dir):
    '''
    Renames and copies the Contigs file and converts the Ray contig length file to csv 

    Arguments: 
    ray_dir (string): path to Ray directory
    composite_dir (string): path to the composite genome directory.

    Returns: None.
    '''

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


def indexCompGenome(composite_dir,threads):
    '''
    Call bowtie2 build on contigs and samtools to index 

    Arguments: 
    composite_dir (string): path to the composite genome directory
    threads (int? string?): number of threads to be used.

    Returns: None.
    '''

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


def beginSetUp(composite_dir):
    '''
    Create file with an entry for each individual site in the alignment
    
    Arguments: 
    composite_dir (string): path to the composite genome directory.

    Returns: None.
    '''

    siteCount=0
    locListFile = open(composite_dir+'/contigs_LocList','w') #originally it was an append mode 'a+'
    with open(composite_dir +"/contigs_SeqLength.tsv","r") as filein:
        for line in iter(filein):
            splitline=line.split()
            for x in range(1, (int(splitline[1])+1)):
                locListFile.write((splitline[0] + '/' + str(x) + '\n'))
                siteCount += 1
    locListFile.close()

    print("==== Site list created: " + str(siteCount) + " total sites ==== \n",flush=True)

def makelinks(link, composite_dir):
    old_composite_dir = link + "/SISRS_Run/Composite_Genome"
    for f in os.listdir(old_composite_dir):
        os.link(old_composite_dir + '/' + f, composite_dir + '/' + f)

def run5(sis, proc, link):
    trim_read_tax_dirs, ray_dir, sisrs, composite_dir = obtainDir(sis)

    if link:
        makelinks(link, composite_dir)
    else:
        fileChanges(ray_dir, composite_dir)
        indexCompGenome(composite_dir, proc)
        beginSetUp(composite_dir)

if __name__ == '__main__':

    # Get arguments
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-d', '--directory', action='store', nargs="?")
    my_parser.add_argument('-p', '--processors', action='store', default=1, nargs="?")
    my_parser.add_argument('--link', action='store', nargs="?", default=None)
    args = my_parser.parse_args()

    sis = args.directory
    proc = args.processors

    run5(sis, proc, args.link)