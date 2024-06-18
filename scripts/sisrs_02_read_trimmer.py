#!/usr/bin/env python3

'''
Trim reads using BBDuk

Paired-end read files must be identically basenamed and end in _1/_2
Input: (Optional) -p/--processors Number of available processors for FastQC
(Default: 1, For 10 processors: python sisrs_02_read_trimmer.py -p 10)
# Output: (1) Trimmed Reads (in <base_dir>/Reads/TrimReads/<Taxon>/<Read_Basename>_Trim.fastq.gz
# Output: (2) Trim Log (in <base_dir>/Reads/RawReads/trimOutput)
# Output: (3) FastQC output for all raw + trimmed read sets (in <base_dir>/Reads/RawReads/fastqcOutput & <base_dir>/Reads/TrimReads/fastqcOutput)
'''

import os
from os import path, listdir
import sys
from glob import glob
import subprocess
from subprocess import check_call
import argparse
from sisrs_08_filter_contigs import get_taxon_list

def findAdapter():
    '''
    This function is designed to find BBDuk + Adapter File.

    Arguments: none.

    Returns:  
    string: the path to the location of the adapter file.
    '''

    cmd = ['which', 'bbduk.sh']
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    o, e = proc.communicate()
    bbduk_adapter = path.dirname(o.decode('ascii')) + "/resources/adapters.fa"

    return bbduk_adapter


def setup(datadir, sisrs_dir, taxa):
    '''
    This function does the remaining folder setup that is needed to do
    trim reads.

    Arguments: 
    datadir (string): folder containing folders containing reads
    sisrs_dir (string): output directory
    taxa (list): taxon list

    Returns: 
    raw read and trimmed data directories.
    '''

    #Set RawRead and TrimRead directories based off of script folder location
    raw_read_dir = datadir
    trim_read_dir = sisrs_dir + "/Reads/TrimReads"

    #Find taxa folders within RawRead folder
    raw_read_tax_dirs = [raw_read_dir + '/' + x for x in taxa]

    #Create folder for BBDuk StdOut
    if not path.isdir(sisrs_dir + "/Reads/RawReads/trimOutput"):
        os.mkdir(sisrs_dir + "/Reads/RawReads/trimOutput")
    trim_output = sisrs_dir + "/Reads/RawReads/trimOutput"

    #Create folder for Raw FastQC output
    if not path.isdir(sisrs_dir + "/Reads/RawReads/fastqcOutput"):
        os.mkdir(sisrs_dir + "/Reads/RawReads/fastqcOutput")
    raw_fastqc_output = sisrs_dir + "/Reads/RawReads/fastqcOutput"

    #Create folder for Trimmed FastQC output
    if not path.isdir(trim_read_dir + "/fastqcOutput"):
        os.mkdir(trim_read_dir + "/fastqcOutput")
    trim_fastqc_output = trim_read_dir + "/fastqcOutput/"

    return raw_read_dir, trim_read_dir, trim_output, raw_fastqc_output, trim_fastqc_output, raw_read_tax_dirs

def newdFastqc(processors, fastqc_output, data_dir, newFiles):
    '''
    This function is designed to run the FastQC command with new data only. Modified
    to run both the raw and the trimmed data, using all available processors

    Arguments: the number of processors, the output directory, and the read directory.

    Returns: none.
    '''

    #Run FastQC on all trimmed files, using all available processors
    fastqc_command = [
        'fastqc',
        '-t',
        '{}'.format(processors),
        '-o',
        '{}'.format(fastqc_output)]

    for item in data_dir:
        for x in glob(item + "/*.fastq.gz"):
            f = x.split('/')

            if f[-1] in newFiles:
                fastqc_command.append(x)

    check_call(fastqc_command)

def fastqcCommand(processors, fastqc_output, read_dir, taxa):
    '''
    Run FastQC on fastq.gz in specified taxon folders, using specified number of processors.

    Arguments: 
    processors (int): the number of processors
    fastqc_output (string): path to fastqc output directory
    read_dir (string): path to raw reads directory
    taxa (list): taxon list

    Returns: none.
    '''

    fastqc_command = [
        'fastqc',
        '-t',
        '{}'.format(processors),
        '-o',
        '{}'.format(fastqc_output)]

    for t in taxa:
        print(read_dir + "/" + t + "/*.fastq.gz")
        flist = glob(read_dir + "/" + t + "/*.fastq.gz")
        if len(flist) == 0:
            sys.exit("No fastq.gz files found - please check your paths and extensions")
        print(flist)
        for x in flist:  #can add recursive here
            print(x)
            fastqc_command.append(x)

    check_call(fastqc_command)

def trimHelper(tax_dir, trim_read_dir, newData):
    '''
    This function obtains all of the
    single end and paired end read files in a folder.

    Arguments: 
    tax_dir (string): path to taxpm folder with fastq.gz files
    trim_read_dir (string): path to trimmed read directory
    newData (list)

    Returns: 
    string: path to trimmed directory for the specific taxon 
    list: left paired end read paths (minus _1.fastq.gz)
    list: right paired end read paths (minus _2.fastq.gz) 
    list: single end read paths (minus .fastq.gz)
    '''

    #List all files and set output dir
    files = sorted(glob(tax_dir + "/*.fastq.gz"))
    taxon_ID = path.basename(tax_dir)
    out_trim_dir = trim_read_dir + "/" + taxon_ID

    left_pairs = list()
    right_pairs = list()
    single_end = list()

    #Find files ending in _1/_2.fastq.gz
    left_files = [s for s in files if "_1.fastq.gz" in s]
    right_files = [s for s in files if "_2.fastq.gz" in s]

    lf = []
    rf = []
    if newData != []:
        for f in left_files:
            o = f.split('/')
            if o[-1] in newData:
                lf += [f]
        left_files = lf
        for f in right_files:
            o = f.split('/')
            if o[-1] in newData:
                rf += [f]
        right_files = rf

    #Strip _1.fastq.gz/_2.fastq.gz and identify pairs based on file name
    left_files = [x.replace('_1.fastq.gz', '') for x in left_files]
    right_files = [x.replace('_2.fastq.gz', '') for x in right_files]

    paired_files = list(set(left_files).intersection(right_files))

    #Reset file names and filter out single-end files
    for pair in paired_files:
        left_pairs.append(pair + "_1.fastq.gz")
        right_pairs.append(pair + "_2.fastq.gz")
    paired_files = sorted(left_pairs + right_pairs)

    single_end = [x for x in files if x not in paired_files]

    se = []
    if newData != []:
        for f in single_end:
            o = f.split('/')
            if o[-1] in newData:
                se += [f]
        single_end = se

    #Remove .fastq.gz from lists to make naming easier
    left_pairs = [x.replace('_1.fastq.gz', '') for x in left_pairs]
    right_pairs = [x.replace('_2.fastq.gz', '') for x in right_pairs]
    single_end = [x.replace('.fastq.gz', '') for x in single_end]

    return out_trim_dir, left_pairs, right_pairs, single_end

def trim(raw_read_tax_dirs, trim_read_dir, bbduk_adapter, trim_output, newData):
    '''
    Trim fastq files provided to
    the program and write to new folder. 

    Arguments: 
    raw_read_tax_dirs (list): paths to raw read directory
    trim_read_dir (string): path to trimmed read directory for output
    bbduk_adapter (string): path to file
    trim_output (string): path to trim output directory
    new_Data (list):

    Returns: none
    '''

    #For each taxa directory...
    for tax_dir in raw_read_tax_dirs:
        out_trim_dir, left_pairs, right_pairs, single_end = trimHelper(tax_dir, trim_read_dir, newData)

        #Trim single-end files if present...
        if len(single_end) > 0:
            for x in single_end:
                se_trim_command = [
                    'bbduk.sh',
                    'maxns=0',
                    'ref={}'.format(bbduk_adapter),
                    'qtrim=w',
                    'trimq=15',
                    'minlength=35',
                    'maq=25',
                    'in={}'.format(x + '.fastq.gz'),
                    'out={}'.format(out_trim_dir + "/" + path.basename(x) + '_Trim.fastq.gz'),
                    'k=23',
                    'mink=11',
                    'hdist=1',
                    'hdist2=0',
                    'ktrim=r',
                    'ow=t',
                    '&>',
                    '{outDir}out_{fileName}_Trim'.format(outDir=trim_output, fileName=path.basename(x))]
                check_call(se_trim_command)

        #Trim paired-end files if present...
        if (len(left_pairs) == len(right_pairs) & len(left_pairs) > 0):
            for x in range(len(left_pairs)):
                file_name = path.basename(left_pairs[x])
                pe_trim_command = [
                    'bbduk.sh',
                    'maxns=0',
                    'ref={}'.format(bbduk_adapter),
                    'qtrim=w',
                    'trimq=15',
                    'minlength=35',
                    'maq=25',
                    'in={}'.format(left_pairs[x] + '_1.fastq.gz'),
                    'in2={}'.format(right_pairs[x] + '_2.fastq.gz'),
                    'out={}'.format(out_trim_dir + "/" + path.basename(left_pairs[x]) + '_Trim_1.fastq.gz'),
                    'out2={}'.format(out_trim_dir + "/" + path.basename(right_pairs[x]) + '_Trim_2.fastq.gz'),
                    'k=23',
                    'mink=11',
                    'hdist=1',
                    'hdist2=0',
                    'ktrim=r',
                    'ow=t',
                    '&>',
                    '{outDir}out_{fileName}_Trim'.format(outDir=trim_output, fileName=file_name)]
                check_call(pe_trim_command)

def makeLinks(data_path, sisrs_dir, taxa_list):

    '''
    This function was designed to make soft links to all the files to avoid repeated trimming

    Arguments:
    data_path (string): the path to the folder with folders of trimmed reads
    sisrs_dir (string): output dir
    taxa_list (list): a list of the taxon/folder names

    Returns: none.
    '''

    for x in taxa_list:
        # Only the files that do not start with '.'
        l = [f for f in listdir(data_path + '/' + x) if not f.startswith('.')]
        for i in l:
            # Creates the soft link to the files
            os.link(data_path + '/' + x + '/' + i,
                    sisrs_dir+"/Reads/TrimReads/" + x + '/' + i)

def run2(datadir, outputdir, processors, link):

    tl = get_taxon_list(outputdir)  # get taxa from TaxonList.txt

    if link:
        makeLinks(datadir, outputdir, tl)
    else:
        bba = findAdapter()  # returns path to adapter file in bbduk
        raw_read_dir, trim_read_dir, trim_output, raw_fastqc_output, trim_fastqc_output, raw_read_tax_dirs = setup(datadir, outputdir, tl)  # list of the folders

        # fastqc raw data
        fastqcCommand(processors, raw_fastqc_output, datadir, tl)  # raw_fastqc_output, raw_read_dir=datadir

        print(raw_read_tax_dirs)
        print(trim_read_dir)
        print(trim_output)
        trim(raw_read_tax_dirs, trim_read_dir, bba, trim_output, [])  # raw_read_tax_dirs, trim_read_dir, path to adapters, trim_output

        # fastqc trimmed data
        fastqcCommand(processors, trim_fastqc_output, trim_read_dir, tl)  # trim_fastqc_output, trim_read_dir folders

if __name__ == "__main__":
    # Get arguments
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-d', '--datadir', action='store', nargs="?")
    my_parser.add_argument('-dir', '--outputdir', action='store', nargs="?")
    my_parser.add_argument('-p', '--processors', action='store', default=1, nargs="?")
    my_parser.add_argument('-l', '--link', action=argparse.BooleanOptionalAction, default=False)
    args = my_parser.parse_args()

    run2(args.datadir, args.outputdir, args.processors, args.link)