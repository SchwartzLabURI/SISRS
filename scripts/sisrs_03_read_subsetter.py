#!/usr/bin/env python3

'''

IMPORTANT: Trimmed reads should be QC checked prior to this step to remove
datasets with low quality or high deplication rates
This script performs read subsetting for a SISRS composite genome assembly
This script calls reformat.sh  (part of BBMap Suite)

All reads for all taxa should be in .fastq.gz format (To change this, find/replace
this script, replacing '.fastq.gz' with your chosen extension)
Paired-end read files must be identically basenamed and end in _Trim_1.fastq.gz/_Trim_2.fastq.gz
Input: (1) -gs/--genomesize Genome size estimate in basepairs (e.g. Human: 3500000000)
Output: For an analysis with X taxa and a genome size estimate of Y bp,  this script
will subset each taxons reads down to (10*Y)/X bases (~10X total coverage when summed)


'''

import os
from os import path
import sys
from glob import glob
import subprocess
from cmdCheck import *
from subprocess import check_call
import pandas as pd
import re
import argparse

'''
This function is designed to do the remaining folder setup for the read subsetting.
It also finishes up the other minor setups needed for this script. The arguments
that are needed for this script are the working sisrs directory and the genomeSize
estimation and output path
'''
def setupDir(sisrs_dir,genomeSize):
    ''' This function does the remaining folder setup for the read subsetting. '''

    # returned list of items
    # trim_read_dir         --> 0
    # subset_output_dir     --> 1
    # subset_log_dir        --> 2
    # trim_read_tax_dirs    --> 3
    # subsetDepth           --> 4
    # df                    --> 5
    # compiled_paired       --> 6
    # compiled_single_end   --> 7
    rtn = []

    #Set TrimRead directories based off of script folder location
    trim_read_dir = sisrs_dir + "/Reads/TrimReads"
    rtn += [trim_read_dir]

    #taxa come from taxonlist file
    with open(sisrs_dir+'/TaxonList.txt') as f:
        taxa = f.readlines()
    taxa = [x.rstrip() for x in taxa]
    taxa = sorted(taxa)
    
    #Find taxa folders within TrimRead folder
    trim_read_tax_dirs = [trim_read_dir+"/"+x for x in taxa]
    print(trim_read_tax_dirs)

    #Create folder for Subset reads
    if(not path.isdir(sisrs_dir +"/Reads/SubsetReads")):
        os.mkdir(sisrs_dir+"/Reads/SubsetReads")
    subset_output_dir = sisrs_dir+"/Reads/SubsetReads/"
    rtn += [subset_output_dir]

    #Create folder for Subset output logs
    if(not path.isdir(trim_read_dir+"/subsetOutput")):
        os.mkdir(trim_read_dir+"/subsetOutput")
    subset_log_dir = trim_read_dir+"/subsetOutput/"
    rtn += [subset_log_dir]

    #Ensure Trim Log/FastQC/subset Log output directories not included in taxa list
    trim_read_tax_dirs = [x for x in trim_read_tax_dirs if not x.endswith('trimOutput/')]
    trim_read_tax_dirs = [x for x in trim_read_tax_dirs if not x.endswith('fastqcOutput/')]
    trim_read_tax_dirs = [x for x in trim_read_tax_dirs if not x.endswith('subsetOutput/')]
    rtn += [trim_read_tax_dirs]

    #Calculate subset depth
    rtn += [int((10*genomeSize)/len(trim_read_tax_dirs))]

    #Initialize Pandas DF to get base counts and lists of paired and single-end files
    rtn += [pd.DataFrame(columns=['Taxon','Dataset','Basecount'])]
    rtn += [list()]
    rtn += [list()]

    return rtn

'''
This function is designed to break up the load that is carried by countBasePair.
'''
def countHelper(tax_dir,compiled_paired,compiled_single_end):
    ''' This is helper function for countBasePair function to break up the load. '''

    #List all files and set output dir
    files = sorted(glob(tax_dir+"/*.fastq.gz"))
    taxon_ID = path.basename(tax_dir)

    left_pairs = list()
    right_pairs = list()
    single_end = list()

    taxon_list = list()

    dataset_list = list()
    basecount_list = list()

    #Find files ending in _1/_2.fastq.gz
    left_files = [s for s in files if "_Trim_1.fastq.gz" in s]
    right_files = [s for s in files if "_Trim_2.fastq.gz" in s]

    #Strip _1_Trim.fastq.gz/_2_Trim.fastq.gz and identify pairs based on file name
    left_files = [x.replace('_Trim_1.fastq.gz', '') for x in left_files]
    right_files = [x.replace('_Trim_2.fastq.gz', '') for x in right_files]
    paired_files = list(set(left_files).intersection(right_files))

    #Reset file names and filter out single-end files
    for pair in paired_files:
        left_pairs.append(pair+"_Trim_1.fastq.gz")
        right_pairs.append(pair+"_Trim_2.fastq.gz")
    paired_files = sorted(left_pairs + right_pairs)

    single_end = [x for x in files if x not in paired_files]

    compiled_paired = compiled_paired + paired_files
    compiled_single_end = compiled_single_end + single_end
    #Remove .fastq.gz from lists to make naming easier
    left_pairs = [x.replace('_1.fastq.gz', '') for x in left_pairs]
    right_pairs = [x.replace('_2.fastq.gz', '') for x in right_pairs]
    single_end = [x.replace('.fastq.gz', '') for x in single_end]

    # This return statement looks this way because of an error that is assosicated
    # with turning it into a list. I will attempt to fix it again later
    return compiled_paired,compiled_single_end,left_pairs,right_pairs,single_end,taxon_list,dataset_list,basecount_list,taxon_ID

'''
Function to execute the first for loop using helper function
'''
def countBasePair(trim_read_tax_dirs,compiled_paired,compiled_single_end,df):
    ''' This function counts bases in single-end files and paired-end files if present. '''

    #For each taxa directory...
    for tax_dir in trim_read_tax_dirs:

        outPut = []
        outPut += countHelper(tax_dir,compiled_paired,compiled_single_end)
        compiled_paired = outPut[0]
        compiled_single_end = outPut[1]
        left_pairs = outPut[2]
        right_pairs = outPut[3]
        single_end = outPut[4]
        taxon_list = outPut[5]
        dataset_list = outPut[6]
        basecount_list = outPut[7]
        taxon_ID = outPut[8]

        #Count bases in single-end files if present...
        if len(single_end) > 0:
            for x in single_end:
                count_command = [
                    'reformat.sh',
                    'in={}'.format(x+'.fastq.gz')]
                proc = subprocess.Popen(count_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                o, e = proc.communicate()

                #Munge output to get base count
                count = re.search('HERE (.*) HERE_2', str(e).replace('(100.00%) \\t','HERE ').replace(' bases (100.00%)\\n\\nTime:',' HERE_2'))
                baseCount = count.group(1)
                taxon_list.append(taxon_ID)
                dataset_list.append(path.basename(x))
                basecount_list.append(baseCount)

        #Count bases in paired-end files if present...
        if(len(left_pairs) == len(right_pairs) & len(left_pairs) > 0):
            for x in range(len(left_pairs)):
                dataset_list.append(path.basename(left_pairs[x]))
                count_command = [
                    'reformat.sh',
                    'in={}'.format(left_pairs[x]+'_1.fastq.gz'),
                    'in2={}'.format(right_pairs[x]+'_2.fastq.gz')]
                proc = subprocess.Popen(count_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                o, e = proc.communicate()

                count = re.search('HERE (.*) HERE_2', str(e).replace('(100.00%) \\t','HERE ').replace(' bases (100.00%)\\n\\nTime:',' HERE_2'))
                baseCount = count.group(1)
                basecount_list.append(int(baseCount))
                taxon_list.append(taxon_ID)
        list_of_tuples = list(zip(taxon_list,dataset_list, basecount_list))

        df = df.append((pd.DataFrame(list_of_tuples, columns = ['Taxon','Dataset', 'Basecount'])))

    return df, compiled_paired, compiled_single_end

'''
Function to execute the second for loop
'''
def checkCoverage(df,subsetDepth,subset_output_dir,compiled_paired,compiled_single_end):
    ''' This function checks taxa total coverage. '''

    df = df.reset_index(drop=True)
    df["Basecount"] = pd.to_numeric(df["Basecount"])
    df = df.sort_values(by=['Taxon'])
    df["ToSample"] = 0
    taxa_list = sorted(df.Taxon.unique())

    for taxa in taxa_list:
        taxaDF = df[df['Taxon']==taxa]

        #Check if taxa has enough total coverage. If not, use all reads and warn user.
        if(taxaDF["Basecount"].sum() < subsetDepth):
            print("NOTE: " + taxa + " samples only have a total of " + str(taxaDF["Basecount"].sum()) + " bp. All " + taxa + " reads will be used.\nWARNING: Having under the required read threshold will have negative consequences on composite genome assembly and site calling.\n")
            df.loc[df["Taxon"] == taxa, ["ToSample"]] = list(df[df["Taxon"] == taxa].Basecount)

        elif(taxaDF["Basecount"].sum() >= subsetDepth):
            #If there is enough total data, check sample values to optimize sampling
            dataset_count = taxaDF.shape[0] #Count samples for taxa
            per_dataset = subsetDepth/dataset_count #Calculate bases needed per sample if evenly distributed
            s = pd.Series(taxaDF["Basecount"]) #Grab basecount colulmn

            if((s < per_dataset).any()):
                subset_countdown =  subsetDepth #Create countdown variable

                while((s < per_dataset).any()):

                    low_datasets = list(taxaDF[taxaDF.Basecount<per_dataset]["Dataset"]) #Find datasets with too few bases
                    high_datasets = list(taxaDF[taxaDF.Basecount>=per_dataset]["Dataset"]) #Find datasets with enough bases (currently)

                    sum_low = taxaDF[taxaDF.Basecount<per_dataset].Basecount.sum() #Sum the bases of the low datasets

                    print("Not all " + taxa + " samples have enough for even read sampling (" + ",".join(low_datasets) + "). Will take all reads from data poor samples and compensate from larger samples...")
                    df.loc[df.Dataset.isin(low_datasets), ["ToSample"]] = list(df[df.Dataset.isin(low_datasets)].Basecount) #Set low datasets to use all reads

                    taxaDF = taxaDF[taxaDF.Dataset.isin(high_datasets)] #Reset taxa DF to remove low samples
                    dataset_count = taxaDF.shape[0] #Count remaining datasets
                    subset_countdown = subset_countdown - sum_low

                    per_dataset = int(subset_countdown/dataset_count) #Reset per_dataset

                    s = pd.Series(taxaDF.Basecount) #Recount basecounts for while loop

                df.loc[df.Dataset.isin(high_datasets), ["ToSample"]] = per_dataset
            else:
                df.loc[df.Taxon == taxa, ["ToSample"]] = per_dataset

    df.to_csv(subset_output_dir+"/Subset_Scheme.csv",index=False)

    compiled_paired = [path.basename(x).replace('_1.fastq.gz', '') for x in compiled_paired]
    compiled_paired = [path.basename(x).replace('_2.fastq.gz', '') for x in compiled_paired]
    compiled_single_end = [path.basename(x).replace('.fastq.gz', '') for x in compiled_single_end]

    out = []
    # Return the three necessary pieces of information
    out += [compiled_paired]
    out += [compiled_single_end]
    out += [df]

    return out

'''
This function is designed to help strip the list of files that are going to be
subsetted. They will be stripped of the file path, file extension, and if there
is paired files the _1 and _2 will also be stripped
'''
def stripFiles(aList):
    ''' This function strip the list of files that are going to be subsetted. '''

    nList = []
    for item in aList:
        fName = os.path.basename(item)
        save = fName.strip('.fastq.gz')
        if '_1' in save:
            nList += [save.strip('_1')]
        elif '_2' in save:
            nList += [save.strip('_2')]
        else:
            nList += [save]
    return nList

'''
Function to do all of the subsetting the is needed by sisrs. This is a general
function and can handle single and paired ends genomes.
'''
def subset(compiledList,df,trim_read_dir,subset_output_dir,subset_log_dir,paired):
    ''' This function performs file subsetting (single and paired-end). '''

    cList = stripFiles(compiledList)

    if len(cList) > 0:

        end_DF = df[df.Dataset.isin(cList)]

        for row in end_DF.itertuples():

            subset_command = []
            if paired:

                subset_command = [
                    'reformat.sh',
                    'in={}'.format(trim_read_dir+"/"+str(row.Taxon)+"/"+str(row.Dataset)+"_1.fastq.gz"),
                    'in2={}'.format(trim_read_dir+"/"+str(row.Taxon)+"/"+str(row.Dataset)+"_2.fastq.gz"),
                    'out={}'.format(subset_output_dir+str(row.Dataset).replace('Trim','GenomeReads_1.fastq.gz')),
                    'out2={}'.format(subset_output_dir+str(row.Dataset).replace('Trim','GenomeReads_2.fastq.gz')),
                    'samplebasestarget={}'.format(int(row.ToSample)),
                    'ow=t',
                    '&>',
                    '{outDir}out_{fileName}'.format(outDir=subset_log_dir,fileName=row.Dataset)]
            else:
                subset_command = [
                    'reformat.sh',
                    'in={}'.format(trim_read_dir+"/"+str(row.Taxon)+"/"+str(row.Dataset)+".fastq.gz"),
                    'out={}'.format(subset_output_dir+str(row.Dataset).replace('Trim','GenomeReads.fastq.gz')),
                    'samplebasestarget={}'.format(int(row.ToSample)),
                    'ow=t',
                    '&>',
                    '{outDir}out_{fileName}'.format(outDir=subset_log_dir,fileName=row.Dataset)]
            check_call(subset_command)

if __name__ == '__main__':

    # Get arguments
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-dir','--directory',action='store',nargs="?")
    my_parser.add_argument('-gs','--genome_size',type=int, action='store',nargs="?")
    args = my_parser.parse_args()

    sis = args.directory
    gs = args.genome_size

    setupInfo = setupDir(sis,gs)

    df, compiled_paired, compiled_single_end = countBasePair(setupInfo[3],setupInfo[6], setupInfo[7],setupInfo[5])

    print("Based on a genome size estimate of " + str(gs) + " bp, and with " + str(len(setupInfo[3])) + " species, the requested subset depth is " + str(setupInfo[4]) + " bp per species")

    out = checkCoverage(df,setupInfo[4],setupInfo[1],compiled_paired, compiled_single_end)

    # Subset Single End
    subset(compiled_single_end,out[2],setupInfo[0],setupInfo[1],setupInfo[2],False)

    #Subset paired ends
    subset(compiled_paired,out[2],setupInfo[0],setupInfo[1],setupInfo[2],True)
