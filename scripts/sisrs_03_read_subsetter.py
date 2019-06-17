#!/usr/bin/env python3

# TRIMMED READS SHOULD BE QC CHECKED PRIOR TO THIS STEP TO REMOVE DATASETS WITH LOW QUALITY OR HIGH DUPLICATION RATES
# This script performs read subsetting for a SISRS composite genome assembly
# This script calls reformat.sh  (part of BBMap Suite), which must be installed and in your path
# All reads for all taxa should be in .fastq.gz format (To change this, find/replace this script, replacing '.fastq.gz' with your chosen extension)
# Paired-end read files must be identically basenamed and end in _Trim_1.fastq.gz/_Trim_2.fastq.gz
# Input: (1) Genome size estimate in basepairs (e.g. Human: 3500000000)
# Output: For an analysis with X taxa and a genome size estimate of Y bp,  this script will subset each taxons reads down to (10*Y)/X bases (~10X total coverage when summed)

import os
from os import path
import sys
from glob import glob
import subprocess
from subprocess import check_call
import pandas as pd
import re

#Set cwd to script location
script_dir = sys.path[0]

#Get Genome size estimate
genomeSize = int(sys.argv[1])

#Set TrimRead directories based off of script folder location
trim_read_dir = path.dirname(path.abspath(script_dir))+"/Reads/TrimReads"

#Find taxa folders within TrimRead folder
trim_read_tax_dirs = sorted(glob(trim_read_dir+"/*/"))

#Create folder for Subset reads
if(not path.isdir(path.dirname(path.abspath(script_dir))+"/Reads/SubsetReads")):
    os.mkdir(path.dirname(path.abspath(script_dir))+"/Reads/SubsetReads")
subset_output_dir = path.dirname(path.abspath(script_dir))+"/Reads/SubsetReads/"

#Create folder for Subset output logs
if(not path.isdir(trim_read_dir+"/subsetOutput")):
    os.mkdir(trim_read_dir+"/subsetOutput")
subset_log_dir = trim_read_dir+"/subsetOutput/"

#Ensure Trim Log/FastQC/subset Log output directories not included in taxa list
trim_read_tax_dirs = [x for x in trim_read_tax_dirs if not x.endswith('trimOutput/')]
trim_read_tax_dirs = [x for x in trim_read_tax_dirs if not x.endswith('fastqcOutput/')]
trim_read_tax_dirs = [x for x in trim_read_tax_dirs if not x.endswith('subsetOutput/')]

#Calculate subset depth
subsetDepth = int((10*genomeSize)/len(trim_read_tax_dirs))

#Initialize Pandas DF to get base counts and lists of paired and single-end files
df = pd.DataFrame(columns=['Taxon','Dataset','Basecount'])
compiled_paired = list()
compiled_single_end = list()

#For each taxa directory...
for tax_dir in trim_read_tax_dirs:
    #List all files and set output dir
    files = sorted(glob(tax_dir+"*.fastq.gz"))
    taxon_ID = path.basename(tax_dir[:-1])

    left_pairs = list()
    right_pairs = list()
    single_end = list()

    taxon_list = list()
    taxon_ID = path.basename(tax_dir[:-1])

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

df = df.reset_index(drop=True)
df["Basecount"] = pd.to_numeric(df["Basecount"])
df = df.sort_values(by=['Taxon'])
df["ToSample"] = 0
taxa_list = sorted(df.Taxon.unique())

print("Based on a genome size estimate of " + str(genomeSize) + " bp, and with " + str(len(trim_read_tax_dirs)) + " species, the requested subset depth is " + str(subsetDepth) + " bp per species")

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

#Subset bases in single-end files if present...
if len(compiled_single_end) > 0:
    single_end_DF = df[df.Dataset.isin(compiled_single_end)]
    for row in single_end_DF.itertuples():
        subset_command = [
            'reformat.sh',
            'in={}'.format(trim_read_dir+"/"+str(row.Taxon)+"/"+str(row.Dataset)+".fastq.gz"),
            'out={}'.format(subset_output_dir+str(row.Dataset).replace('Trim','GenomeReads.fastq.gz')),
            'samplebasestarget={}'.format(int(row.ToSample)),
            'ow=t',
            '&>',
            '{outDir}out_{fileName}_Trim'.format(outDir=subset_log_dir,fileName=row.Dataset)]
        check_call(subset_command)

if len(compiled_paired) > 0:
    paired_DF = df[df.Dataset.isin(compiled_paired)]
    for row in paired_DF.itertuples():
        subset_command = [
            'reformat.sh',
            'in={}'.format(trim_read_dir+"/"+str(row.Taxon)+"/"+str(row.Dataset)+"_1.fastq.gz"),
            'in2={}'.format(trim_read_dir+"/"+str(row.Taxon)+"/"+str(row.Dataset)+"_2.fastq.gz"),
            'out={}'.format(subset_output_dir+str(row.Dataset).replace('Trim','GenomeReads_1.fastq.gz')),
            'out2={}'.format(subset_output_dir+str(row.Dataset).replace('Trim','GenomeReads_2.fastq.gz')),
            'samplebasestarget={}'.format(int(row.ToSample)),
            'ow=t',
            '&>',
            '{outDir}out_{fileName}_Trim'.format(outDir=subset_log_dir,fileName=row.Dataset)]
        check_call(subset_command)
