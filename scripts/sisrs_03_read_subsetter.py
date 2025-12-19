#!/usr/bin/env python3

'''
Subset reads to 10x total for assembly into a composite genome

IMPORTANT: Trimmed reads should be QC checked prior to this step to remove
datasets with low quality 
This script performs read subsetting for a SISRS composite genome assembly
This script calls reformat.sh  (part of BBMap Suite)

Paired-end read files must be identically basenamed and end in _Trim_1.fastq.gz/_Trim_2.fastq.gz
Single-end read files must be end in _Trim.fastq.gz

Input: (1) -gs/--genomesize Genome size estimate in basepairs (e.g. Human: 3500000000)
Output: For an analysis with X taxa and a genome size estimate of Y bp,  this script
will subset each taxons reads down to (10*Y)/X bases (~10X total coverage when summed)

'''

from os import path
from glob import glob
import pandas as pd
import re
from pathlib import Path
import shutil
import argparse
import os
import subprocess
from sisrs_08_filter_contigs import get_taxon_list

def setupDir(sisrs_dir,genomeSize):
    '''
    This function is designed to do the remaining folder setup for the read subsetting.
    It also finishes up the other minor setups needed for this script.

    Arguments: 
      sisrs_dir (string): path to the output directory
      genomeSize (int): genome size estimate in basepairs.

    Returns: 
      trim_read_dir (string)
      subset_output_dir (string)
      subset_log_dir (string) 
      trim_read_tax_dirs (list)
      subsetDepth (int)
    '''

    #Set TrimRead directories based off of script folder location
    if path.isdir(sisrs_dir + "/Reads/TrimReads_nocl"):
        trim_read_dir = sisrs_dir + "/Reads/TrimReads_nocl"
    else:
        trim_read_dir = sisrs_dir + "/Reads/TrimReads"

    #taxa come from taxonlist file
    taxa = get_taxon_list(sisrs_dir+'/TaxonList.txt')

    #Find taxa folders within TrimRead folder
    trim_read_tax_dirs = [trim_read_dir+"/"+x for x in taxa]
    print('The following folders contained trimmed reads that will be subset for composite genome assembly')
    print(trim_read_tax_dirs)

    #Create folder for Subset reads
    if(not path.isdir(sisrs_dir +"/Reads/SubsetReads")):
        os.mkdir(sisrs_dir+"/Reads/SubsetReads")
    subset_output_dir = sisrs_dir+"/Reads/SubsetReads/"

    #Create folder for Subset output logs
    if(not path.isdir(trim_read_dir+"/subsetOutput")):
        os.mkdir(trim_read_dir+"/subsetOutput")
    subset_log_dir = trim_read_dir+"/subsetOutput/"

    #Calculate subset depth
    subsetDepth = int((10*genomeSize)/len(trim_read_tax_dirs))

    return trim_read_dir, subset_output_dir, subset_log_dir, trim_read_tax_dirs, subsetDepth

def countBases(f, paired): 
    '''
    This function counts bases in single-end files or paired-end files

    Arguments: 
    f (string): path to a SE file or a pair of files (minus "_Trim_1.fastq.gz")
    paired (bool): whether PE or SE

    Returns: 
    int: basecount

    '''
    
    if not paired:
        count_command = ['reformat.sh', f'in={f}_Trim.fastq.gz']
    else:
        count_command = ['reformat.sh', f'in={f}_Trim_1.fastq.gz', f'in2={f}_Trim_2.fastq.gz']

    result = subprocess.run(count_command, capture_output=True, text=True, check=True)

    count = re.search(r'HERE (.*) HERE_2', 
                  result.stderr.replace('(100.00%) \t', 'HERE ')
                                     .replace(' bases (100.00%)\n\nTime:', ' HERE_2'))
    baseCount = count.group(1)

    return baseCount
  
def countHelper(tax_dir):
    '''
    This helper function is designed to break up the load that is carried by countBasePair().

    Arguments: 
    tax_dir (string): path to taxa directory

    Returns: 
    dataframe with 'Taxon','Dataset','Basecount', 'Paired' for one taxon
    
    '''

    files = sorted(glob(tax_dir+"/*.fastq.gz")) #all fastqs
    taxon_ID = path.basename(tax_dir)
    data_info_df = pd.DataFrame(columns=['Taxon','Dataset','Basecount','Paired'])
    df_index = 0 #need to keep track of row in df

    #Find files ending in _1/_2.fastq.gz
    left_files = [s for s in files if "_Trim_1.fastq.gz" in s]
    right_files = [s for s in files if "_Trim_2.fastq.gz" in s]

    #Strip _1_Trim.fastq.gz/_2_Trim.fastq.gz and identify pairs based on file name
    left_files = [x.replace('_Trim_1.fastq.gz', '') for x in left_files]
    right_files = [x.replace('_Trim_2.fastq.gz', '') for x in right_files]
    paired_files = list(set(left_files).intersection(right_files)) #one fastq w/o extension per pair, includes path
    
    for pair in paired_files:
        baseCount = countBases(pair, True)  #currently has full path but no extension
        df_new_row = pd.DataFrame({
            'Taxon': taxon_ID, 
            'Dataset': path.basename(pair), 
            'Basecount': baseCount, 
            'Paired': 'True',
            }, index=[df_index]) # Create the new row as its own dataframe
        data_info_df = pd.concat([data_info_df, df_new_row]) #concatenate dfs
        df_index+=1

    #Reset file names and filter out single-end files
    left_pairs = list()
    right_pairs = list()
    for pair in paired_files:
        left_pairs.append(pair+"_Trim_1.fastq.gz") #actual file names of lefts with pairs
        right_pairs.append(pair+"_Trim_2.fastq.gz") #actual file names of rights with pairs

    single_end = [re.sub(r"_Trim.*",'',x) for x in files if x not in left_pairs+right_pairs] #all remaining files without a match
    for s in single_end:
        baseCount = countBases(s, False)  #currently has full path but no extension
        df_new_row = pd.DataFrame({
            'Taxon': taxon_ID, 
            'Dataset': path.basename(s), 
            'Basecount': baseCount, 
            'Paired': 'False',
            }, index=[df_index]) # Create the new row as its own dataframe
        data_info_df = pd.concat([data_info_df, df_new_row]) #concatenate dfs
        df_index+=1

    return data_info_df

def countBasePair(trim_read_tax_dirs): 
    '''
    This function counts bases in single-end files and paired-end files if present.

    Arguments: 
    trim_read_tax_dirs (string): path to trimmed reads taxa directories

    Returns: 
    dataframe with 'Taxon','Dataset','Basecount', 'Paired' for all taxa

    '''

    data_info_df = pd.DataFrame(columns=['Taxon','Dataset','Basecount','Paired'])
  
    #For each taxon directory get data and add to dataframe
    for tax_dir in trim_read_tax_dirs:
        df_new_row = countHelper(tax_dir) #each has file names without dir or _Trim etc
        data_info_df = pd.concat([data_info_df, df_new_row]) #concatenate dfs
  
    return data_info_df


def checkCoverage(trim_read_dir,df,subsetDepth,subset_output_dir):
    '''
    This function figures out how many bases should be sampled from each sample from each taxon (some don't have enough so we need to compensate from other samples).

    Arguments: 
    df: dataframe with 'Taxon','Dataset','Basecount', 'Paired'
    subsetDepth (int): number of reads needed per taxon
    subset_output_dir (string): path to subset output directory

    Returns: 
    df: dataframe with 'Taxon','Dataset','Basecount', 'Paired',"ToSample"

    '''

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

    return df

def subset(df,trim_read_dir,subset_output_dir,subset_log_dir):
    '''
    Subsets fastq files based on specified values. This is a general
    function and can handle single and paired ends reads.

    Arguments: 
    df: dataframe with 'Taxon','Dataset','Basecount', 'Paired',"ToSample"
    trim_read_dir (string): path to trimmed reads directory
    subset_output_dir (string): path to subset output directory
    subset_log_dir (string): path to subset log directory

    Returns: None.

    '''

    for i in range(len(df)):
        ToSample = int(df.loc[i, "ToSample"])
        Basecount = int(df.loc[i, "Basecount"])
        Taxon = str(df.loc[i, "Taxon"])
        Dataset = str(df.loc[i, "Dataset"])
        if ToSample < Basecount:
            if df.loc[i, "Paired"] == 'True':
                subset_command = [
                    'reformat.sh',
                    f'in={trim_read_dir}/{Taxon}/{Dataset}_Trim_1.fastq.gz',
                    f'in2={trim_read_dir}/{Taxon}/{Dataset}_Trim_2.fastq.gz',
                    f'out={subset_output_dir}{Dataset}_GenomeReads_1.fastq.gz',
                    f'out={subset_output_dir}{Dataset}_GenomeReads_1.fastq.gz',
                    f'samplebasestarget={ToSample}',
                    'ow=t']
            else:
                subset_command = [
                    'reformat.sh',
                    f'in={trim_read_dir}/{Taxon}/{Dataset}_Trim.fastq.gz',
                    f'out={subset_output_dir}{Dataset}_GenomeReads.fastq.gz',
                    f'samplebasestarget={ToSample}',
                    'ow=t']
            output_file = f'{subset_log_dir}out_{Dataset}'
            with open(output_file, 'w') as f:
                subprocess.run(subset_command, stdout=f, stderr=subprocess.STDOUT, check=True)
        else:
            trim_dir = Path(trim_read_dir) / Taxon
            for file in trim_dir.glob(f'{Dataset}_Trim*'):
                shutil.copy2(file, subset_output_dir)

def run3(sis, gs):
    trim_read_dir, subset_output_dir, subset_log_dir, trim_read_tax_dirs, subsetDepth = setupDir(sis, gs)

    data_info_df = countBasePair(trim_read_tax_dirs)

    print("Based on a genome size estimate of " + str(gs) + " bp, and with " + str(
        len(trim_read_tax_dirs)) + " species, the requested subset depth is " + str(subsetDepth) + " bp per species")

    data_info_df = checkCoverage(trim_read_dir, data_info_df, subsetDepth, subset_output_dir)

    print(data_info_df)

    subset(data_info_df, trim_read_dir, subset_output_dir, subset_log_dir)

if __name__ == '__main__':

    # Get arguments
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-d','--directory',action='store',nargs="?")
    my_parser.add_argument('-gs','--genome_size',type=int, action='store',nargs="?")
    args = my_parser.parse_args()

    sis = args.directory
    gs = args.genome_size

    run3(sis, gs)
