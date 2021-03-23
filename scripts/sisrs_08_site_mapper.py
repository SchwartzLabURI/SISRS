#!/usr/bin/env python3

'''
Argument 1: <>.bed (BED file with whole genome mapping information)
Argument 2: List of desired locs to locate
Argument 3: Dataset name for output files
'''

import sys
import os
import time
import pandas as pd

#Read in whole genome map file
allData = pd.read_csv(sys.argv[1],header=None,dtype={0: object, 1: int, 2: int, 3: object},sep="\t")
allData.columns=['Chrom','Start','End','Name']
allSites = allData.shape[0]
dataPath=os.path.dirname(os.path.abspath(sys.argv[1]))
sisrs_output_dir=dataPath.replace('Whole_Genome_Mapping','SISRS_Sites')
genomeName = os.path.basename(sys.argv[1]).replace("WholeGenome_","").replace("_Mapped.bed","")
siteID = sys.argv[3]
#Read in loc list of desired sites to locate in  reference genome
with open(sys.argv[2],'r') as f:
    sitesOfInterest= f.read().splitlines()

#Filter out sites from original genome data that map to the same location
print("\nFiltering out sites that map to the same genomic location...")

nonDupData = allData.drop_duplicates(subset=['Chrom', 'Start'], keep=False)
nonDupSites = nonDupData.shape[0]
nonDupData.to_csv(sys.argv[1].replace("_Mapped.bed","_Mapped_NonDup.bed"),sep='\t',header=False,index=False)

print("\n- Of " + "{:,}".format(allSites) + " sites that could be genomically mapped, " + "{:,}".format(allSites - nonDupSites) + " sites mapped to overlapping reference positions and were removed.\n")
print("A genome map file for " + genomeName + " without overlapping sites was created at " + dataPath + "/WholeGenome" + genomeName + "_Mapped_NonDup.bed\n")

#Delete whole dataset
allData = None
del allData

#Filter all data based on provided site list
print("Beginning site selection...\n")
nonDupData = nonDupData[nonDupData['Name'].isin(sitesOfInterest)]
foundSites = nonDupData.shape[0]

#Print output files
nonDupData.to_csv(sisrs_output_dir + "/" + genomeName + "_" + siteID + "_Mapped_NonDup.bed",sep='\t',header=False,index=False)
nonDupData['Name'].to_csv(sisrs_output_dir + "/" + genomeName + "_" + siteID + "_Mapped_NonDup_LocList.txt",sep='\t',header=False,index=False)

print("- Of " + "{:,}".format(len(sitesOfInterest)) + " requested sites from the " + siteID + " dataset, " + "{:,}".format(foundSites) + " sites were able to be genomically mapped based off non-overlapping data from " + str(sys.argv[1]) + ".")
print("- Created BED file for these sites at " + sisrs_output_dir + "/" + genomeName + "_" + siteID + "_Mapped_NonDup.bed")
print("- Created a loc list for these sites at " + sisrs_output_dir + "/" + genomeName + "_" + siteID + "_Mapped_NonDup_LocList.txt")
