#!/usr/bin/env python3

'''
Argument 1: <GenomeName>_MapData.tsv (From Bowtie Script)
'''

import sys
import os
import time
import pandas as pd

def blowUpCIGAR(cigar):

    expandedCIGAR = ""
    tempCount = ""
    for c in cigar:
        if c.isdigit():
            tempCount+=c
        else:
            expandedCIGAR+=int(tempCount)*c
            tempCount=""
    return expandedCIGAR

def getReferencePosition(commandList_Set,genomeName,output_dir):
    for commandList in commandList_Set:
        humanPosition = commandList[2]
        sisrsForwardPosition = 0
        sisrsReversePosition = commandList[5] - 1

        locList = []
        locDict = {}

        for x in range(1,commandList[5]+1):
            locList.append(commandList[0] +'/'+str(x))
        with open(output_dir + "/WholeGenome_" + genomeName + "_Mapped.bed",'a+') as printFile:
            if commandList[3] == 0:
                for c in commandList[4]:
                    if c=='M':
                        printFile.write(str(commandList[1])+"\t"+str(humanPosition-1)+"\t"+str(humanPosition)+"\t"+str(locList[sisrsForwardPosition])+"\n")
                        sisrsForwardPosition+=1
                        humanPosition+=1
                    elif c=='I':
                        sisrsForwardPosition+=1
                    elif c=='D':
                        humanPosition+=1

            elif commandList[3] == 16:
                for c in commandList[4]:
                    if c=='M':
                        printFile.write(str(commandList[1])+"\t"+str(humanPosition-1)+"\t"+str(humanPosition)+"\t"+str(locList[sisrsReversePosition])+"\n")
                        humanPosition+=1
                        sisrsReversePosition-=1
                    elif c=='I':
                        sisrsReversePosition-=1
                    elif c=='D':
                        humanPosition+=1

#Check arguments
assert len(sys.argv) == 2, "This script accepts one argument " + str(len(sys.argv)-1) + " were given."
assert sys.argv[1].split('_')[-1] == 'MapData.tsv', "The argument for this script must be the <GenomeName>_MapData.tsv file."

startTime = time.time()
#Read in output of Bowtie script (<>_MapData.tsv)
mapLength = pd.read_csv(sys.argv[1],dtype={0: object, 1: object, 2: int, 3: int, 4: object},sep='\t',header=None)
mapLength.columns=['Contig','Chrom','Start','Direction','CIGAR']
scafCount = mapLength.shape[0]

#Strip path and name of genome file
output_dir=os.path.dirname(os.path.abspath(sys.argv[1]))+'/Whole_Genome_Mapping'
genomeName = str(os.path.basename(sys.argv[1]).replace('_MapData.tsv',''))

#Expand CIGAR
mapLength['CIGAR'] = mapLength['CIGAR'].apply(blowUpCIGAR)

#Get length  of contig based off CIGAR
mapLength['Length'] = mapLength['CIGAR'].apply(lambda x: (x.count('M')+x.count('I')))
totalLength = sum(mapLength['Length'])

#Convert mapLength to commandList
commandList = mapLength.values.tolist()

#Create whole-genome map file (print directly to file)
print("Creating map dataset for whole genome...")
sys.stdout.flush()
startTime = time.time()
getReferencePosition(commandList,genomeName,output_dir)
print("\nComplete genome map file for " + genomeName + " was created at " + output_dir + "/WholeGenome_" + genomeName + "_Mapped.bed in " + str(format((time.time() - startTime),'0.0f')) + "s")
