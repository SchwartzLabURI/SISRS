#!/usr/bin/env python3

"""Take multiple mpileup files and determine the base at each position of the reference genome

    arguments:
        tax_dir: folder containing mpileup files ending in pileups (SISRS_Run/SPP)
        contig_dir: folder containing assembled composite genome (SISRS_Run/Composite_Genome)
        minread: number of reads at a position required to call a base (3 is a good start)
        thresh: proportion of reads that must be one base for calling to occur (1 means no intraspecies variation)

    output:
        tax_dir/<SP>_LocList: list of bases for each position
"""

import sys
from collections import Counter
from glob import glob
import string
import re
import os
from os import path
from collections import defaultdict

#get combined pileup info
def getallbases(tax_dir,contig_dir,minread,thresh):
    basePath=path.dirname(tax_dir)
    assert len(glob(tax_dir+"/*.pileups"))==1,'More than one pileup file in'+tax_dir
    speciesDict=defaultdict(lambda: 'N')
    minPenalty=0
    threshPenalty=0
    bothPenalty=0
    with open (tax_dir+"/"+path.basename(tax_dir)+'.pileups',"r") as filein:
        for line in iter(filein):
            splitline=line.split()
            if len(splitline)>4 & int(splitline[3])>0:
                node,pos,ref,num,bases,qual=line.split()
                loc=node+'/'+pos
                cleanBases=getCleanList(ref,bases)  #Get clean bases where * replaced with -
                finalBase,minPenalty,threshPenalty,bothPenalty=getFinalBase_Pruned(cleanBases,minread,thresh,minPenalty,threshPenalty,bothPenalty)
                speciesDict[loc] = finalBase

    printSpecies = open(tax_dir+"/"+path.basename(tax_dir)+'_LocList', 'w')
    with open(contig_dir+"/Composite_Genome_LocList") as f:
        for line in f:
            printSpecies.write(speciesDict[line.strip()]+'\n')
    f.close()
    printSpecies.close()

    c = Counter(speciesDict.values())
    nCount = c['N']
    siteCount = len(speciesDict) - nCount
    sitePercent = format((float(siteCount)/len(speciesDict))*100,'.2f')
    nPercent = format((float(nCount)/len(speciesDict))*100,'.2f')
    print("Of "+ str(len(speciesDict)) + " positions, " + path.basename(tax_dir) + " has good calls for " + str(siteCount) + " sites (" + sitePercent +"%). There were " + str(nCount) + " N calls ("+ nPercent + "%).",flush=True)
    print("Of " + str(nCount) + " Ns, " + path.basename(tax_dir) + " lost " + str(threshPenalty) + " via homozygosity threshold, " + str(minPenalty)  +" from low coverage, and " + str(bothPenalty) + " from both. "+ str(nCount - threshPenalty - minPenalty - bothPenalty) + " sites had no pileup data.\n",flush=True)

    return siteCount

def getCleanList(ref,bases):
    bases=bases.replace('.',ref) #insert ref base
    bases=bases.replace(',',ref)
    bases=bases.upper() #everything in uppercase
    bases=list(bases)

    okbases=['A','C','G','T','*']
    indels=['+','-']

    new_base_list=[]
    ibase_list = iter(bases)
    for b in ibase_list:
        if b in okbases:
            new_base_list.append(b)         #Get base
        elif b in indels:                   #skip indels
            i = int(next(ibase_list))
            j = str(next(ibase_list))
            if str.isdigit(j):
                skip=int(str(i)+j)
                while skip>0:
                    z=next(ibase_list)
                    skip=skip-1
            else:
                while i>1:
                    z=next(ibase_list)
                    i = i-1
        elif b=='^':                        #skip read qual noted at end of read
            z=next(ibase_list)

    return new_base_list

def getFinalBase_Pruned(cleanBases,minread,thresh,minPenalty,threshPenalty,bothPenalty):
    singleBase=(Counter(cleanBases).most_common()[0][0])
    if singleBase == '*':
        singleBase = '-'
    counts=int((Counter(cleanBases).most_common()[0][1]))

    if counts >= minread and counts/float(len(cleanBases)) >= thresh:
        finalBase=singleBase
    else:
        finalBase='N'
        if counts < minread and counts/float(len(cleanBases)) < thresh:
            bothPenalty+=1
        elif counts < minread:
                minPenalty+=1
        elif counts/float(len(cleanBases)) < thresh:
                threshPenalty+=1

    return finalBase,minPenalty,threshPenalty,bothPenalty
###############################################
def main(tax_dir, contig_dir, minread, thresh):

    allbases=getallbases(tax_dir,contig_dir,minread,thresh)      #dictionary of combined pileups - locus/pos:bases(as list)
    if allbases==0:
        print('No data for '+tax_dir,flush=True)
        sys.exit(1)

if __name__ == '__main__':
    tax_dir = sys.argv[1]
    contig_dir = sys.argv[2]

    #Ensure no trailing /
    if tax_dir[-1]=="/":
        tax_dir=tax_dir[0:-1]
    if contig_dir[-1]=="/":
        contig_dir=tax_dir[0:-1]
    minread = int(sys.argv[3])
    thresh = float(sys.argv[4])
    main(tax_dir, contig_dir, minread, thresh)
