#!/usr/bin/env python3
import os
import sys
from collections import Counter
from Bio import SeqIO
import glob

#get combined pileup info
def getallbases(path):
    ''' 
    This function gets combined pileup information. 

    Arguments:
    path (string): path to folder containing pileup file

    Returns:
    dict: for each location in the composite genome we get one final base call
    '''

    assert len(glob.glob1(path,"*.pileups")) == 1, 'More than one pileup file in'+path
    allbases=dict()
    with open(path+'/'+os.path.basename(path)+'.pileups',"r") as filein:
        for line in filein:
            splitline = line.split()
            if len(splitline) > 4 and int(splitline[3]) > 0:
                node, pos, ref, num, bases, qual = line.split()
                loc = node+'/'+pos
                cleanBases = getCleanList(ref,bases)
                assert len(cleanBases) == int(num), 'bases are being counted incorrectly: '+ str(bases) + ' should have '+str(num)+' bases, but it is being converted to '+"".join(cleanBases)
                finalBase = getFinalBase_Specific(cleanBases)
                allbases[loc]=finalBase
    return allbases

def getCleanList(ref,bases):
    ''' 
    This function builds and returns clean list of bases without indels. 

    Arguments:
    ref (chr): the reference base
    bases (string): bases in the reads at the site from the pileup

    Returns: 
    list: bases from the pileup in uppercase, either A, C, G, T, no indels
    '''

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

def getFinalBase_Specific(cleanBases):
    ''' 
    This function returns the base if the site is homozygous or finds bases with highest read support. 

    Arguments:
    cleanBases (list): bases from the pileup in uppercase, either A, C, G, T, no indels 

    Returns:
    chr: final base call given list of input 
    '''

    baseCount = Counter(cleanBases)
    #If the site is homozygous, return the base
    if len(baseCount) == 1:
        finalBase = baseCount.most_common()[0][0]
        if finalBase == '*':
            finalBase = 'N'
    else:
        #Given possible bases from getCleanList, find bases with highest read support
        possibleBases = ['A','T','C','G','*']
        finalBases = []
        maxCount=int(baseCount.most_common()[0][1])
        for base in possibleBases:
            if baseCount[base] == maxCount:
                finalBases.append(base)

        #If one base remains, return final base
        if len(finalBases) == 1:
            if finalBases[0] == '*':
                finalBase = 'N'
            else:
                finalBase = finalBases[0]

        #If more than one base remains, return alpha sorted base with maxCount support
        if len(finalBases) > 1:
            if '*' in finalBases:
                finalBases.remove('*') #If * vs. base, return base
            finalBase = (sorted(finalBases))[0]
    return finalBase

    #If we can find a mapper to handle degenerate bases...
    #degenDict = {"ACGT":"N","GT":"K","AC":"M","AG":"R","CT":"Y","CG":"S","AT":"W","CGT":"B","AGT":"D","ACT":"H"}
    #counter=Counter(cleanBases)
    #baseList = sorted([key for key, _ in counter.most_common()])
    #if '*' in baseList:
    #    finalBase = 'N'
    #else:
    #    finalBase = degenDict[baseList]
    #return finalBase

###############################################
def getbases_main(path, contig_file):
    '''
    Write a contigs.fa file specif to the species

    Arguments:
    path (string): path to folder containing pileup file 
    contig_file (string): path to fa file containing contigs

    Returns: None

    '''
    allbases=getallbases(path)      #dictionary of combined pileups - locus/pos:bases(as list)
    if len(allbases)==0:
        print('No data for '+path, flush=True)
        sys.exit(1)

    # Read contig fasta file into dictionary with sequence ID as the key
    contig_handle = open(contig_file, "r")
    fasta_seq = SeqIO.parse(contig_handle, 'fasta')
    fasta_dict = {read.id:list(str(read.seq)) for read in fasta_seq}
    contig_handle.close()

    for locus_pos,base in sorted(allbases.items()):
        locus,pos = locus_pos.split('/')
        fasta_dict[locus][int(pos)-1] = base

    output = open(path+'/contigs.fa', 'w')
    for l,seq in sorted(fasta_dict.items()):
        output.write('>'+str(l)+"\n"+"".join(seq)+"\n")
    output.close()

if __name__ == '__main__':
    path = sys.argv[1]
    contig_file = sys.argv[2]
    getbases_main(path, contig_file)
