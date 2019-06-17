#!/usr/bin/env python3
'''
Argument 1: Unrooted tree file
Argument 2: Alignment file
Argument 3: Loc list for all sites in alignment
Argument 4: User-defined name for output file prefix
Argument 5: Path to BadSplit dir
'''

from Bio import AlignIO,Align
from Bio.Align import AlignInfo
import pandas as pd
import os
import sys
from ete3 import Tree
from natsort import natsorted
from collections import Counter

def getSplits(tree_path):
    #Read in reference tree and define splits
    refTree = Tree(tree_path, format=1)

    speciesList = list()
    #Get species list
    for leaf in refTree:
        speciesList.append(leaf.name)

    #Get split data from tree
    nodeList=list()
    for node in refTree.traverse("preorder"):
        tempList=list()
        if len(node) > 1 and len(node) < len(refTree):
            for leaf in node:
               tempList.append(leaf.name)
               tempList=sorted(tempList)
            nodeList.append(tempList)
    #nodeList.pop(0)

#Create paired list for other side of split
    partnerList = list()
    for group in nodeList:
        partnerList.append(sorted(list(set(speciesList) - set(group))))

    #Create two-way dictionary defining each tree split
    split1Dict={}
    split2Dict={}
    for i in range(len(nodeList)):
        split1Dict[','.join(sorted(set(nodeList[i])))]=i+1
        split1Dict[i+1]=','.join(sorted(set(nodeList[i])))
        split2Dict[','.join(sorted(set(partnerList[i])))]=i+1
        split2Dict[i+1]=','.join(sorted(set(partnerList[i])))

    goodSplits = []
    for k in natsorted(split1Dict.keys()):
        if isinstance(k,int):
            goodSplits.append(split1Dict[k])
            goodSplits.append(split2Dict[k])

    goodSplits_Dict = {}
    for i in goodSplits:
        if i in split1Dict.keys():
            goodSplits_Dict[i] = split1Dict[i]
        else:
            goodSplits_Dict[i] = split2Dict[i]
    return speciesList, goodSplits_Dict,split1Dict,split2Dict

def goodOrBad(columnList):

    global alignmentSpecies_Dict
    global goodSplits_Dict
    global allLocs

    column = list(columnList[0])
    loc = allLocs[columnList[1]]

    base1 = list(set(column))[0]
    base2 = list(set(column))[1]

    speciesWithBase1 = ','.join(sorted(map(lambda sp: alignmentSpecies_Dict[sp],[j for j,x in enumerate(column) if x == base1])))
    speciesWithBase2 = ','.join(sorted(map(lambda sp: alignmentSpecies_Dict[sp],[j for j,x in enumerate(column) if x == base2])))

    returnList = [loc,speciesWithBase1,speciesWithBase2]
    return returnList

if __name__ == "__main__":

    #SYS.ARGV[1]
    #Read in tree file
    tree_path=sys.argv[1]

    treeSpecies,goodSplits_Dict,split1Dict,split2Dict = getSplits(tree_path)

    #SYS.ARGV[2]
    #Read in alignment
    fformat = os.path.basename(sys.argv[2].split('.')[-1])
    formats = {'nex':'nexus', 'phy':'phylip-relaxed', 'fa':'fasta','phylip-relaxed':'phylip-relaxed'}
    alignment = AlignIO.read(sys.argv[2], formats[fformat])

    #Get alignment info (length, summary, species)
    align_length = alignment.get_alignment_length()
    summary = AlignInfo.SummaryInfo(alignment)

    alignmentSpecies = []
    alignmentSpecies_Dict = {}
    i=0
    for record in alignment:
        alignmentSpecies.append(record.id)
        alignmentSpecies_Dict[i]=record.id
        i+=1
    assert set(alignmentSpecies) == set(treeSpecies), " - Species from tree do NOT match species in alignment!"

    #SYS.ARGV[3]
    outPath = os.path.abspath(os.path.dirname(sys.argv[3]))

    with open(sys.argv[3],'r') as f:
        allLocs = f.read().splitlines()
    assert len(allLocs) == align_length, " - Length of loc list does NOT equal alignment length!"

    siteID = str(sys.argv[4])

    outPath = sys.argv[5]

    columnList = []
    for column in range(0,align_length):
        columnList.append([alignment[:,column],column])
    results = map(goodOrBad,columnList)

    df = pd.DataFrame.from_records(results, columns=['Site','One','Two'])
    df2=df.drop(['Site'],axis=1)
    df2['Set']=df2.values.tolist()
    df2.apply(lambda row: row.Set.sort(),axis=1)
    df2['String']=df2.apply(lambda row:"...".join(row.Set),axis=1)
    badCount=Counter(df2.String)
    splitList=[]
    countList=[]
    for x in badCount:
        splitList.append(x)
        countList.append(badCount[x])
    df3 = pd.DataFrame({'Split': splitList,'Split_Count': countList})
    df3.sort_values(['Split_Count'],ascending=False,inplace=True)
    df3['Split_Number'] = range(1,len(splitList)+1)
    cols = ['Split_Number','Split','Split_Count']
    df3 = df3[cols]

    df3[['One','Two']] = df3.Split.str.split(pat = "\.\.\.",expand=True)
    split_lookup_dict = {**df3.set_index('One')['Split_Number'].to_dict(),**df3.set_index('Two')['Split_Number'].to_dict()}
    df3=df3[['Split_Number','Split','Split_Count']]
    df['Split_Number'] = df.replace({"One": split_lookup_dict})['One']
    df = df[['Site','Split_Number']]

    for value in sorted(set(df['Split_Number'])):
        df[df['Split_Number']==value]['Site'].to_csv(outPath+"/bad_split_"+str(value)+"_LocList.txt",header=None,index=None)

    df.to_csv(outPath + '/' + siteID + ".tsv",sep='\t',index=None,header=None)
    df3.to_csv(outPath + '/' + siteID + "_Counts.tsv",sep='\t',index=None,header=None)
