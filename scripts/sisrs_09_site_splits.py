#!/usr/bin/env python3
'''
Argument 1: Unrooted tree file
Argument 2: Alignment file
Argument 3: Loc list for all sites in alignment
Argument 4: User-defined name for output file prefix
Argument 5: Directory for site splits
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
    outPath = os.path.abspath(os.path.dirname(tree_path))

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

    #Print split data
    with open(outPath+"/RefTree_SplitKey.txt", 'w') as outfile:
        goodSplits = []
        for k in natsorted(split1Dict.keys()):
            if isinstance(k,int):
                outfile.write('Split '+ str(k) + ': ' + split1Dict[k] + '...' + split2Dict[k]+'\n')
                goodSplits.append(split1Dict[k])
                goodSplits.append(split2Dict[k])
    print(" - Site split ID key written to " + outPath + "/RefTree_SplitKey.txt")

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

    if base1=='-' or base2=='-':
        gap = 1
    else:
        gap = 0

    speciesWithBase1 = ','.join(sorted(map(lambda sp: alignmentSpecies_Dict[sp],[j for j,x in enumerate(column) if x == base1])))
    speciesWithBase2 = ','.join(sorted(map(lambda sp: alignmentSpecies_Dict[sp],[j for j,x in enumerate(column) if x == base2])))

    if speciesWithBase1 in set(goodSplits_Dict.keys()):
        returnList = [loc,goodSplits_Dict[speciesWithBase1],gap]
    else:
        returnList = [loc,0,gap]
    return returnList

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
with open(sys.argv[3],'r') as f:
    allLocs = f.read().splitlines()
assert len(allLocs) == align_length, " - Length of loc list does NOT equal alignment length!"

columnList = []
for column in range(0,align_length):
    columnList.append([alignment[:,column],column])
results = map(goodOrBad,columnList)

outPath = sys.argv[5]

df = pd.DataFrame.from_records(results, columns=['Site','Split','Gap'])
siteID = str(sys.argv[4])
df.to_csv(outPath + '/' + siteID + "_AllSiteSplits.tsv",sep='\t',index=None,header=True)
print(" - Complete site split data written to " + outPath + "/" + siteID + "_AllSiteSplits.tsv")

dnaSplits = df[df['Gap']==0]
gapSplits = df[df['Gap']==1]

allDNASplits = dnaSplits.shape[0]
allGapSplits = gapSplits.shape[0]

if allDNASplits > 0:
    goodDNASplits = dnaSplits[dnaSplits['Split'] > 0].shape[0]
    badDNASplits = dnaSplits[dnaSplits['Split']  == 0].shape[0]
    goodDNAPercent = format((float(goodDNASplits)/allDNASplits)*100,'.2f')

    print("\n - Of " + str(align_length) + " total sites in the alignment, " + str(allDNASplits) + " alignment columns were gappless.")
    print(" - " + str(goodDNASplits) + " of these supported a split in the supplied tree, while " + str(badDNASplits) + " sites supported a non-canonical split. (" + str(goodDNAPercent) + "% good sites)")

else:
    print(" - No sites in the alignment were gappless.")

if allGapSplits > 0:
    goodGapSplits = gapSplits[gapSplits['Split'] > 0].shape[0]
    badGapSplits = gapSplits[gapSplits['Split']  == 0].shape[0]
    goodGapPercent = format((float(goodGapSplits)/allGapSplits)*100,'.2f')

    print("\n - Of " + str(align_length) + " total sites in the alignment, " + str(allGapSplits) + " columns contain indels/gaps.")
    print(" - " + str(goodGapSplits) + " gapped columns supported a split in the supplied tree, while " + str(badGapSplits) + " sites supported a non-canonical split. (" + str(goodGapPercent) + "% good sites)")

else:
    print(" - No sites in the alignment were gapped.")

print("\n - Split Support:\n")
splitCount = Counter(df['Split'])
splitCount_df =  pd.DataFrame(list(splitCount.items()), columns=['Split', 'Count']).sort_values(by=['Split'])

bad_dir = sys.argv[6]
good_dir = sys.argv[7]

for split in splitCount_df['Split']:
    if split==0:
        badDNASplits = list(dnaSplits[dnaSplits['Split']  == 0]['Site'])
        with open(bad_dir+"/Bad_Split_LocList.txt","w") as bad_file:
            for bad_split in badDNASplits:
                bad_file.write(bad_split+"\n")
        print("Non-canonical splits: " + str(splitCount[split]) + "\n")
    else:
        tempDNASplits = list(dnaSplits[dnaSplits['Split']  == split]['Site'])
        with open(good_dir+"/split_"+str(split)+"_LocList.txt","w") as good_file:
            for good_split in tempDNASplits:
                good_file.write(good_split+"\n")
        print("Split " + str(split) + " - "+split1Dict[split] + '...' + split2Dict[split] + " - " + str(splitCount[split]))
