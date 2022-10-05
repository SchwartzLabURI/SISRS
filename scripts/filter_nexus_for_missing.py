#!/usr/bin/env python3

"""
    Take a nex file, outputs phylip-relaxed with X missing data

    arguments:
    data  -- extension should be nex/phy/fa
    missing -- the number of species in the alignment allowed to have missing data
    gap_allowed - True/False

    output:
    phylip formatted file ending with _mX.phylip-relaxed where X is the number missing
"""

from __future__ import division
import sys
from os import path
import linecache
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment, AlignInfo
from Bio import AlignIO, SeqIO

def filter_nexus(alignment_filename, missing_list):
    '''
    Reads nexus, passes to filtering function that Writes phylip file reducing data to max missing per site

    Arguments:
    alignment_filename (string): path to alignment file (in theory can use different types but since it gets loc IDs from nexus format.... (fix this)
    missinglist (list): number of taxa allowed to be missing

    Returns: none

    '''
    formats = {'nex':'nexus', 'phy':'phylip-relaxed', 'fa':'fasta'}
    fformat = formats[alignment_filename.split('.')[-1]]
    data = SeqIO.to_dict(SeqIO.parse(alignment_filename, fformat))

    #Extract loc IDs from nexus
    locline = linecache.getline(alignment_filename, 7) # ASSUMES NEXUS (fix?)
    locs = locline.split()[1:-1] #Remove brackets from top and bottom

    for k in data:
        data[k] = list(data[k].seq)

    for missing in missing_list:
        for bases in [['A', 'C', 'G', 'T', 'a', 'c', 'g', 't'], ['A', 'C', 'G', 'T', 'a', 'c', 'g', 't', '-']]:
            filter_nexus1(alignment_filename, data, locs, int(missing), bases)

def filter_nexus1(alignment_filename, data, locs, missing, bases):
    '''
    Writes phylip file reducing data to max missing per site
    Prints info about output

    Arguments:
    alignment_filename (string): path to alignment file - not reading data in this function, just used for path
    data (dict): species: sequence as list
    missing (int): number of taxa allowed to be missing
    bases (list): list of bases - either includes gap or not

    Returns: none

    '''

    species = list(data.keys())
    minsp = len(species)-missing
    newlocs = list()

    newdata = {sp: list() for sp in species}
    for i in range(len(data[species[0]])):
        site = [data[sp][i] for sp in species if data[sp][i] in bases]
        if len(set(site)) > 1 and len(site) >= minsp:
            for sp in species:
                newdata[sp].append(data[sp][i])
            newlocs.append(locs[i])

    datalist = []
    for k,v in sorted(newdata.items()):
        seq = SeqRecord(Seq(''.join(v)), id=k)
        datalist.append(seq)

    if '-' in bases:
        SeqIO.write(datalist, path.dirname(alignment_filename) + '/'+ path.basename(alignment_filename).split('.')[0] + '_m'+str(missing)+'.phylip-relaxed', "phylip-relaxed")
        locfile = open(path.dirname(alignment_filename)+'/'+ path.basename(alignment_filename).replace('.nex','') + '_locs_m'+str(missing)+'.txt', 'w')    
    else:
        SeqIO.write(datalist, path.dirname(alignment_filename) + '/'+ path.basename(alignment_filename).split('.')[0] + '_m'+str(missing)+'_nogap.phylip-relaxed', "phylip-relaxed")
        locfile = open(path.dirname(alignment_filename)+'/'+ path.basename(alignment_filename).replace('.nex','') + '_locs_m'+str(missing)+'_nogap.txt', 'w')

    locfile.write("\n".join(newlocs))
    locfile.close()
    origLength = len(data[species[0]])
    newLength = len(newlocs)

    if '-' in bases:
        print('With '+str(missing)+' taxa allowed to be missing, '+str(origLength)+' sites from '+path.basename(alignment_filename)+' ('+str(len(species)-2)+' allowed missing) are reduced to '+str(len(newlocs))+' sites ('+str(origLength-newLength)+' sites or '+str('%.2f' % (((origLength-newLength)/origLength)*100))+'% lost)')
    else:
        print('With '+str(missing)+' taxa allowed to be missing, and with no gaps allowed, '+str(origLength)+' sites from '+path.basename(alignment_filename)+' ('+str(len(species)-2)+' allowed missing) are reduced to '+str(len(newlocs))+' sites ('+str(origLength-newLength)+' sites or '+str('%.2f' % (((origLength-newLength)/origLength)*100))+'% lost)')


if __name__ == '__main__':
    # Get arguments
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-f','--filename',action='store',nargs="?")
    my_parser.add_argument('-m', '--missing',type=int, action='store',default=[0],nargs="*")
    args = my_parser.parse_args()

    alignment_filename = args.filename
    missinglist = args.missing #contains a list of missing

    filter_nexus(alignment_filename, missinglist)

