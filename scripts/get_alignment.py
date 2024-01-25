#!/usr/bin/env python3

"""
    output alignment of sites useful for phylogenetics (nexus format)
    can specify max number of species missing for each site

    arguments:
        num_missing -- the number of species allowed to be missing at a site so that site will show in output alignment
        sisrs_dir -- Path to SISRS run directory (no traililng /)
        composite_dir -- Path to SISRS composite genome directory (no trailing /)

    output:
        alignment.nex : nexus formatted alignment including position in composite reference genome; each site as up to num_missing missing data
        alignment_bi.nex : above but only biallelic sites
        alignment_pi.nex : above but only phylogenetically informative sites (no singletons)

"""

import os
import re
import glob
from collections import Counter,defaultdict
import sys

#########################
def numsnps(site_dict_labels):
    '''
    print variable sites and how many singletons

    site_dict_labels: dict with label site as singleton, biallelic, or number of alleles
    '''

    print(str(len(site_dict_labels)) + ' total variable sites (alignment.nex)')

    singleton_count = sum(x == 1 for x in site_dict_labels.values())
    print(str(singleton_count)+' variable sites are singletons')

def get_base_info(linelist):
    '''

    linelist: data input as contig/site bases e.g. SISRS_contig-2/64 AAAAAAGGAAAN

    return:
        int: 0 for 2 or fewer samples or no variation; 1 for singletons; 2 for biallelic; or more bases
    '''
    bases = [b for b in linelist[1] if b in ['A', 'C', 'G', 'T', '-']]  # bases for that site
    if len(bases) > 2:  # looking for at least 3 species with data
        c = Counter(bases).most_common(5)
        if len(c) > 1:  # only count variable sites
            if c[1][1] == 1:  # there's only one second most common base ie the site has a singletone
                return 1  # note singleton
            elif len(c) == 2:  # biallelic but not singleton
                return 2
            else:
                return len(c)
        else:
            return 0
    else:
        return 0

def get_phy_sites(sisrs_dir,composite_dir,num_missing):
    ''' 
    This function saves to file all the sites with their sample info.
    
    Arguments:
    sisrs_dir (string):
    composite_dir (string): path to folder containing composite genome
    num_missing (int): max number of species allowed to be missing data for inclusion in output alignment 

    Returns:
    dict: site_dict_labels: each contig/site has a flag (0 for 2 or fewer samples or no variation; 1 for singletons; 2 for biallelic; or more bases)
    dict: species_data: each species has a sequence including only sites with flag>0

    '''

    #Fetch contig data
    contigList=glob.glob(composite_dir +'/contigs_LocList')

    assert len(contigList) > 0, 'Total site list not found in assembly folder'

    #Fetch sorted species data
    dataLists = sorted(glob.glob(sisrs_dir + '/*/*_LocList'))
    dataLists = [x for x in dataLists if 'contigs_LocList' not in x]
    splist=[os.path.basename(os.path.dirname(path)) for path in dataLists]
    speciesCount=len(dataLists)
    assert len(dataLists) > 0, 'No species had data from the pileup'

    allLists = contigList+dataLists

    site_dict_labels = {}
    species_data = {species: [] for species in splist}

    files = [open(i, "r") for i in allLists]
    outfile = open('sitedata', 'w')
    for rows in zip(*files):
        rowList = list(map(lambda foo: foo.replace('\n', ''), list(rows)))
        speciesData = rowList[1:(speciesCount+1)]
        newrow = filter_sites(speciesData, rowList,num_missing)
        if len(newrow)==2:  #we have data that is useful
            outfile.write(" ".join(newrow))
            outfile.write('\n')

            flag = get_base_info(newrow)
            if flag > 0:
                site_dict_labels[newrow[0]] = flag
                speciesData = list(newrow[1])
                for j in range(0, (len(splist))):
                    species_data[splist[j]].append(speciesData[j])
    outfile.close()
    return site_dict_labels,species_data

def filter_sites(speciesData, rowList,num_missing):
    newrow = []
    if speciesData.count("N") <= num_missing and len(set(filter(lambda a: a != "N", speciesData))) > 1: #enough data points & variation
        newrow.append(rowList[0])
        d = []
        for j in range(0, (len(speciesData))):
            d.append(speciesData[j])
            s = "".join(d)
        newrow.append(s)
    return newrow

def generic_write_alignment(fi, ntax, varsites, alignment_locations, species_data):
    spp = sorted(species_data.keys())
    ALIGNMENT=open(fi,'w')
    ALIGNMENT.write('#NEXUS\n\nBEGIN DATA;\nDIMENSIONS NTAX='+str(ntax)+' NCHAR='+str(varsites)+';\nFORMAT MISSING=? GAP=- DATATYPE=DNA;\nMATRIX\n')
    ALIGNMENT.write('[ '+ " ".join(alignment_locations)+' ]'+"\n")
    for species in spp: #write sequences for each species
        ALIGNMENT.write(species+"\t"+"".join(species_data[species])+"\n")
    ALIGNMENT.write(';\nend;')
    ALIGNMENT.close()

def write_alignment(fi,site_dict_labels, species_data):
    ''' 
    Writes an Alignment to nexus files - all data; bi locs only; pi locs only. 

    Arguments:
    fi (string) : filename to write to
    alignment (Alignment): Alignment to be written

    Returns: none
    '''

    spp = species_data.keys()
    ntax = len(spp)
    varsites = len(site_dict_labels)
    alignment_locations = sorted(site_dict_labels.keys())

    #Process alignment.nex
    generic_write_alignment(fi, ntax, varsites, alignment_locations, species_data)

    #Process alignment_bi.nex
    biallelic_num = sum(x == 2 for x in site_dict_labels.values())
    loc = [k for k,v in site_dict_labels.items() if v == 2 or v==1]
    flags = []
    for i in sorted(site_dict_labels.keys()):
        if site_dict_labels[i] == 2:
            flags.append(1)
        elif site_dict_labels[i] == 1:
            flags.append(1)
        else:
            flags.append(0)
    sp_data = {}
    for species in spp:
        sp_data[species] = [species_data[species][i] for i in range(len(flags)) if flags[i] == 1]
    generic_write_alignment(fi.replace('.nex','_bi.nex'), ntax, len(loc), loc, sp_data)
    print(str(len(loc)) + ' total biallelic sites including singletons (alignment_bi.nex)')

    #Process alignment_pi.nex
    singletons = sum(x == 1 for x in site_dict_labels.values())
    loc = [k for k, v in site_dict_labels.items() if v > 1] #keep non singletons
    flags = []
    for i in sorted(site_dict_labels.keys()):
        if site_dict_labels[i] >1:
            flags.append(1)
        else:
            flags.append(0)
    sp_data = {}
    for species in spp:
        sp_data[species] = [species_data[species][i] for i in range(len(flags)) if flags[i] == 1]

    generic_write_alignment(fi.replace('.nex', '_pi.nex'), ntax, len(loc), loc, sp_data)
    print(str(len(loc))+' total variable sites excluding singletons (alignment_pi.nex)')


#########################

def main(num_missing, sisrs_dir, composite_dir):
    site_dict_labels,species_data = get_phy_sites(sisrs_dir,composite_dir,num_missing)
    numsnps(site_dict_labels) #prints numbers of snps, biallelic snps, and singletons
    write_alignment(sisrs_dir+'/alignment.nex',site_dict_labels, species_data)

if __name__ == '__main__':
    num_missing = int(sys.argv[1])
    sisrs_dir = sys.argv[2]
    composite_dir = sys.argv[3]
    main(num_missing, sisrs_dir, composite_dir)

