'''
Script Developed by Devin J. McConnell
Schwartz Lab May 2019
'''

from sisrs_01_folder_setup import *
import subprocess

if __name__ == '__main__':
    cmdln = sys.argv

    taxon_ID_file = ""
    data_path = ""
    sisrs_dir = ""
    taxa_list = ""

    # Let the user specify where to place the sisrs output
    if '-dir' in cmdln:
        sisrs_dir = cmdln[cmdln.index('-dir') + 1]
    else:
        sisrs_dir = '.'

    ############################################################################
    '''
                    FLAGS AND PROCCESING FOR SISRS_01
            Modify this to look cleaner... I am thinking of doing All
            Flags in one function and return the results of the flags
            after that...
    '''

    # Flag if you plan on using a TaxonID txt file and
    # wish to move files by hand
    if '-id' in cmdln:
        taxon_ID_file = cmdln[cmdln.index('-id') + 1]
        taxa_list = devTaxaList(taxon_ID_file)

    # Flag to tell us were the raw data is stored
    if '-rd' in cmdln:
        data_path = cmdln[cmdln.index('-rd') + 1]
        taxa_list = readTaxaFromFolders(data_path)

    # Calls the filestructure function in sisrs_01
    fileStructure(sisrs_dir, taxa_list)

    # Determine if the data has been pretrimed or not
    if '-t' not in cmdln:
        if data_path != "" :
            makeLinks(data_path, sisrs_dir, taxa_list, False)
    else:
        if data_path != "" :
            makeLinks(data_path, sisrs_dir, taxa_list, True)
