'''
Script Developed by Devin J. McConnell
Schwartz Lab May 2019
'''

from sisrs_01_folder_setup import *
import subprocess

if __name__ == '__main__':
    cmdln = sys.argv

    taxon_ID_file = ""
    raw_data_path = ""

    # Flag if you plan on using a TaxonID txt file and
    # wish to move files by hand
    if '-id' in cmdln:
        taxon_ID_file = cmdln[cmdln.index('-id') + 1]

    # Flag to tell us were the raw data is stored
    if '-rd' in cmdln:
        raw_data_path = cmdln[cmdln.index('-rd') + 1]

    subprocess.call(["python scripts/sisrs_01_folder_setup.py %s"%raw_data_path],shell=True)
