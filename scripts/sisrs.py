'''
Script Developed by Devin J. McConnell
Schwartz Lab May 2019
'''

from sisrs_01_folder_setup import *
from sisrs_02_read_trimmer import *
from cmdCheck import *

'''
Function to run all of the first script
'''
def sisrs01(taxon_ID_file,data_path,sisrs_dir,trimed):
    taxa_list = []
    if taxon_ID_file != "":
        taxa_list = devTaxaList(taxon_ID_file)
    else:
        taxa_list = readTaxaFromFolders(data_path)

    fileStructure(sisrs_dir, taxa_list)

    if data_path != "" :
        if trimed:
            makeLinks(data_path, sisrs_dir, taxa_list, True)
        else:
            makeLinks(data_path, sisrs_dir, taxa_list, False)

'''
Function to run all of the second script
'''
def sisrs2(trimed,processors,sisrs_dir):
    if not trimed:
        bbduk_adapter = findAdapter()
        out = setup(sisrs_dir)

        # raw_fastqc_command
        fastqcCommand(processors,out[3],out[0])

        trim(out[5],out[1],bbduk_adapter,out[2])

        # trim_fastqc_command
        fastqcCommand(processors,out[4],out[1])
    else:
        print("FILES ALREADY TRIMMED")


if __name__ == '__main__':
    cmdln = sys.argv
    rtn = commandLine(cmdln)
    sisrs01(rtn[1],rtn[2],rtn[0],rtn[3])
    sisrs2(rtn[3],rtn[4],rtn[0])
