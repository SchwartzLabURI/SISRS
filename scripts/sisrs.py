'''
Script Developed by Devin J. McConnell
Schwartz Lab May 2019
'''

from sisrs_01_folder_setup import *
from sisrs_02_read_trimmer import *
import subprocess

'''
Validation checks for string to int and float to prevent any errors
'''
def isInt(number):
    try:
        int(number)
        return True
    except:
        return False

def isFloat(number):
    try:
        float(number)
        return True
    except:
        return False

'''
Function that will parse the command line arguments. This might be better served
as a seperate file so this file does not get to large.
'''
def commandLine(cmdln):

    if len(cmdln) < 7:
        print("Requiered flags to run: -id or -rd, -gs, and -ms")
        exit()

    # Variables that will get passed to the scripts
    # Set defaults are shown
    # sisrs_dir = ""      --> 0
    # taxon_ID_file = ""  --> 1
    # data_path = ""      --> 2
    # trimed = False      --> 3
    # threads = 1         --> 4
    # genomeSize = 0      --> 5
    # threshold = 1       --> 6
    # minRead = 3         --> 7
    # missing = 0         --> 8
    rtn = [ "" for i in range(9)]

    # Variablese used throughout if's
    bool = True

    # Let the user specify where to place the sisrs output
    if '-dir' in cmdln:
        rtn[0] = cmdln[cmdln.index('-dir') + 1]
    else:
        rtn[0] = '.'

    # Flag if you plan on using a TaxonID txt file and
    # wish to move files by hand
    if '-id' in cmdln:
        rtn[1] = cmdln[cmdln.index('-id') + 1]
    # Flag to tell us were the raw data is stored
    elif '-rd' in cmdln:
        rtn[2] = cmdln[cmdln.index('-rd') + 1]
    else:
        print("MISSING TAXA INFORMATION")
        exit()

    # Determine if the data has been pretrimed or not
    rtn[3] = True if '-trm' in cmdln else False

    # Determine the number of the threads to use
    if '-th' in cmdln:
        bool = isInt(cmdln[cmdln.index('-th') + 1])
        rtn[4] = int(cmdln[cmdln.index('-th') + 1]) if bool else 1
    else:
        rtn[4] = 1

    # Obtain the genomeSize estimation for script 3
    if '-gs' in cmdln:
        bool = isInt(cmdln[cmdln.index('-gs') + 1])
        rtn[5] = int(cmdln[cmdln.index('-gs') + 1]) if bool else 0
    else:
        if bool == False:
            print("GENOME SIZE ESTIMATION IS NOT A VALID NUMBER")
        else:
            print("GENOME SIZE ESTIMATION IS REQUEIRED FOR SISRS")
        exit()

    # Get threshold for script 5
    if '-thrs' in cmdln:
        bool = isFloat(cmdln[cmdln.index('-thrs') + 1])
        rtn[6] = float(cmdln[cmdln.index('-thrs') + 1]) if bool else 1
        if bool == False:
            print("MINIMUM SITE HOMOZYGOSITY FOR SISRS SITES WAS NOT" +
                        "A VALID NUMBER --> SWITCHED TO 1")

    # Gets the minread need to run the 5 script
    if '-mr' in cmdln:
        bool = isInt(cmdln[cmdln.index('-mr') + 1])
        rtn[7] = int(cmdln[cmdln.index('-mr') + 1]) if bool else 3

        if bool == False:
            print("MINREAD COVERAGE MUST BE A VALID INTEGER --> SWITCHED TO 3")
        elif int(cmdln[cmdln.index('-mr') + 1]) < 1 or int(cmdln[cmdln.index('-mr') + 1]) > 3:
            print("MINREAD COVERAGE MUST BE BETWEEN 1 AND 3 --> SWITCHED TO 3")
    else:
        rtn[7] = 3

    # Gets the allowed amount of missing data from the files
    if '-ms' in cmdln:
        bool = isInt(cmdln[cmdln.index('-ms') + 1])
        rtn[8] = int(cmdln[cmdln.index('-ms') + 1]) if bool else 0
        if bool == False:
            print("THE ALLOWED AMOUNT OF MISSING DATA MUST BE A VALID NUMBER")
            exit()
    else:
        print("THE ALLOWED AMOUNT OF MISSING DATA IS REQUEIERED TO RUN SISRS")
        exit()

    return rtn

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


if __name__ == '__main__':
    cmdln = sys.argv
    rtn = commandLine(cmdln)
    sisrs01(rtn[1],rtn[2],rtn[0],rtn[3])
    sisrs2(rtn[3],rtn[4],rtn[0])
