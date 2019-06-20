#!/usr/bin/env python3
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
def commandLine(cmdln,script_dir):

    if len(cmdln) < 5:
        print("Requiered flags to run: -rd and -gs")
        exit()

    # Variables that will get passed to the scripts
    # Set defaults are shown
    # sisrs_dir = ""      --> 0
    # data_path = ""      --> 1
    # trimed = False      --> 2
    # threads = 1         --> 3
    # genomeSize = 0      --> 4
    # threshold = 1       --> 5
    # minRead = 3         --> 6
    # missing = 0         --> 7
    # addTaxon = False    --> 8
    #addData = False      --> 9
    rtn = [ "" for i in range(10)]

    # Variablese used throughout if's
    bool = True

    # Let the user specify where to place the sisrs output
    if '-dir' in cmdln:
        rtn[0] = cmdln[cmdln.index('-dir') + 1]
    else:
        rtn[0] = script_dir

    # Flag to tell us were the raw data is stored
    if '-rd' in cmdln:
        rtn[1] = cmdln[cmdln.index('-rd') + 1]
    else:
        print("MISSING TAXA INFORMATION")
        exit()

    # Determine if the data has been pretrimed or not
    rtn[2] = True if '-trm' in cmdln else False

    # Determine the number of the threads to use
    if '-th' in cmdln:
        bool = isInt(cmdln[cmdln.index('-th') + 1])
        rtn[3] = int(cmdln[cmdln.index('-th') + 1]) if bool else 1
    else:
        print("DEFAULT NUMBER OF THREADS BEING USED: 1")
        rtn[3] = 1

    # Obtain the genomeSize estimation for script 3
    if '-gs' in cmdln:
        bool = isInt(cmdln[cmdln.index('-gs') + 1])
        rtn[4] = int(cmdln[cmdln.index('-gs') + 1]) if bool else 0
    else:
        if bool == False:
            print("GENOME SIZE ESTIMATION IS NOT A VALID NUMBER")
        else:
            print("GENOME SIZE ESTIMATION IS REQUEIRED FOR SISRS")
        exit()

    # Get threshold for script 5
    if '-thrs' in cmdln:
        bool = isFloat(cmdln[cmdln.index('-thrs') + 1])
        rtn[5] = float(cmdln[cmdln.index('-thrs') + 1]) if bool else 1
        if bool == False:
            print("MINIMUM SITE HOMOZYGOSITY FOR SISRS SITES WAS NOT" +
                        "A VALID NUMBER --> SWITCHED TO 1")
        elif rtn[5] <= 0 and rtn[5] >= 1:
            print("MINIMUM SITE HOMOZYGOSITY FOR SISRS SITES MUST BE " +
                        "BETWEEN 0 AND 1  --> SWITCHED TO 1")
    else:
        print("DEFAULT THRESHOLD IS BEING USED: 1")
        rtn[5] = 1

    # Gets the minread need to run the 5 script
    if '-mr' in cmdln:
        bool = isInt(cmdln[cmdln.index('-mr') + 1])
        rtn[6] = int(cmdln[cmdln.index('-mr') + 1]) if bool else 3

        if bool == False:
            print("MINREAD COVERAGE MUST BE A VALID INTEGER --> SWITCHED TO 3")
        elif rtn[6] < 1 or rtn[6] > 3:
            print("MINREAD COVERAGE MUST BE BETWEEN 1 AND 3 --> SWITCHED TO 3")
    else:
        print("DEFAULT MINREAD IS BEING USED: 3")
        rtn[6] = 3

    # Gets the allowed amount of missing data from the files
    if '-ms' in cmdln:
        bool = isInt(cmdln[cmdln.index('-ms') + 1])
        rtn[7] = int(cmdln[cmdln.index('-ms') + 1]) if bool else 0
        if bool == False:
            print("THE ALLOWED AMOUNT OF MISSING DATA MUST BE A VALID NUMBER")
            exit()
    else:
        print("THE DEFAULT AMOUNT OF MISSING DATA IS BEING USED TO RUN SISRS: 0")
        rtn[7] = 0

    # Deal with adding a taxon
    rtn[8] = True if '-aT' in cmdln else False

    # Deal with adding new data to an existing taxon
    rtn[9] = True if '-aD' in cmdln else False

    return rtn
