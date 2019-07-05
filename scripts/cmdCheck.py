#!/usr/bin/env python3
import datetime
from os import path
'''
This function is designed to create a log file that contains all of the information
about the comand line from the previous run.
'''
def logFile(cmdln,rtn):
    name = '/sisrsLog.log'
    print(rtn[0] + name)
    f = open(rtn[0] + name, 'a+')

    currentDT = datetime.datetime.now()
    f.write("NEW SISRS RUN ON: {}\n".format(str(currentDT)))
    f.write("Command Line: {}\n".format(cmdln))
    f.write("SISRS Directory: {0} --> Full Path {1}\n".format(rtn[0],path.abspath(rtn[0])))
    f.write("Path to Data: {0} --> Full Path {1} \n".format(rtn[1],path.abspath(rtn[1])))
    f.write("Was the data trimmed (T=Y,F=N): {}\n".format(rtn[2]))
    f.write("Number of Processors Used: {}\n".format(rtn[3]))
    f.write("Genome Size: {}\n".format(rtn[4]))
    f.write("Threshold: {}\n".format(rtn[5]))
    f.write("Minimum Read Limit: {}\n".format(rtn[6]))
    f.write("Allowed Number of Missing Taxon (Single or Range): {}\n".format(rtn[7]))
    f.write("Adding New Taxons (T=Y,F=N): {}\n".format(rtn[8]))
    f.write("Adding New Sequences to Existing Taxons: {}\n".format(rtn[9]))

    f.close()

'''
Validation checks for string to int and float to prevent any errors
'''
def isInt(number):
    try:
        int(number)
        return True
    except:
        return False

def isRange(i):
    if '-' not in i:
        return False
    items = i.split('-')
    l = isInt(items[0])
    r = isInt(items[1])

    if l == False or r == False:
        return False
    else:
        return True

def isFloat(number):
    try:
        float(number)
        return True
    except:
        return False

'''
This function is used to determine which version of a flag was used
'''
def isFound(first,second,cmdln):
    data = ''
    try:
        data = cmdln[cmdln.index(first) + 1]
    except:
        data = cmdln[cmdln.index(second) + 1]
    return data
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
    # trimmed = False      --> 2
    # threads = 1         --> 3
    # genomeSize = 0      --> 4
    # threshold = 1       --> 5
    # minRead = 3         --> 6
    # missing = 0         --> 7
    # addTaxon = False    --> 8
    # addData = False     --> 9
    rtn = [ "" for i in range(10)]

    # Variablese used throughout if's
    bool = True

    # Let the user specify where to place the sisrs output
    if '-dir' in cmdln or '--directory' in cmdln:
        rtn[0] = isFound('-dir','--directory',cmdln)
    else:
        rtn[0] = script_dir

    # Flag to tell us were the data is stored
    if '-d' in cmdln or '--data' in cmdln:
        rtn[1] = isFound('-d','--data',cmdln)
    else:
        print("MISSING PATH TO DATA (-d,--data)")
        exit()

    # Determine if the data has been pretrimmed or not
    if '-trm' in cmdln or '--trimmed' in cmdln:
        rtn[2] = True
    else:
        rtn[2] = False

    # Determine the number of the threads to use
    if '-p' in cmdln or '--processors' in cmdln:
        bool = isInt(isFound('-p','--processors',cmdln))
        rtn[3] = int(isFound('-p','--processors',cmdln)) if bool else 1
    else:
        print("DEFAULT NUMBER OF PROCESSORS (-p,--processors) BEING USED: 1")
        rtn[3] = 1

    # Obtain the genomeSize estimation for script 3
    if '-gs' in cmdln or '--genomeSize' in cmdln:
        bool = isInt(isFound('-gs','--genomeSize',cmdln))
        rtn[4] = int(isFound('-gs','--genomeSize',cmdln)) if bool else 0
    else:
        if bool == False:
            print("GENOME SIZE ESTIMATION (-gs,--genomeSize) IS NOT A VALID NUMBER")
        else:
            if '-aT' not in cmdln and '--addTaxon' not in cmdln and '-aD' not in cmdln and '--addData' not in cmdln:
                print("GENOME SIZE ESTIMATION (-gs,--genomeSize) IS REQUEIRED FOR SISRS")
                exit()

    # Get threshold for script 5
    if '-thresh' in cmdln or '--threshold' in cmdln:
        bool = isFloat(isFound('-thresh','--threshold',cmdln))
        rtn[5] = float(isFound('-thresh','--threshold',cmdln)) if bool else 1
        if bool == False:
            print("MINIMUM SITE HOMOZYGOSITY FOR SISRS SITES WAS NOT A VALID NUMBER --> SWITCHED TO 1")
        elif rtn[5] <= 0 or rtn[5] > 1:
            print("MINIMUM SITE HOMOZYGOSITY FOR SISRS SITES MUST BE BETWEEN 0 AND 1  --> SWITCHED TO 1")
            rtn[5] = 1
    else:
        print("DEFAULT THRESHOLD IS BEING USED: 1")
        rtn[5] = 1

    # Gets the minread need to run the 5 script
    if '-mr' in cmdln or '--minread' in cmdln:
        bool = isInt(isFound('-mr','--minread',cmdln))
        rtn[6] = int(isFound('-mr','--minread',cmdln)) if bool else 3

        if bool == False:
            print("MINREAD COVERAGE (-mr,--minread) MUST BE A VALID INTEGER --> SWITCHED TO 3")
        elif rtn[6] < 1:
            print("MINREAD COVERAGE (-mr,--minread) MUST BE GREATER THAN 1 --> SWITCHED TO 3")
            rtn[6]=3
    else:
        print("DEFAULT MINREAD IS BEING USED: 3")
        rtn[6] = 3

    # Gets the allowed amount of missing data from the files
    if '-ms' in cmdln or '--missing' in cmdln:
        bool = isInt(isFound('-ms','--missing',cmdln))
        rtn[7] = int(isFound('-ms','--missing',cmdln)) if bool else 0
        if bool == False:
            bool2 = isRange(isFound('-ms','--missing',cmdln))
            rtn[7] = isFound('-ms','--missing',cmdln) if bool2 else 0
            if bool2 == False:
                print("THE ALLOWED AMOUNT OF MISSING DATA (-ms/--missing) MUST BE A VALID NUMBER --> SWITCHED TO 0")
                rtn[7] = 0
    else:
        print("THE DEFAULT AMOUNT OF MISSING DATA(-ms/--missing) IS BEING USED TO RUN SISRS: 0")
        rtn[7] = 0

    # Deal with adding a taxon
    if '-aT' in cmdln or '--addTaxon' in cmdln:
        rtn[8] = True
    else:
        rtn[8] = False

    # Deal with adding new data to an existing taxon
    if '-aD' in cmdln or '--addData' in cmdln:
        rtn[9] = True
    else:
        rtn[9] = False

    # Call to create the needed log file
    logFile(cmdln,rtn)

    return rtn
