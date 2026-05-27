#!/usr/bin/env python3

'''
Trim reads using BBDuk

Paired-end read files must be identically basenamed and end in _1/_2
Input: (Optional) -p/--processors Number of available processors for FastQC
(Default: 1, For 10 processors: python sisrs_02_read_trimmer.py -p 10)
# Output: (1) Trimmed Reads (in <base_dir>/Reads/TrimReads/<Taxon>/<Read_Basename>_Trim.fastq.gz
# Output: (2) Trim Log (in <base_dir>/Reads/RawReads/trimOutput)
# Output: (3) FastQC output for all raw + trimmed read sets (in <base_dir>/Reads/RawReads/fastqcOutput & <base_dir>/Reads/TrimReads/fastqcOutput)
'''

import os
from os import path, listdir
import sys
from glob import glob
import subprocess
import argparse
from sisrs_08_filter_contigs import get_taxon_list
from pathlib import Path
import resource
import psutil
import time
import gc


def findAdapter():
    '''
    This function is designed to find BBDuk + Adapter File.

    Arguments: none.

    Returns:  
    string: the path to the location of the adapter file.
    '''

    bbduk_cmd = ['which', 'bbduk.sh']
    result = subprocess.run(bbduk_cmd, capture_output=True, text=True, check=True)
    bbduk_adapter = path.dirname(result.stdout.strip()) + "/resources/adapters.fa"

    return bbduk_adapter

def setup(datadir, sisrs_dir, sp):
    '''
    This function does the remaining folder setup that is needed to do
    trim reads.

    Arguments: 
    datadir (string): folder containing folders containing reads
    sisrs_dir (string): output directory
    sp (list): taxon 

    Returns: 
    raw read and trimmed data directories.
    '''

    #Set RawRead and TrimRead directories based off of script folder location
    raw_read_dir = datadir + '/' + sp
    trim_read_dir = sisrs_dir + "/Reads/TrimReads/" + sp

    #Create folder for BBDuk StdOut
    if not path.isdir(sisrs_dir + "/Reads/RawReads/trimOutput"):
        os.mkdir(sisrs_dir + "/Reads/RawReads/trimOutput")
    trim_output = sisrs_dir + "/Reads/RawReads/trimOutput"

    #Create folder for Raw FastQC output
    if not path.isdir(sisrs_dir + "/Reads/RawReads/fastqcOutput"):
        os.mkdir(sisrs_dir + "/Reads/RawReads/fastqcOutput")
    raw_fastqc_output = sisrs_dir + "/Reads/RawReads/fastqcOutput"

    #Create folder for Trimmed FastQC output
    if not path.isdir(trim_read_dir + "/fastqcOutput"):
        os.mkdir(trim_read_dir + "/fastqcOutput")
    trim_fastqc_output = trim_read_dir + "/fastqcOutput/"

    return raw_read_dir, trim_read_dir, trim_output, raw_fastqc_output, trim_fastqc_output

def newdFastqc(processors, fastqc_output, data_dir, newFiles):
    '''
    This function is designed to run the FastQC command with new data only. Modified
    to run both the raw and the trimmed data, using all available processors

    Arguments: the number of processors, the output directory, and the read directory.

    Returns: none.
    '''

    #Run FastQC on all trimmed files, using all available processors
    fastqc_command = [  #leave as list for subprocess.run, which treats list items as args
        'fastqc',
        '-t', f'{processors}',
        '-o', f'{fastqc_output}']

    for item in data_dir:
        for x in glob(item + "/*.fastq.gz"):
            f = x.split('/')

            if f[-1] in newFiles:
                fastqc_command.append(x)

    subprocess.run(fastqc_command, check = True)

def fastqcCommand(processors, fastqc_output, read_dir):
    '''
    Run FastQC on fastq.gz in specified taxon folders, using specified number of processors.

    Arguments: 
    processors (int): the number of processors
    fastqc_output (string): path to fastqc output directory
    read_dir (string): path to reads directory

    Returns: none.
    '''

    fastqc_command = [ #leave as list for subprocess.run, which treats list items as args
        'fastqc',
        '-t', f'{processors}',
        '-o', f'{fastqc_output}'
    ]

    print(read_dir + "/*.fastq.gz")
    flist = glob(read_dir + "/*.fastq.gz")
    if len(flist) == 0:
        sys.exit("No fastq.gz files found - please check your paths and extensions")
    print(flist)
    for x in flist:  #can add recursive here
        print(x)
        fastqc_command.append(x)

    subprocess.run(fastqc_command, check=True)

def trimHelper(tax_dir, trim_read_dir):
    '''
    This function obtains all of the
    single end and paired end read files in a folder.

    Arguments: 
    tax_dir (string): path to taxpm folder with fastq.gz files
    trim_read_dir (string): path to trimmed read directory

    Returns: 
    string: path to trimmed directory for the specific taxon 
    list: left paired end read paths (minus _1.fastq.gz)
    list: right paired end read paths (minus _2.fastq.gz) 
    list: single end read paths (minus .fastq.gz)
    '''

    #List all files and set output dir
    files = sorted(glob(tax_dir + "/*.fastq.gz"))
    taxon_ID = path.basename(tax_dir)
    out_trim_dir = trim_read_dir

    left_pairs = list()
    right_pairs = list()
    single_end = list()

    #Find files ending in _1/_2.fastq.gz
    left_files = [s for s in files if "_1.fastq.gz" in s]
    right_files = [s for s in files if "_2.fastq.gz" in s]

    lf = []
    rf = []

    #Strip _1.fastq.gz/_2.fastq.gz and identify pairs based on file name
    left_files = [x.replace('_1.fastq.gz', '') for x in left_files]
    right_files = [x.replace('_2.fastq.gz', '') for x in right_files]

    paired_files = list(set(left_files).intersection(right_files))

    #Reset file names and filter out single-end files
    for pair in paired_files:
        left_pairs.append(pair + "_1.fastq.gz")
        right_pairs.append(pair + "_2.fastq.gz")
    paired_files = sorted(left_pairs + right_pairs)

    single_end = [x for x in files if x not in paired_files]

    se = []

    #Remove .fastq.gz from lists to make naming easier
    left_pairs = [x.replace('_1.fastq.gz', '') for x in left_pairs]
    right_pairs = [x.replace('_2.fastq.gz', '') for x in right_pairs]
    single_end = [x.replace('.fastq.gz', '') for x in single_end]

    return out_trim_dir, left_pairs, right_pairs, single_end

def trim(raw_read_tax_dirs, trim_read_dir, bbduk_adapter, trim_output):
    '''
    Trim fastq files with cluster-aware memory management.
    '''
    
    # Memory monitoring configuration
    MEMORY_THRESHOLD_PERCENT = 85  # Don't let system exceed 85% RAM usage
    CHECK_INTERVAL = 2  # seconds between memory checks
    
    def get_available_memory_percent():
        """Get current system memory usage percentage"""
        return psutil.virtual_memory().percent
    
    def wait_for_memory(target_percent=70):
        """Wait until system memory drops below target"""
        while get_available_memory_percent() > target_percent:
            current = get_available_memory_percent()
            print(f"Memory usage at {current}%. Waiting for cleanup (target: {target_percent}%)...")
            time.sleep(CHECK_INTERVAL)
        print(f"Memory available. Proceeding with trim.")
    
    def run_trim_command(command, output_file, timeout=3600):
        """
        Execute trim command with streaming I/O and memory monitoring.
        """
        # Check memory before starting
        if get_available_memory_percent() > MEMORY_THRESHOLD_PERCENT:
            wait_for_memory(target_percent=70)
        
        try:
            print(f"Starting trim: {Path(output_file).name}")
            
            # Stream output with minimal buffering
            with open(output_file, 'w', buffering=1) as f:
                process = subprocess.Popen(
                    command,
                    stdout=f,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1,  # line buffering
                    env={**os.environ, 'OMP_NUM_THREADS': '1'}  # prevent nested parallelism
                )
                
                # Monitor process memory during execution
                process_start = time.time()
                while process.poll() is None:
                    if get_available_memory_percent() > MEMORY_THRESHOLD_PERCENT:
                        print(f"System memory exceeded {MEMORY_THRESHOLD_PERCENT}% during trim. Terminating process.")
                        process.terminate()
                        try:
                            process.wait(timeout=30)
                        except subprocess.TimeoutExpired:
                            process.kill()
                        raise RuntimeError("Trim terminated: system memory threshold exceeded")
                    time.sleep(CHECK_INTERVAL)
                
                # Final return code check
                if process.returncode != 0:
                    raise subprocess.CalledProcessError(
                        process.returncode, command
                    )
                
                elapsed = time.time() - process_start
                print(f"Trim completed in {elapsed:.1f}s: {Path(output_file).name}")
                
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Trim failed (exit code {e.returncode}): {output_file}")
        except Exception as e:
            raise RuntimeError(f"Trim failed: {str(e)}")
    
    def build_base_command(adapter_ref):
        """Build reusable base command parameters"""
        return {
            'maxns': '0',
            'ref': adapter_ref,
            'qtrim': 'w',
            'trimq': '15',
            'minlength': '35',
            'maq': '25',
            'k': '23',
            'mink': '11',
            'hdist': '1',
            'hdist2': '0',
            'ktrim': 'r',
            'ow': 't',
            'threads': '4'  # Limit threads per process to prevent memory spike
        }
    
    base_params = build_base_command(bbduk_adapter)
    
    print(f"Starting memory monitoring (threshold: {MEMORY_THRESHOLD_PERCENT}%)")
    
    tax_dir = raw_read_tax_dirs
    out_trim_dir, left_pairs, right_pairs, single_end = trimHelper(tax_dir, trim_read_dir)
        
    # Trim single-end files
    if single_end:
        print(f"Processing {len(single_end)} single-end files from {tax_dir}")
        for se_file in single_end:
            se_params = base_params.copy()
            se_params.update({
                'in': f'{se_file}.fastq.gz',
                'out': f'{out_trim_dir}/{Path(se_file).name}_Trim.fastq.gz'
            })
            
            se_command = ['bbduk.sh'] + [
                f'{k}={v}' for k, v in se_params.items()
            ]
            
            output_file = f'{trim_output}out_{Path(se_file).name}_Trim'
            run_trim_command(se_command, output_file)
            
            # Clean up after each file to free memory
            gc.collect()
        
    # Trim paired-end files
    if left_pairs and right_pairs and len(left_pairs) == len(right_pairs):
        print(f"Processing {len(left_pairs)} paired-end file pairs from {tax_dir}")
        for left, right in zip(left_pairs, right_pairs):
            pe_params = base_params.copy()
            left_name = Path(left).name
            right_name = Path(right).name
            
            pe_params.update({
                'in': f'{left}_1.fastq.gz',
                'in2': f'{right}_2.fastq.gz',
                'out': f'{out_trim_dir}/{left_name}_Trim_1.fastq.gz',
                'out2': f'{out_trim_dir}/{right_name}_Trim_2.fastq.gz'
            })
            
            pe_command = ['bbduk.sh'] + [
                f'{k}={v}' for k, v in pe_params.items()
            ]
            
            output_file = f'{trim_output}out_{left_name}_Trim'
            run_trim_command(pe_command, output_file)
            
            # Clean up after each file
            gc.collect()
    
    print("All trim operations completed successfully")

def makeLinks(data_path, sisrs_dir, sp):

    '''
    This function was designed to make soft links to all the files to avoid repeated trimming

    Arguments:
    data_path (string): the path to the folder with folders of trimmed reads
    sisrs_dir (string): output dir
    sp (str): taxon/folder name

    Returns: none.
    '''

    x = sp
    # Only the files that do not start with '.'
    data_path2 = data_path + "/Reads/TrimReads/" + x
    filestocopy = [f for f in listdir(data_path2) if f.endswith('.fastq.gz')]
    for f in filestocopy:
        # Creates the soft link to the files
        os.link(data_path2 + '/' + f,
                sisrs_dir+"/Reads/TrimReads/" + x + '/' + f)

def run2(datadir, outputdir, sp, processors, link):

    if link:
        makeLinks(datadir, outputdir, sp)
    else:
        bba = findAdapter()  # returns path to adapter file in bbduk
        raw_read_dir, trim_read_dir, trim_output, raw_fastqc_output, trim_fastqc_output = setup(datadir, outputdir, sp)  # list of the folders

        # fastqc raw data
        fastqcCommand(processors, raw_fastqc_output, raw_read_dir)  

        print(trim_read_dir)
        print(trim_output)
        trim(raw_read_dir, trim_read_dir, bba, trim_output)  # raw_read_tax_dirs, trim_read_dir, path to adapters, trim_output

        # fastqc trimmed data
        fastqcCommand(processors, trim_fastqc_output, trim_read_dir)  # trim_fastqc_output, trim_read_dir folders

if __name__ == "__main__":
    # Get arguments
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-d', '--datadir', action='store', nargs="?")
    my_parser.add_argument('-dir', '--outputdir', action='store', nargs="?")
    my_parser.add_argument('-sp', '--sp', action='store', nargs="?")
    my_parser.add_argument('-p', '--processors', action='store', default=1, nargs="?")
    my_parser.add_argument('--link', action=argparse.BooleanOptionalAction, default=False)
    args = my_parser.parse_args()

    run2(args.datadir, args.outputdir, args.sp, args.processors, args.link)
