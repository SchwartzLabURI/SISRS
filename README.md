# **SISRS - PhyloSignal (Literman and Schwartz, 2019)**
##### This GitHub repo contains all the necessary scripts to run analyses from *McConnell, Literman, and Schwartz* (June 2019)

## Getting Started

### **Software Requirements** (Manuscript Version):
* Python3 (3.6.1)
  * Add packages  
* Bowtie2 (2.3.4)  
* Samtools (1.3.1)  
* BBMap Suite (37.40)  
* Genome assembler of your choosing (Ray - 2.3.2-devel)  
  *  Ray requires a working version of MPI, even if running on a single node
* R (3.5.3)  
  * Add packages  
* FastQC (0.11.5)
* RAxML (8.2.11)
* MAFFT (7.310)
* BEDTools (2.26.0)

### **Data Requirements:**
1) Approximate genome size estimate for group  
* Note: If genome size estimates vary greatly within your group, use an estimate on the larger side    

2) Raw or Trimmed NGS reads for each taxa  
* SISRS default parameters require 3X depth to call a site, so higher per-taxa coverage is ideal. In the associated manuscript, we required 20X coverage per taxa. The pipeline can be run with less, but site recovery will become  reduced as coverage drops.
* All read files should be of the same 'type' (e.g. Don't mix DNA-seq + RNA-seq)
* Ensure high sequence data quality prior to analysis (low read quality and high sequence duplication levels are both red flags for analysis). Scripts are provided to run FastQC.

3) Reference genome (.fasta) with associated annotations (.bed)

4) All RawReads must be labeled accordingly:
    1) (Taxon)_1/_2.fastq.gz for paired end reads
    2) (Taxon)_SE.fastq.gz for single end reads

5) Sample Data Set:
    1) SISRS_Small.zip
    2) It has the correct naming convention and folder setup

## **Instructions for Running SISRS**

There are two way to run SISRS

1) Run SISRS as a continuous program
    * Only script that needs to be ran directly: sisrs.py
    * WARNING: DEPENDING ON DATA SIZE THIS CAN TAKE A LONG AMOUNT OF TIME
2) Run SISRS as individual scripts
    * Scripts that need to be run directly
        * sisrs_01_folder_setup.py
        * sisrs_02_read_trimmer.py --> (OPTIONAL BASED ON DATA)
        * sisrs_03_read_subsetter.py
        * sisrs_04_ray_composite.py
        * sisrs_05_setup_sisrs.py
        * sisrs_06_run_sisrs.py
        * sisrs_07_output_sisrs.py

### Flags to run SISRS
#### Some flags can only be used for the continuous run of SISRS

* -dir; --directory --> Specify where the SISRS output will go. DEFAULT: The directory preceding where the scripts are located
    * This option is only available for the continuous run. The default is used in the individual scripts
* -id --> Specify the file that has all of the Taxon IDs. Providing this file requires the user to manually move the files.
    * This option is only available for the individual run
    * FILE MUST LOOK LIKE
```
> cat TaxonIDs

AotNan
CalJac
ColAng
GorGor
HomSap
HylMol
MacMul
MacNem
PanPan
PanTro
PapAnu
PapCyn
```
* -rd; --rawData --> Specify the path where all the data is. Providing this path will make links to the correct locations for you
* -trm; --trimmed --> Tell's SISRS the data has been trimmed if it is present. It will link the files files in Reads/RawReads if this flag is not present and if it is present the links go in Reads/TrimReads. In the continues run it skips the trimming step. DEFAULT: untrimmed
* -p; --processors --> Specify the number of threads. DEFAULT: 1
* -gs; --genomeSize --> Specify the approximate genome size estimate for group  
* -thresh; --threshold --> Specify the threshold for the minimum site homozygosity for SISRS sites, must be between 0 and 1. DEFAULT: 1
* -mr; --minread --> Specify the minimum read coverage to call a SISRS site, must be 1-3. Default: 3
* -ms; --missing --> Specify the number of taxa to leave out. DEFAULT: 0

### REQUIREMENTS FOR CONTINUOUS SISRS RUN
In order to run SISRS you need to use the minimum of two flags, -id or -rd and -gs, without the use of these two flags the software will not run. All other flags will move to the default.

Data organization:

Your data must be stored in a folder consisting of sub-directories that are labeled after the taxon names and nothing else. Inside each folder there can only be files labeled (taxa)_1,_2,or _SE.fastq.gz (*which is required for every version of SISRS*).

Ex:

```
> cd SISRS/ ; ls
> README.md	SISRS_Small	TaxonIDs scripts
> cd SISRS_Small/ ; ls
> AotNan CalJac	ColAng	GorGor	HomSap	HylMol	MacMul	MacNem	PanPan	PanTro	PapAnu	PapCyn
> cd AotNan ; ls
> AotNan_1.fastq.gz	AotNan_2.fastq.gz	AotNan_SE.fastq.gz
```

Example Commands:

```
> python sisrs.py -rd ./SISRS_Small/ -gs 100000000

> python sisrs.py -rd ./SISRS_Small/ -gs 100000000 -th 20 -trm

> python sisrs.py -rd ./SISRS_Small/ -gs 100000000 -th 20 -thrs .66 -mr 2 -ms 1
```

### REQUIREMENTS FOR SPLIT UP SISRS RUN

##### 1) Set up folder structure using sisrs_01_folder_setup.py

* In the main directory (e.g. /home/dmcconnell/SISRS) edit the included blank text file (TaxonIDs), adding each of your taxon IDs on a new line or place the directory with all of your data here organized as specified above in the Data organization section of the requirements for continuous run of SISRS.  
  * Note: Six-letter taxon IDs and folder names may be easier to work with (e.g. HomSap vs. Homo_sapiens)
* sisrs_01_folder_setup.py will create taxon folders in the RawReads, TrimReads, and SISRS_Run directories  

Flags needed: -rd and -trm is optional. It is important to note that if you use the -id flag you are expected to manually move your data into the correct folder.

Example:
```
> python scripts/sisrs_01_folder_setup.py -rd ./SISRS_Small/
> python scripts/sisrs_01_folder_setup.py -rd ./SISRS_Small/ -trm
```

##### 2) Trimming data with sisrs_02_read_trimmer.py

* Trimming can certainly be done using your preferred methods, but this script is a wrapper around the BBDuk program (part of the BBMap Suite).

If your reads are already trimmed, place all read files into the appropriate taxon folder in TrimReads.  

  * Multiple read files per taxon is fine, as are mixes of single-end and paired-end reads.

* 1: Put all raw reads in (**.fastq.gz format**) in the appropriate taxon folder in RawReads.  
    * Raw data files are not modified or used after trimming, so avoid duplication if possible by making links to the fastq.gz files as opposed to copying or moving raw data

* 2: Run sisrs_02_read_trimmer.py, which will:  
    * Trim all reads in RawReads
    * Output trimmed reads to the TrimReads directory  
    * Run FastQC on both trimmed and untrimmed datasets  

Flags (Optional): -th

Example:
```
python scripts/sisrs_02_read_trimmer.py -th <# processors>
```  

----
#### OPTIONAL MIDPOINT STEP (*THIS OPTION IS NOT EMBEDDED IN THE CONTINUOUS RUN*)
* Before continuing, **check FastQC output** for trimmed data. If using 'random' DNA-seq data, be especially wary of **high sequence duplication levels** which could indicate non-random sequencing.
* Data will eventually be pooled, so **best to remove low-quality data early** to prevent it from being incorporated into the genome assembly
----

##### 3) Subset trimmed samples with sisrs_03_read_subsetter.py

* Based on the requested depth [10*Genome Size Estimate/Number of Species], if a species has fewer total reads all reads will be sampled and the user will be notified with a warning message  
* For species that do have sufficient total coverage, an attempt is made to sample evenly from each read set.
* If certain read sets lack the required coverage, they are also sampled completely and the deficit is covered by a deeper read set
* Reads must be in the appropriate TrimReads subfolder by Taxon
* Single-end should all end in "_Trim.fastq.gz"
* Paired-end reads should end in "_Trim_1.fastq.gz/_Trim_2.fastq.gz"  

Flags needed: -gs

```
#For a SISRS run on ape species (with ~3.5Gb genome)
> python scripts/sisrs_03_read_subsetter.py -gs 3500000000
```

##### 4) Composite Genome Assembly with sisrs_04_ray_composite.py

* As with the trimming, you can use your preferred genome assembler, but we provide a script that wraps around Ray, which is fast but requires MPI even on one node

* SISRS identifies orthologs through a composite genome assembly step, where reads from all taxa are combined and assembled together  
* The subsampling step above ensures that no one species dominates the assembly, and also limits the assembly of loci that are present in few species   

Flags (Optioanl): -th

```
#1 node, 20 processors per node
> python scripts/sisrs_04_ray_composite.py -th <# processors>
```

##### 5) Setting up the SISRS run with sisrs_05_setup_sisrs.py

* This step will place all the files where they need to be for a SISRS run, including:
  * Renaming scaffolds with 'SISRS_' at the front for downstream data handling  
  * Linking all trim reads to the SISRS analysis folders  
  * Indexing and processing the composite genome  
  * Creating SISRS runs scripts for each species  

 Flags (Optional): -th, -mr, and -thrs

```
# 20 processors, 3 reads required to call a site, 100% homozygosity per species
> python scripts/sisrs_05_setup_sisrs.py -th 20 -mr 3 -thrs 1  

# 10 processors, 5 reads required to call a site, 99% homozygosity per species
> python scripts/sisrs_05_setup_sisrs.py
```

##### 6) Running Taxon-specific SISRS Steps with sisrs_06_run_sisrs.py

* sisrs_05_setup_sisrs.py generates a bash script in each taxon folder  
* These scripts are independent and can be run in parallel (e.g. on an HPC machine or cluster as separate jobs, but be sure to specify processors accordingly above)
* This script will run them serially on one machine (one after another)
* Individual log files are created in each taxon folder

Flags: NONE

```
> python scripts/sisrs_06_run_sisrs.py
```

##### 7) Complete SISRS Run and output alignments with sisrs_07_output_sisrs.py

* This script takes the taxon files and finds all variable sites  
* Sites are filtered into all sites, parsimony-informative sites, and biallelic sites
* This script also filters the biallelic site alignment down to only data with 0 missing taxa (or whatever number you choose), both with and without biallelic gaps
* This script also creates a final log with all mapping data and SISRS output in **SISRS_Run/out_SISRS_Log**

Flags (Optional): -ms

```
> python scripts/sisrs_07_output_sisrs.py -ms 3

> python scripts/sisrs_07_output_sisrs.py
```
