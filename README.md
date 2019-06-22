# **SISRS**
##### Dr. Rachel Schwartz, Dr. Robert Literman
##### Code Development: Devin McConnell  
##### University of Rhode Island

## Getting Started
**Note:** This repo contains a small sample set (SISRS_Small.zip). We provide this as a model for setting up data for a SISRS run, and also as a way to ensure all necessary software is installed prior to running with your data.

### **Software Requirements**:
* Python3 (3.6.4)
  * Biopython
  * Numpy
* Bowtie2 (2.3.4)  
* Samtools (1.3.1)  
* BBMap Suite (37.40)  
*  Ray - 2.3.2-devel  
  *  Ray requires a working version of MPI, even if running on a single node
* FastQC (0.11.5) (If trimming)
* RAxML (8.2.11)

### **Data Requirements:**
1) Approximate genome size estimate for group  
* Note: If genome size estimates vary greatly within your group, use an estimate on the larger side    

2) Raw or Trimmed NGS reads for each taxa  
* Currently, starting NGS reads should either be **all untrimmed or all trimmed**. We plan to implement mixed trim status in an upcoming release.  
* Paired-end reads must have the same basename and end in _1/_2.fastq.gz to be recognized as paired end. Any other naming convention will trigger single-end assignment
* SISRS default parameters require 3X depth to call a site, so higher per-taxa coverage is ideal. As a **minimum we recommend at least 10X coverage**. The pipeline can be run with less, but site recovery will become  reduced as coverage drops.
* All read files should be of the same 'type' (e.g. Don't mix DNA-seq + RNA-seq)
* Ensure high sequence data quality prior to analysis (low read quality and high sequence duplication levels are both red flags for analysis). Note that the built-in trimming scripts are fairly conservative.

## **Instructions for Running SISRS**

Currently there are two way to run SISRS:  

1) Run SISRS as a **continuous program**  
  * **WARNING:** Depending on the size of the dataset, running SISRS in serial mode may take days or weeks to run.

2) Run SISRS as **individual scripts** (allowing independent steps to be run on multiple cores if available, e.g. in an **HPC environment**)  

#### **Continuous Run Instruction**  
Continuous SISRS is run as
```python scripts/sisrs.py```, and requires at least two arguments:

1) **Genome Size Estimate** (-gs):  Specify the approximate genome size estimate for group in basepairs (e.g. 3500000000 for primates)

2) **Specifying starting data** (-rd): Path to directory where reads are already split into folders by taxon (these can be linked files).  
  - **Note 1:** No spaces or special characters are allowed when naming taxon directories.  
  -**Note 2:** If using -rd option with **pre-trimmed reads**, you should also use the **-trm** flag, which tells SISRS to skip the trimming step

**Optional Flags**  

**Output Directory** (-dir): Path to desired SISRS output directory. DEFAULT: The directory preceding where the scripts are located  

**Processors** (*-p*): Specify the number of available processors. DEFAULT: 1  
  - Note: If running this on a multi-core machine, specify the number of processors per node here.  

**Homozygosity Threshold** (-thresh): Specify the minimum site homozygosity for SISRS sites, must be between 0 and 1. DEFAULT: 1 (SISRS sites have support for only a single base within taxa)  

**Minimum Coverage Threshold** (-minread): Specify the minimum read coverage to call a SISRS site. DEFAULT: 3  (Three reads required to call a site)

**Missing Taxa Allowed** (-missing): When creating the final SISRS alignment, specify the maximum number of missing taxa allowed per column. DEFAULT: 0 (Coverage for all taxa for all sites)  

Example Commands:

```
#Run SISRS with a genome size estimate of 100 MB
> python sisrs.py -rd ./SISRS_Small/ -gs 100000000

#Run SISRS with a genome size estimate of 2GB, 20 processors, and pre-trimmed reads
> python sisrs.py -rd ./SISRS_Small/ -gs 2000000000 -p 20 -trm

#Run SISRS with a genome size estimate of 100bp, 10 processors, allowing 2/3 homozygosity, a minimum read coverage of two reads, and allowing one taxon to be missing for any given site in the final alignment
> python sisrs.py -rd ./SISRS_Small/ -gs 100 -p 10 -thresh .66 -minread 2 -missing 1
```

#### **Split Run Instruction**  

 1) Clone repo  
 2) Edit TaxonIDs file so that includes each Taxon ID on a new line. ex)  
 ```
 > cat TaxonIDs
 Homo_sapiens
 Pan_troglodytes
 Mus_musculus
 etc.
 ```
  * Note: Abbreviated taxon IDs may be easier to work with downstream (e.g. HomSap vs. Homo_sapiens)



3) When TaxonIDs file is ready, running ```python scripts/sisrs_01_folder_setup.py Taxon``` will create taxon folders in the RawReads, TrimReads, and SISRS_Run directories  

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
python scripts/rsisrs_02_read_trimmer.py -th <# processors>
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
