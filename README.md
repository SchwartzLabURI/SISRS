# **SISRS v2.0**
#### Dr. Rachel Schwartz, Dr. Robert Literman
##### Code Development: Devin McConnell  

### **Software Requirements**  
The following programs must be installed and in your path prior to running SISRS:
* Python3 (3.6.4)
  * Biopython
  * Numpy
* Bowtie2 (2.3.4)  
* Samtools (1.3.1)  
* BBMap Suite (37.40)  
*  Ray - 2.3.2-devel  
  *  Ray requires a working version of MPI, even if running on a single node
* FastQC (0.11.5) (Only if trimming)

### **Data Requirements:**
1) Approximate genome size estimate for group  
* Note: If genome size estimates vary greatly within your group, use an estimate on the larger side    

2) Raw or Trimmed NGS reads for each taxa  
* Currently, starting NGS reads should either be **all untrimmed or all trimmed**. We plan to implement mixed trim status in an upcoming release.  
* Paired-end reads must have the same basename and end in _1/_2.fastq.gz to be recognized as paired end. Any other naming convention will trigger single-end assignment
* SISRS default parameters require 3X depth to call a site, so higher per-taxa coverage is ideal. As a **minimum we recommend at least 10X coverage**. The pipeline can be run with less, but site recovery will become  reduced as coverage drops.
* All read files should be of the same 'type' (e.g. Don't mix DNA-seq + RNA-seq)
* Ensure high sequence data quality prior to analysis (low read quality and high sequence duplication levels are both red flags for analysis). Note that the built-in trimming scripts are fairly conservative.

### **Instructions for Running SISRS**

Currently there are two way to run SISRS:  

1) Run SISRS as a **continuous program**  
  * **WARNING:** Depending on the size of the dataset, running SISRS in serial mode may take days or weeks to run.

2) Run SISRS as **individual scripts** (allowing independent steps to be run on multiple cores if available, e.g. in an **HPC environment**)  

### **Continuous Run Instruction**  
Continuous SISRS is run as
```python scripts/sisrs.py```, and requires at least two arguments:  

**Required Arguments**  

**1) Genome Size Estimate** (*-gs, --genomeSize*):  Specify the approximate genome size estimate for group in basepairs (e.g. 3500000000 for primates)

**2) Specifying starting data** (*-rd, --rawData*): Path to directory where reads are already split into folders by taxon (these can be linked files).  
  - **Note 1:** No spaces or special characters are allowed when naming taxon directories.  
  - **Note 2:** If using -rd option with **pre-trimmed reads**, you should also use the **-trm** flag, which tells SISRS to skip the trimming step (See SISRS_Small.zip for data structure)  

**Optional Flags**  

**Output Directory** (*-dir,--directory*): Path to desired SISRS output directory. DEFAULT: The directory preceding where the scripts are located  

**Processors** (*-p,--processors*): Specify the number of available processors. DEFAULT: 1  
  - Note: If running this on a multi-core machine, specify the number of processors per node here.  

**Homozygosity Threshold** (*-thresh,--threshold*): Specify the minimum site homozygosity for SISRS sites, must be between 0 and 1. DEFAULT: 1 (SISRS sites have support for only a single base within taxa)  

**Minimum Coverage Threshold** (*-mr,--minread*): Specify the minimum read coverage to call a SISRS site. DEFAULT: 3  (Three reads required to call a site)

**Missing Taxa Allowed** (*-ms,--missing*): When creating the final SISRS alignment, specify the maximum number of missing taxa allowed per column. DEFAULT: 0 (Coverage for all taxa for all sites)  

Example Commands:

```
#Run SISRS with a genome size estimate of 100 MB
> python sisrs.py -rd ./SISRS_Small/ -gs 100000000
> python sisrs.py --rawData ./SISRS_Small/ -genomeSize 100000000


#Run SISRS with a genome size estimate of 2GB, 20 processors, and pre-trimmed reads
> python sisrs.py -rd ./SISRS_Small/ -gs 2000000000 -p 20 -trm
> python sisrs.py --rawData ./SISRS_Small/ -genomeSize 2000000000 --processors 20 -trm


#Run SISRS with a genome size estimate of 100bp, 10 processors, allowing 2/3 homozygosity, a minimum read coverage of two reads, and allowing one taxon to be missing for any given site in the final alignment
> python sisrs.py -rd ./SISRS_Small/ -gs 100 -p 10 -thresh .66 -mr 2 -ms 1
> python sisrs.py --rawData ./SISRS_Small/ --genomeSize 100 --processors 10 --threshold .66 --minread 2 --missing 1

```

### **Split Run Instruction**  

**1A) Specifying starting data with directory** (*-rd, --rawData*): Path to directory where reads are already split into folders by taxon (these can be linked files).  
  - **Note 1:** No spaces or special characters are allowed when naming taxon directories.  
  -**Note 2:** If using -rd option with **pre-trimmed reads**, you should also use the **-trm** flag, which tells SISRS to skip the trimming step (See SISRS_Small.zip for data structure)  

```
python scripts/sisrs_01_folder_setup.py -rd /SISRS_Small/

#If data in already trimmed...
python scripts/sisrs_01_folder_setup.py -rd /SISRS_Small/ -trm
```  

**1B) Specifying starting data with taxon list** (*-id*): Edit TaxonIDs file so that includes each Taxon ID on a new line. ex for SISRS_Small.zip)  
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
  * Note: Abbreviated taxon IDs may be easier to work with downstream (e.g. HomSap vs. Homo_sapiens)  

```python scripts/sisrs_01_folder_setup.py -id TaxonIDs``` will create taxon folders in the RawReads, TrimReads, and SISRS_Run directories  

**Note:** It is important to note that if you use the -id flag you are expected to manually move your data into the correct folder (Raw vs. Trimmed).  

----

**2) Trimming Data**  


Running ```python scripts/sisrs_02_read_trimmer.py``` will:  

  * Trim all reads in RawReads
  * Output trimmed reads to the TrimReads directory  
  * Run FastQC on both trimmed and untrimmed datasets  

Flags (Optional): -p,--processors

Example:
```
python scripts/sisrs_02_read_trimmer.py -p <# processors>
```  
* Trimming can certainly be done using your preferred methods, but this script is a wrapper around the BBDuk program (part of the BBMap Suite).

* **If your reads are already trimmed**, ensure they are in the TrimReads directory and move on to Step 3 (Subsetting)

* Multiple read files per taxon is fine, as are mixes of single-end and paired-end reads.

* Put all raw reads in (**.fastq.gz format**) in the appropriate taxon folder in RawReads.  
    * Raw data files are not modified or used after trimming, so avoid duplication if possible by making links to the fastq.gz files as opposed to copying or moving raw data  

**Note: Check FastQC output** for trimmed data. If using 'random' DNA-seq data, be especially wary of **high sequence duplication levels** which could indicate non-random sequencing.
* Data will eventually be pooled, so **best to remove low-quality data early** to prevent it from being incorporated into the genome assembly  

----

##### 3) Subset trimmed samples with sisrs_03_read_subsetter.py
* SISRS identifies orthologs through a composite genome assembly step, where reads from all taxa are combined and assembled together  
* The subsampling step ensures that no one species dominates the assembly, and also limits the assembly of loci that are present in few species   
* Single-end should all end in "_Trim.fastq.gz"
* Paired-end reads should end in "_Trim_1.fastq.gz/_Trim_2.fastq.gz"  
* Reads must be in the appropriate TrimReads subfolder by Taxon  

Flags needed: -gs,--genomeSize

```
#For a SISRS run on ape species (with ~3.5Gb genome)
> python scripts/sisrs_03_read_subsetter.py -gs 3500000000
```

* Based on the requested depth [10*Genome Size Estimate/Number of Species], if a species has fewer total reads all reads will be sampled and the user will be notified with a warning message  
* For species that do have sufficient total coverage, an attempt is made to sample evenly from each read set.
* If certain read sets lack the required coverage, they are also sampled completely and the deficit is covered by a deeper read set


----

##### 4) Composite Genome Assembly with sisrs_04_ray_composite.py

* We provide a genome assembly script that wraps around Ray, which is fast but requires MPI even on one node. We plan to offer more assembly options in later releases.
*```mpirun``` must be in your path

Flags (Optional): -p,--processors

```
#1 node, 20 processors per node
> python scripts/sisrs_04_ray_composite.py -p 20

#3 nodes, 50 processors per node
> python scripts/sisrs_04_ray_composite.py -p 150
```  

----

##### 5) Setting up the SISRS run with sisrs_05_setup_sisrs.py

* This step will place all the files where they need to be for a SISRS run, including:
  * Renaming scaffolds with 'SISRS_' at the front for downstream data handling  
  * Linking all trim reads to the SISRS analysis folders  
  * Indexing and processing the composite genome  
  * Creating SISRS runs scripts for each species  

 Flags (Optional): -p/--processors, -mr/--minread, and -thrs/--threshold

```
# 20 processors, 3 reads required to call a site, 100% homozygosity per species
> python scripts/sisrs_05_setup_sisrs.py -p 20 -mr 3 -thrs 1  

# 10 processors, 5 reads required to call a site, 99% homozygosity per species
> python scripts/sisrs_05_setup_sisrs.py --processors 10 --minread 5 --threshold .99
```  

----

##### 6) Running Taxon-specific SISRS Steps with sisrs_06_run_sisrs.py

* sisrs_05_setup_sisrs.py generates a bash script in each taxon folder  
* This script will run them serially on one machine (one after another)

* **If multiple nodes are available, you likely want to skip this step** as these scripts are independent and can be run in parallel (e.g. on an HPC machine or cluster as separate jobs). Be sure to specify processors accordingly above.
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

Flags (Optional): -ms,--missing

```
> python scripts/sisrs_07_output_sisrs.py        #0 missing taxa allowed
> python scripts/sisrs_07_output_sisrs.py -ms 3  #3 missing taxa allowed

```
