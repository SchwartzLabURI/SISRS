## **PhyloSignal (Literman and Schwartz, 2019)**  
#### This GitHub repo contains all the necessary scripts to run analyses from *Literman and Schwartz* (2019)

#### **Software Requirements** (Manuscript Version):
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

#### **Data Requirements:**
1) Approximate genome size estimate for group  
* Note: If genome size estimates vary greatly within your group, use an estimate on the larger side    

2) Raw or Trimmed NGS reads for each taxa  
* SISRS default parameters require 3X depth to call a site, so higher per-taxa coverage is ideal. In the associated manuscript, we required 20X coverage per taxa. The pipeline can be run with less, but site recovery will become  reduced as coverage drops.
* All read files should be of the same 'type' (e.g. Don't mix DNA-seq + RNA-seq)
* Ensure high sequence data quality prior to analysis (low read quality and high sequence duplication levels are both red flags for analysis). Scripts are provided to run FastQC.


3) Reference genome (.fasta) with associated annotations (.bed)

#### **Instructions for Running Scripts**

##### 1) Clone Repo locally
```
cd /home
git clone https://github.com/BobLiterman/Literman_PhyloSignal.git
```
* Repo contains:  
  * All necessary scripts
  * Base folder architecture for a SISRS run and subsequent analyses
  * Empty Taxon ID file

##### 2) Set up folder structure using folder_setup.py
* In the main directory (e.g. /home/Literman_PhyloSignal/) edit the included blank text file (TaxonIDs), adding each of your taxon IDs on a new line.  
  * Note: Six-letter taxon IDs may be easier to work with (e.g. HomSap vs. Homo_sapiens)
* folder_setup.py will create taxon folders in the RawReads, TrimReads, and SISRS_Run directories  

```
> cd /home/Literman_PhyloSignal
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

> python scripts/folder_setup.py TaxonIDs
```

##### 3) If starting from raw (untrimmed) NGS reads, and trimming with our wrapper...  
Trimming can certainly be done using your preferred methods, but we provide a wrapper script for the BBDuk program (part of the BBMap Suite).  

If your reads are already trimmed, place all read files into the appropriate taxon folder in TrimReads.  

  * Multiple read files per taxon is fine, as are mixes of single-end and paired-end reads.

To use our wrapper:

* 1) Put all raw reads in (**.fastq.gz format**) in the appropriate taxon folder in RawReads.  
  * Raw data files are not modified or used after trimming, so avoid duplication if possible by making links to the fastq.gz files as opposed to copying or moving raw data  

```
> cd /home/Literman_PhyloSignal/Reads/RawReads/HomSap
> cp -as /Path/to/your/raw/data/HomSap*.gz .
> cd ../
> cd GorGor/
> cp -as /Path/to/your/raw/data/GorGor*.gz .
```

* 2) Run read_trimmer.py, which will:  
  * Trim all reads in RawReads
  * Output trimmed reads to the TrimReads directory  
  * Run FastQC on both trimmed and untrimmed datasets  

```
python scripts/read_trimmer.py <# processors>
```  

##### !!! Data Check !!!
* Before continuing, **check FastQC output** for trimmed data. If using 'random' DNA-seq data, be especially wary of **high sequence duplication levels** which could indicate non-random sequencing.
* Data will eventually be pooled, so **best to remove low-quality data early** to prevent it from being incorporated into the genome assembly

##### 4) Once all trimmed data is in place and QC'd, samples are subset with the following scheme:  
* read_subsetter.py takes one argument: An estimate of the group genome size in basepairs  

```
#For a SISRS run on ape species (with ~3.5Gb genome)
> python scripts/read_subsetter.py 3500000000
```
* Based on requested depth [10*Genome Size Estimate/Number of Species], if a species has fewer total reads all reads will be sampled and the user will be notified with a warning message  
* For species that do have sufficient total coverage, an attempt is made to sample evenly from each read set.
* If certain read sets lack the required coverage, they are also sampled completely and the deficit is covered by a deeper read set
* Reads must be in the appropriate TrimReads subfolder
* Single-end should all end in "_Trim.fastq.gz"
* Paired-end reads should end in "_Trim_1.fastq.gz/_Trim_2.fastq.gz"  

##### 5) Composite Genome Assembly  
* SISRS identifies orthologs through a composite genome assembly step, where reads from all taxa are combined and assembled together  
* The subsampling step above ensures that no one species dominates the assembly, and also limits the assembly of loci that are present in only few species  
* As with the trimming, you can use any genome assembler that you're comfortable with, but we provide  a wrapper for Ray, which is fast but requires MPI even on one node  
* **ray_composite.py** takes two arguments: (1) Number of nodes available (2) Number of processors per node  
```
#1 node, 20 processors per node
> python scripts/ray_composite.py 1 20
```

##### 6) Setting up the SISRS run  
* **setup_sisrs.py** will place all the files where they need to be for a SISRS run, including:
  * Renaming scaffolds with 'SISRS_' at the front for downstream data handling  
  * Linking all trim reads to the SISRS analysis folders  
  * Indexing and processing the composite genome  
  * Creating SISRS runs scripts for each species  
  * The bulk of the user-settings take place here including arguments in this order:  
    1) Number of available processors  
    2) Minimum read coverage to call a SISRS site (Default: 3)   
    3) Minimum site homozygosity for SISRS sites (Default: 1)  

```
# 20 processors, 3 reads required to call a site, 100% homozygosity per species
> python scripts/setup_sisrs.py 20 3 1  

# 10 processors, 5 reads required to call a site, 99% homozygosity per species
> python scripts/setup_sisrs.py 10 5 .99
```

##### 6) Running Taxon-specific SISRS Steps
* setup_sisrs.py generates a bash script in each taxon folder  
* These scripts are independent and can be run in parallel (e.g. on an HPC machine or cluster as separate jobs, but be sure to specify processors accordingly above)  
* To run them serially on one machine (one after another), you can run **run_sisrs.py**  
* Individual log files are created in each taxon folder
```
> python scripts/run_sisrs.py
```

##### 7) Complete SISRS Run and output alignments
* **output_sisrs.py** takes the taxon files and finds all variable sites  
* Sites are filtered into all sites, parsimony-informative sites, and biallelic sites
* This script also filters the biallelic site alignment down to only data with 0 missing taxa (or whatever number you choose), both with and without biallelic gaps
* This script also creates a final log with all mapping data and SISRS output in **SISRS_Run/out_SISRS_Log**
```
> python scripts/output_sisrs.py 0
```
