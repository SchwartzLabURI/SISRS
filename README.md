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

* -dir --> Specify where the SISRS output will go. DEFAULT: The directory preceding where the scripts are located
    * This option is only available for the continuous run. The default is used in the individual scripts
* -id --> Specify the file that has all of the Taxon IDs. Providing this file requires the user to manually move the files
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
* -rd --> Specify the path where all the data is. Providing this path will make links to the correct locations for you
* -trm --> Tell's SISRS the data has been trimmed if it is present. It will link the files files in Reads/RawReads if this flag is not present and if it is present the links go in Reads/TrimReads. In the continues run it skips the trimming step. DEFAULT: untrimmed
* -th --> Specify the number of threads. DEFAULT: 1
* -gs --> Specify the approximate genome size estimate for group  
* -thrs --> Specify the threshold for the minimum site homozygosity for SISRS sites, must be between 0 and 1. DEFAULT: 1
* -mr --> Specify the minimum read coverage to call a SISRS site, must be 1-3. Default: 3
* -ms --> Specify the number of taxa to leave out. DEFAULT: 0

### REQUIREMENTS FOR CONTINUOUS SISRS RUN
In order to run SISRS you need to use the minimum of two flags, -id or -rd and -gs, without the use of these two flags the software will not run. All other flags will move to the default.

Data organization: Your data must be stored in a folder consisting of sub-directories that are labeled after the taxon names. Inside each folder there can only be files labeled (taxa)_1,_2,or _SE.fastq.gz.
Ex:
```
> cd SISRS/ ; ls
> README.md	SISRS_Small	TaxonIDs scripts
> cd SISRS_Small/ ; ls
> AotNan CalJac	ColAng	GorGor	HomSap	HylMol	MacMul	MacNem	PanPan	PanTro	PapAnu	PapCyn
> cd AotNan ; ls
> AotNan_1.fastq.gz	AotNan_2.fastq.gz	AotNan_SE.fastq.gz

```
### **IT IS IMPORTANT TO NOTE THAT THIS FILE NAMING CONVENTION IS NEEDED FOR ALL SISRS RUNS**

EXAMPLE COMMANDS:
```
> python sisrs.py -id TaxonIDs -gs 100000000 -th 20

> python sisrs.py -rd ./SISRS_Small/ -gs 100000000

> python sisrs.py -rd ./SISRS_Small/ -gs 100000000 -th 20 -trm

> python sisrs.py -id TaxonIDs -gs 100000000 -th 20 -thrs .66 -mr 2 -ms 1

```
