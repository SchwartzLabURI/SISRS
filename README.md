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
    * WARNING: DEPENDING ON DATA SIZE THIS CAN TAKE A LONG AMOUNT OF TIME
2) Run SISRS as individual scripts

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
* -trm --> Tell's SISRS the data has been trimmed if it is present. It will link the files files in Reads/RawReads if this flag is not present and if it is present the links go in Reads/TrimReads. In the continues run it skips the trimming step.
* -th --> Specify the number of threads. DEFAULT: 1
* -gs --> Specify the approximate genome size estimate for group  
* -thrs --> Specify the threshold for the minimum site homozygosity for SISRS sites. DEFAULT: 1
* -mr --> Specify the minimum read coverage to call a SISRS site. Default: 3
* -ms --> Specify the number of taxa to leave out. DEFAULT: 0
