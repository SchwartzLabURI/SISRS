**Documentation**
=================

Table of Contents:
    * `Command Line Arguments`_
    * `About the Scripts`_
       * `sisrs_01_folder_setup`_
       * `sisrs_02_read_trimmer`_
       * `sisrs_03_read_subsetter`_
       * `sisrs_04_ray_composite`_
       * `sisrs_05_setup_sisrs`_
       * `sisrs_06_run_sisrs`_
       * `sisrs_07_output_sisrs`_


**********************
Command Line Arguments
**********************

1. Genome Size Estimate (-gs, --genomeSize)
    * Specify the approximate genome size estimate for group in basepairs (e.g. 3500000000 for primates)

2. Specifying starting data (-d, --data)
    * Path to directory where reads are already split into folders by taxon (these can be linked files).
    * **Note 1**: No spaces or special characters are allowed when naming taxon directories.
    * **Note 2**: If using -d option with pre-trimmed reads, you should also use the -trm flag, which tells SISRS to skip the trimming step (See SISRS_Small.zip for data structure)

3. Taxon ID File (-id)
    * Specify the taxon names with the TaxonIDs file. One taxon per line as follows:

        .. code-block:: bash

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

4. Output Directory (-dir,--directory)
    * Path to desired SISRS output directory. DEFAULT: The directory preceding where the scripts are located

5. Trimmed Data (-trm,--trimmed)
        * Specify if the data has already been trimmed
        * DEFAULT: untrimmed

6. Processors (-p,--processors)
    * Specify the number of available processors.
    * DEFAULT: 1
    * **Note**: If running this on a multi-core machine, specify the number of processors per node here.

7. Homozygosity Threshold (-thresh,--threshold)
    * Specify the minimum site homozygosity for SISRS sites, must be between 0 and 1.
    * DEFAULT: 1 (SISRS sites have support for only a single base within taxa)

8. Minimum Coverage Threshold (-mr,--minread)
    * Specify the minimum read coverage to call a SISRS site.
    * DEFAULT: 3 (Three reads required to call a site)

9. Missing Taxa Allowed (-ms,--missing)
    * When creating the final SISRS alignment, specify the maximum number of missing taxa allowed per column.
    * You can give it a single number or give it a range of numbers, such as 0-6, and it will do a final alignment for 0, 1, 2, 3, 4, 5, and 6 missing taxa allowed per column.
    *  It will also separate all the data out into folders labeled missing_(#).
    * DEFAULT: 0 (Coverage for all taxa for all sites)

10. Existing Run This feature will auto detect if a previous SISRS run has been done based on the file structure and if specific files are present.

    i. Adding a Taxon (-aT,--addTaxon): Adding a new taxon to the previous *SISRS* run
    ii. Adding Additional Sequences (-aD,--addData): Adding a new file to existing taxons

    **Note**: This relies on the -d/--data folder to specify were the new data is located


*****************
About the Scripts
*****************

This section will give a brief over view of what is needed to run each script. This will include the expected input (including flags) and it will also discuss the expected output.

**Note**: *SISRS* runs on an organized system in which you need to organize your data by taxon. Prior to running our software we ask that you make separate taxon folders for all taxa and move the data to the correct folder.

For examples please see the tutorial page.

sisrs_01_folder_setup
#####################

* This script will start the setup of the folder structure needed to run *SISRS*.
* **Required Arguments**: -d/--date or -id
* **Optional Arguments**: -trm/--trimmed

sisrs_02_read_trimmer
#####################

Running this script will:

    * Trim all reads in RawReads
    * Output trimmed reads to the TrimReads directory
    * Run FastQC on both trimmed and untrimmed datasets

* **Optional Arguments**: -p,--processors

Other Trimming Options:

    * Trimming can certainly be done using your preferred methods, but this script is a wrapper around the BBDuk program (part of the BBMap Suite).

    * If your reads are already trimmed, ensure they are in the TrimReads directory and move on to Step 3 (Subsetting)

    * Multiple read files per taxon is fine, as are mixes of single-end and paired-end reads.

Using our Wrapper:

    * Put all raw reads in (.fastq.gz format) in the appropriate taxon folder in RawReads.

    * Raw data files are not modified or used after trimming, so avoid duplication if possible by making links to the fastq.gz files as opposed to copying or moving raw data

    **Note**: Check FastQC output for trimmed data. If using 'random' DNA-seq data, be especially wary of high sequence duplication levels which could indicate non-random sequencing.

    * Data will eventually be pooled, so best to remove low-quality data early to prevent it from being incorporated into the genome assembly

sisrs_03_read_subsetter
#######################

Running this script will:

    * SISRS identifies orthologs through a composite genome assembly step, where reads from all taxa are combined and assembled together
    * The subsampling step ensures that no one species dominates the assembly, and also limits the assembly of loci that are present in few species
    * Based on the requested depth [10*Genome Size Estimate/Number of Species], if a species has fewer total reads all reads will be sampled and the user will be notified with a warning message
    * For species that do have sufficient total coverage, an attempt is made to sample evenly from each read set.
    * If certain read sets lack the required coverage, they are also sampled completely and the deficit is covered by a deeper read set

Expectations of data:

    * Single-end should all end in "_Trim.fastq.gz"
    * Paired-end reads should end in "_Trim_1.fastq.gz/_Trim_2.fastq.gz"
    * Reads must be in the appropriate TrimReads subfolder by Taxon

* **Required Arguments**: -gs,--genomeSize

sisrs_04_ray_composite
######################

This script provides:

    * A genome assembly script that wraps around Ray, which is fast but requires MPI even on one node. We plan to offer more assembly options in later releases. *mpirun must be in your path*

* **Optional Arguments**: -p,--processors

sisrs_05_setup_sisrs
####################

This script will:
    * Place all the files where they need to be for a *SISRS* run, including:
        * Renaming scaffolds with ``SISRS_`` at the front for downstream data handling
        * Linking all trim reads to the *SISRS* analysis folders
        * Indexing and processing the composite genome
        * Creating *SISRS* runs scripts for each species

* **Optional Arguments**: -p/--processors, -mr/--minread, and -thrs/--threshold

sisrs_06_run_sisrs
##################

This script will:

    * sisrs_05_setup_sisrs generates a bash script in each taxon folder
    * This script will run them serially on one machine (one after another)
    * If multiple nodes are available, you likely want to skip this step as these scripts are independent and can be run in parallel (e.g. on an HPC machine or cluster as separate jobs).
        * Be sure to specify processors accordingly when running sisrs_05_setup_sisrs
    * Individual log files are created in each taxon folder

* No Arguments needed

sisrs_07_output_sisrs
#####################

This script will:
    * Take the taxon files and find all variable sites
    * Sites are filtered into all sites, parsimony-informative sites, and biallelic sites
    * This script also filters the biallelic site alignment down to only data with 0 missing taxa (or whatever number or range you choose), both with and without biallelic gaps
    * This script also creates a final log with all mapping data and SISRS output in SISRS_Run/out_SISRS_Log

* **Optional Arguments**: -ms,--missing
