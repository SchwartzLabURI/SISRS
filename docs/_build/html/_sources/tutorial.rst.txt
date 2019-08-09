**SISRS2.0 Tutorial**
=====================

SISRS v2.0 can be run in two different ways:
    1. Continuous Run
        * This option may take a longer amount of time to run
    2. Split Run
        * This option allows the user to go through and check the output at each step

**Note**: All examples are working on the SISRS_Small.zip that was provided in the github. We are assuming that you have already done a git pull from `SISRS <github.com/SchwartzLabURI/SISRS>`_ and unzipped the small data set.

**********************
Continuous Run Example
**********************

This section will give you an example on how to run *SISRS* using the continues run option that is implemented.
As you move down in the list of commands more arguments get added to show more complex examples:

    .. code-block:: bash

       #Run SISRS with a genome size estimate of 100 MB
       > python sisrs.py -d ./SISRS_Small/ -gs 100000000
       > python sisrs.py --data ./SISRS_Small/ --genomeSize 100000000


       #Run SISRS with a genome size estimate of 2GB, 20 processors, and pre-trimmed reads
       > python sisrs.py -d ./SISRS_Small/ -gs 2000000000 -p 20 -trm
       > python sisrs.py --data ./SISRS_Small/ --genomeSize 2000000000 --processors 20 -trm


       #Run SISRS with a genome size estimate of 100bp, 10 processors, allowing 2/3 homozygosity,
       # a minimum read coverage of two reads, and allowing one taxon to be missing for any given
       # site in the final alignment
       > python sisrs.py -d ./SISRS_Small/ -gs 100 -p 10 -thresh .66 -mr 2 -ms 1
       > python sisrs.py --data ./SISRS_Small/ --genomeSize 100 --processors 10 --threshold .66 --minread 2 --missing 1

       #Run SISRS with a genome size estimate of 100bp, 10 processors, allowing 2/3 homozygosity,
       # a minimum read coverage of two reads, and allowing 1-6 taxon to be missing for any given
       # site in the final alignment
       > python sisrs.py -d ./SISRS_Small/ -gs 100 -p 10 -thresh .66 -mr 2 -ms 1-6
       > python sisrs.py --data ./SISRS_Small/ --genomeSize 100 --processors 10 --threshold .66 --minread 2 --missing 1-6


*****************
Split Run Example
*****************

The section is split up into subsections based on each of the seven scripts that defines the most important parts of the *SISRS* software. Each subsection will give you an example on how to run that individual script based on the some general command line options and some more complex options.

sisrs_01_folder_setup
#####################

    .. code-block:: bash

       > python scripts/sisrs_01_folder_setup.py -d /SISRS_Small/
       > python scripts/sisrs_01_folder_setup.py -data /SISRS_Small/

       #If data in already trimmed...
       > python scripts/sisrs_01_folder_setup.py -d /SISRS_Small/ -trm
       > python scripts/sisrs_01_folder_setup.py -data /SISRS_Small/ --trimmed

    You also have the option of using the -id flag which will require you yo manually move your data to the correct location after the script has completed running.

    .. code-block:: bash

       > python scripts/sisrs_01_folder_setup.py -id TaxonIDs

sisrs_02_read_trimmer
#####################

    .. code-block:: bash

       > python scripts/sisrs_02_read_trimmer.py -p <# processors>
       > python scripts/sisrs_02_read_trimmer.py --processors <# processors>

sisrs_03_read_subsetter
#######################

    .. code-block:: bash

       #For a SISRS run on ape species (with ~3.5Gb genome)
       > python scripts/sisrs_03_read_subsetter.py -gs 3500000000
       > python scripts/sisrs_03_read_subsetter.py --genomeSize 3500000000

sisrs_04_ray_composite
######################

    .. code-block:: bash

       #1 node, 20 processors per node
       > python scripts/sisrs_04_ray_composite.py -p 20
       > python scripts/sisrs_04_ray_composite.py --processors 20

       #3 nodes, 50 processors per node
       > python scripts/sisrs_04_ray_composite.py -p 150
       > python scripts/sisrs_04_ray_composite.py --processors 150

sisrs_05_setup_sisrs
####################

    .. code-block:: bash

       # 20 processors, 3 reads required to call a site, 100% homozygosity per species
       > python scripts/sisrs_05_setup_sisrs.py -p 20 -mr 3 -thrs 1

       # 10 processors, 5 reads required to call a site, 99% homozygosity per species
       > python scripts/sisrs_05_setup_sisrs.py --processors 10 --minread 5 --threshold .99

sisrs_06_run_sisrs
##################

    .. code-block:: bash

       > python scripts/sisrs_06_run_sisrs.py

sisrs_07_output_sisrs
#####################

    .. code-block:: bash

       #0 missing taxa allowed
       > python scripts/sisrs_07_output_sisrs.py

       #3 missing taxa allowed
       > python scripts/sisrs_07_output_sisrs.py -ms 3

       # 3, 4, 5, and 6 missing taxa allowed
       > python scripts/sisrs_07_output_sisrs.py --missing 3-6
