.. SISRS2.0 documentation master file, created by
   sphinx-quickstart on Mon Aug  5 12:10:22 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SISRS: SNP Identification from Short Read Sequences
===================================================

.. figure:: sisrs.png
   :width: 500px
   :align: center
   :height: 250px
   :alt: alternate text
   :figclass: align-center

   *SISRS* — pronounced “scissors” — is a program for identifying phylogenetically informative sites from next-generation whole-genome sequencing of multiple species. It identifies homologous sites without the need to do de novo assembly, annotation, and alignment. It identifies conserved regions by doing joint de novo assembly on multiple species. Sequencing reads are then aligned back to the contigs to identify variable sites.

**Authors**:
    This software is written by member of:
        `Schwartz Lab @ URI <schwartzlaburi.github.io/index.html>`_
    Developers:
        * Dr. Rachel Schwartz (Original developer)
        * Dr. Robert Literman
        * Devin J. McConnell

**Associated Publications**:
    1. `Literman, R. and R.S. Schwartz. 2019. Genome-scale profiling reveals higher proportions of phylogenetic signal in non-coding data. BioRxiv Preprint. <www.biorxiv.org/content/10.1101/712646v2>`_
    2. `Harkins, K.M., R.S. Schwartz, R.A. Cartwright, and A.C. Stone. 2016. Phylogenomic reconstruction supports supercontinent origins for Leishmania. Infection, Genetics and Evolution 38: 101-109. <www.biorxiv.org/content/10.1101/028969v1.full>`_
    3. `Schwartz, R.S., K.M Harkins, A.C. Stone, and R.A. Cartwright. 2015. A composite genome approach to identify phylogenetically informative data from next-generation sequencing. BMC Bioinformatics. 16:193. <bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0632-y>`_

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   softwareSetup
   dataReq
   tutorial

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
