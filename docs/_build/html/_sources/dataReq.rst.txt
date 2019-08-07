**Data Requirements**
=====================


1. You must have an approximate genome size estimate for group

    * Note: If genome size estimates vary greatly within your group, use an estimate on the larger side
2. Raw or Trimmed NGS reads for each taxa

    * Currently, starting NGS reads should either be all untrimmed or all trimmed.
    * Paired-end reads must have the same basename and end in _1/_2.fastq.gz to be recognized as paired end. Any other naming convention will trigger single-end assignment
    * SISRS default parameters require 3X depth to call a site, so higher per-taxa coverage is ideal. As a minimum we recommend at least 10X coverage. The pipeline can be run with less, but site recovery will become reduced as coverage drops.
    * All read files should be of the same 'type' (e.g. Don't mix DNA-seq + RNA-seq)
    * Ensure high sequence data quality prior to analysis (low read quality and high sequence duplication levels are both red flags for analysis). Note that the built-in trimming scripts are fairly conservative.
