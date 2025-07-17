#!/usr/bin/env python3

import os
import argparse
import VCF

def get_taxon_list(tpath):
    '''

    Args:
        tpath (string): path to the file that contains the liast of species

    Returns:
        list: list of species (folders)

    '''
    with open(tpath) as f:
        taxon_list = f.readlines()
        taxon_list = sorted([x.rstrip() for x in taxon_list])
    return (taxon_list)

def format_consensus_output(output_path, taxa_threshold, full_seqs):
    '''
    This function outputs each contig to a file (note N's have been removed so needs aligning).

    Arguments: path to the output directory, taxa threshold, dict of contigs.

    Returns: none.
    '''

    #where to put contigs
    p = output_path + 'SISRS_Run/contigs_outputs/'
    if not os.path.exists(p):
        os.makedirs(p)

    for contig, seq_dict in full_seqs.items():  #full_seqs[contig][taxon] = seq string
        if len(seq_dict) >= taxa_threshold:  # check how many taxa passed
            # if passing the threshold, create a file
            contigs_file = output_path + 'SISRS_Run/contigs_outputs/' + contig + '.fasta'
            out_file = open(contigs_file,'w')

            # populate the file
            for dir in sorted(seq_dict):
                out_file.write(">" + dir + '\n')
                out_file.write(seq_dict[dir] + '\n')

            # close the output handle
            out_file.close()

def read_vcf_gz(file_path):
    """
    Process the samples from a .vcf.gz file and print information about each variant's samples.

    This function opens a gzipped VCF file, iterates over each variant record, and prints
    detailed sample-specific information including genotype, depth of coverage, and allelic depth
    for each sample.

    Assumes that ref = N and only one sample.

    Parameters:
    file_path (str): The path to the gzipped VCF file. This file should be in .vcf.gz format.

    Returns:
    Dict: This function returns a dict where keys are contigs and values are dicts where keys are pos and values are
        - Alleles (list)
        - Allelic depth (D) (list) - provides counts for reference and alternative alleles' supporting reads
    """

    # Nested dictionary to store VCF data by chromosome and position
    vcf_data = {}

    vcf_reader = VCF(file_path)
    for record in vcf_reader:
        alt = [str(a) for a in record.ALT if str(a) != '.']
        if alt:  # Skip records where ALT is empty after filtering out '.'
            chrom = record.CHROM[13:]  # contig number
            pos = record.POS  # Position
            allelic_depth = str(record).strip().split(':')[-1].split(',')[1:] #get the line as a string and manually split out AD
            # Initialize chromosome key if not already present
            if chrom not in vcf_data:
                vcf_data[chrom] = {}
            # Store information in the nested dictionary structure
            vcf_data[chrom][pos] = {
                'A': alt,
                'D': allelic_depth
            }
    return vcf_data

def vcf_consensus(output_path, coverage_threshold, hz_threshold):
    '''
    This function produces the consensus sequences from VCF data.
    The function assumes all sites (including non variable) were
    output.
    First, get contig names from the header, this is needed since
    even with all sites output, loci without any reads are completely
    missing from the VCF body.
    Then iterate over VCF body to get sites, keep track of current
    contig in the body vs in the header, fill in empty contigs (loci),
    for non empty, filter sites with low allelic coverage and
    positions marked with '.', as well as filter out (replace with N)
    sequences with high heterozygosity.

    Arguments:
        output_path (string): path to the output directory
        coverage_threshold (int): coverage threshold
        heterozygosity threshold (float)

    Returns: none
    '''

    # list of taxa
    taxon_list = get_taxon_list(output_path + 'TaxonList.txt')

    full_seqs = {}
    for directory in taxon_list:
        taxon_vcf_file = f"{output_path}SISRS_Run/{directory}/{directory}.vcf.gz" #/work/pi_rsschwartz_uri_edu/jmartins/solanum/SISRS_Run
        vcf_dict = read_vcf_gz(taxon_vcf_file)

        taxon_consensus_file = f"{output_path}SISRS_Run/{directory}/{directory}_single_line_format_consensus.fa"
        if os.path.exists(taxon_consensus_file): # in cases where the analysis ran previously and we are linking, this file would already exist but permissions may be wrong to overwrite
            os.remove(taxon_consensus_file)

        outfile = open(taxon_consensus_file, "w")

        #data on heterozygosity - for debugging purposes
        hzfile = output_path + 'SISRS_Run/' + directory + '/' + directory +'_hztable.csv'
        if os.path.exists(hzfile): # in cases where the analysis ran previously and we are linking, this file would already exist but permissions may be wrong to overwrite
            os.remove(hzfile)
        hzhandle = open(hzfile, 'w')

        #iterate over each site in vcf
        for contig, contig_info in vcf_dict.items():
            loclen = len(contig_info) #number of pos - ie length of contig (excluding N)
            print(loclen)
            hznum = 0 #track hets per locus
            outputseq = ""

            for pos, details in contig_info.items(): #each contig has a dict where the pos is a key and details dict: {'A': ['T'], 'D': [7, 4]}
                if len(details['A']) > 1:
                    hznum += 1
                max_value = max(details['D']) #allelic depth of most common allele
                if max_value >= coverage_threshold: #sufficient coverage for that allele
                    alt_index = details['D'].index(max_value)
                    outputseq += details['A'][alt_index]

            hzval = hznum / loclen #check for too many heterozygotes in the locus which could mean paralogy
            print(contig + "," + str(hzval), file=hzhandle)
            if hzval <= hz_threshold:  #include locus seq in the output if it's reasonable
                print(">" + contig, file=outfile)
                print(outputseq, file=outfile)
                if contig not in full_seqs:
                    full_seqs[contig] = {}
                full_seqs[contig][directory] = outputseq

        #close handles
        outfile.close()
        hzhandle.close()

    return full_seqs

def contig_driver(output_path, proc):
    '''
    This function perform pileups on each of the bam files and then outputs a vcf.gz based on bcftools call.

    Arguments: path to the output directory.

    Returns: none.

    '''

    # list of taxa
    taxon_list = get_taxon_list(output_path + '/TaxonList.txt')

    for tax in taxon_list:
        bam_file = output_path+'SISRS_Run/'+tax+'/'+tax+'.bam'
        vcf_zipped = output_path+'SISRS_Run/'+tax+'/'+tax+'.vcf.gz'

        if os.path.exists(vcf_zipped): # in cases where the analysis ran previously and we are linking, this file would already exist but permissions may be wrong to overwrite
            os.remove(vcf_zipped)

        #do pileup using no reference and report allelic depth
        #then call variants using all sites (M) and zip the output
        bcf_command = f"bcftools mpileup -Ou --no-reference -a FORMAT/AD {bam_file} | bcftools call -Oz --threads {proc} -mM -o {vcf_zipped}"
        os.system(bcf_command)

def run7a(output_path, taxa_threshold, cov_thresh, hz_thresh, proc):
    contig_driver(output_path, proc)

    #call consensus and filter by allelic coverage and ratio of heterozygous sites
    full_seqs = vcf_consensus(output_path, coverage_threshold=cov_thresh, hz_threshold=hz_thresh)

    #recompile sequence of each taxon by locus
    format_consensus_output(output_path, taxa_threshold, full_seqs)

if __name__ == '__main__':

    # python3 sisrs_07_a_contigs_processing.py -t $T -dir $OUTFOLDER
    # Get arguments
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-t', '--threshold', action='store', nargs="?", type=int)
    my_parser.add_argument('-d', '--dir', action='store', nargs="?")
    my_parser.add_argument('-p', '--processors', type=int, action='store', default=1, nargs="?")
    my_parser.add_argument('-c', '--cov', action='store', default=3, nargs="?", type=int)
    my_parser.add_argument('-z', '--hz', action='store', default=0.01, nargs="?", type=float)

    args = my_parser.parse_args()

    run7a(args.dir + '/', args.threshold, args.cov, args.hz, args.processors)
