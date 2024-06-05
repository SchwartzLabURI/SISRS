from count_snps import *
import os

def test_count_snps():
    fa = '../test_data/contigs_for_probes.fa.fasta'
    aligned_contigs2 = '../test_data/aligned/'
    outf = '../test_data/variation_in_contigs.txt'
    num_pi_sites = count_snps(fa, aligned_contigs2, outf)

    assert num_pi_sites == 199

    os.remove(outf)