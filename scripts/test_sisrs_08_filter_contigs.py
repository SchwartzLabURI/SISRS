# run python3 -m pytest

import os
import pytest
import shutil

from sisrs_08_filter_contigs import *

@pytest.fixture
def setup_data():
    path_to_seq_lengths = '../test_data/contigs_SeqLength.tsv'
    path_to_sites = '../test_data/alignment_pi_locs_m4.txt'
    min_threshold = 1
    length_of_locus = 1000

    taxon_list = get_taxon_list('../test_data/TaxonList.txt')

    composite = SeqIO.to_dict(SeqIO.parse("../test_data/contigs.fa", "fasta"))

    high_count_contigs, contigcounts = get_high_count_contigs(path_to_seq_lengths, path_to_sites, min_threshold, length_of_locus)

    contigPath = '../test_data/'
    num_sp = 4
    new_contig_folder = '../test_data/temp/'
    if not path.exists(new_contig_folder):
        os.mkdir(new_contig_folder)
    else:
        shutil.rmtree(new_contig_folder)
        os.mkdir(new_contig_folder)

    contigs_remaining = write_alignment_plus_composite(high_count_contigs, contigPath, num_sp, composite, new_contig_folder)

    return(high_count_contigs, contigcounts, composite, taxon_list, contigs_remaining)


def test_get_high_count_contigs_basic(setup_data):
    assert len(setup_data[0]) == 2

def test_write_alignment_plus_composite_basic(setup_data):

    assert setup_data[4] == 2

def test_filter_contigs_distance_basic(setup_data):

    a = '../test_data/aligned/'

    contigs_remaining08 = filter_contigs_distance(setup_data[0], a, setup_data[3], 0.08, setup_data[1])
    contigs_remaining11 = filter_contigs_distance(setup_data[0], a, setup_data[3], 0.11, setup_data[1])

    assert len(contigs_remaining08) == 0
    assert len(contigs_remaining11) == 2