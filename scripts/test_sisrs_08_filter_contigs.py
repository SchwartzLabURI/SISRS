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

    composite = SeqIO.to_dict(SeqIO.parse("../test_data/contigs.fa", "fasta"))

    high_count_contigs, contigcounts = get_high_count_contigs(path_to_seq_lengths, path_to_sites, min_threshold, length_of_locus)

    return(high_count_contigs, contigcounts, composite)


def test_get_high_count_contigs_basic(setup_data):
    assert len(setup_data[0]) == 2

def test_write_alignment_plus_composite_basic(setup_data):
    high_count_contigs = setup_data[0]
    contigPath = '../test_data/'
    num_sp = 4
    composite = setup_data[2]

    new_contig_folder = '../test_data/temp/'
    if not path.exists(new_contig_folder):
        os.mkdir(new_contig_folder)
    else:
        shutil.rmtree(new_contig_folder)
        os.mkdir(new_contig_folder)

    contigs_remaining = write_alignment_plus_composite(high_count_contigs, contigPath, num_sp, composite, new_contig_folder)
    assert contigs_remaining == 2

    shutil.rmtree(new_contig_folder)