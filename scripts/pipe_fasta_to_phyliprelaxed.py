#!/usr/bin/env python3
import sys
from Bio import SeqIO

SeqIO.convert(sys.stdin, "fasta", sys.argv[1], "phylip-relaxed")
