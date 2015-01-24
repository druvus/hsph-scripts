#!/usr/bin/env
from Bio import SeqIO
import sys

count = SeqIO.convert(sys.argv[1], "genbank", sys.argv[2], "fasta")