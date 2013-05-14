#!/usr/bin/env python

#from BCBio import GFF
import GFF
from Bio import SeqIO
import sys

# Arguments #
usage = """     %prog genbank_in_file gff_out_name"""
description = """Description:   Use Biopython to convert a genbank to a gff."""
epilog = """"""
version = "0.1"


in_file = sys.argv[1]
out_file = sys.argv[2]
in_handle = open(in_file)
out_handle = open(out_file, "w")
GFF.write(SeqIO.parse(in_handle, "genbank"), out_handle)
in_handle.close()
out_handle.close()
