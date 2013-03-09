#!/usr/bin/env python

from Bio import AlignIO
import sys
from optparse import OptionParser
from pprint import pprint

# Arguments #
usage = """	%prog -i -o align_in align_out"""
description = """Description:	Use Biopython to convert the format of an alignment."""
epilog = """"""
version = "0.1"

parser = OptionParser(usage=usage, version = version, description=description, epilog=epilog)
parser.add_option("-i", "--input", dest="input", default="fasta", help="Input format of alignment [fasta]")
parser.add_option("-o", "--output", dest="output", default="phylip-relaxed", help="Output format of alignment [phylip-relaxed]")
parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False, help="verbose ouput")
(options, args) = parser.parse_args()

# I/O #
input_handle = open(sys.argv[1], "rb")
output_handle = open(sys.argv[2], "wb")
 
# conversion #
alignments = AlignIO.parse(input_handle, options.input)
AlignIO.write(alignments, output_handle, options.output)
 
output_handle.close()
input_handle.close()
