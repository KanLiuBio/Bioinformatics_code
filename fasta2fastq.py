#!/usr/bin/python
Usage = """
Convert the FASTA sequences to static FASTQ format
Quality values are encoded as H for 39
needs input file and ouputfile names

Usage: -version 1.0 (Python3)
  fasta2fastq.py inputfile.fasta output.fastq

Kan Liu
liukan.big@gmail.com
11/04/2018
"""

import sys
from Bio.SeqIO.FastaIO import SimpleFastaParser

if len(sys.argv)<3:
    print (Usage)
else:
   cmdargs = str(sys.argv)

   with open(str(sys.argv[1]),'r') as in_handle:
      with open(str(sys.argv[2]), "w") as out_handle:
         for title, seq in SimpleFastaParser(in_handle):
            out_handle.write("@%s\n%s\n+\n%s\n" \
                             % (title, seq, "H" * len(seq)))