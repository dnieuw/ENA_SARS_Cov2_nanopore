#!/usr/bin/env python
    
from Bio import SeqIO
from Bio.Align import PairwiseAligner
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-i',
                '--infile',
                metavar='File',
                help='Input file to be aligned',
                type=str,
                required=True)

parser.add_argument('-o',
                '--outfile',
                metavar='File',
                help='Output file of aligned sequence',
                type=str,
                required=True)

parser.add_argument('-r', 
                    '--reference',
                    help='Reference to be aligned to',
                    type=str,
                    required=True)

parser.add_argument('-n', 
                    '--seq_name',
                    help='Name of the aligned sequence',
                    type=str,
                    required=True)

args = parser.parse_args()

aligner = PairwiseAligner()
aligner.mode = 'global'
aligner.match_score = 1
aligner.mismatch_score = 0
aligner.open_gap_score = -2
aligner.extend_gap_score = -1

ref = SeqIO.read(args.reference,"fasta")
ref.seq = str(ref.seq.upper()).replace('-','N')
cons = SeqIO.read(args.infile,"fasta")
aln = aligner.align(ref.seq,cons.seq)
with open(args.outfile,'w') as out:
    print(">",args.seq_name, file=out)
    print(str(aln[0]).strip().split('\n')[2], file=out)