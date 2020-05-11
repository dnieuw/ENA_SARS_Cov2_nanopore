#!/usr/bin/env python 

import pysam
import sys
from collections import Counter
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-i',
                '--infile',
                metavar='File',
                help='Input bamfile; should be sorted and indexed',
                type=str,
                required=True)

parser.add_argument('-o',
                '--outfile',
                metavar='File',
                help='Output fasta file for consensus (default stdout)',
                default=sys.stdout,
                type=argparse.FileType('w'),
                required=False)

parser.add_argument('-d', 
                    '--mindepth',
                    help='Minimal depth needed (default: 1)',
                    default=1,
                    type=int,
                    required=False)

parser.add_argument('-g', 
                    '--keepgap',
                    help='Keep gaps in the alignment as "N" regions',
                    default=0,
                    type=int,
                    required=False)

parser.add_argument('-l', 
                    '--keepdel',
                    help='Keep deletions in the reads as "-" (adhere more to the reference sequence)',
                    default=0,
                    type=int,
                    required=False)

args = parser.parse_args()

def makeConsensus(bamfile, covlim):
    with pysam.AlignmentFile(bamfile, "rb") as bamfile:
        all_ref = []
        for ref in bamfile.references:
            consensus = []
            prevpos = -1
            for n, pileupcolumn in enumerate(bamfile.pileup(ref, ignore_orphans=False, min_mapping_quality=0, min_base_quality=0)):
                #Keep track of gaps, and different starting position compared to reference
                if (pileupcolumn.pos != prevpos+1 and args.keepgap == 1):
                    gapsize = pileupcolumn.pos - prevpos - 1
                    consensus.append('N'*gapsize)
                prevpos=pileupcolumn.pos
                
                cov = pileupcolumn.get_num_aligned()
                if (cov < covlim):
                    consensus.append('N')
                    continue
                pos = Counter()
                for pileupread in pileupcolumn.pileups:
                    if pileupread.is_del:
                        pos['*'] += 1
                    elif not pileupread.is_refskip:
                        nt = pileupread.alignment.query_sequence[pileupread.query_position]
                        pos[nt] += 1
                best = max(pos, key=pos.get)
                if best=='*':
                    #frac_del = pos['*']/sum(pos.values())
                    #print("Warning, deletion at: "+str(n)+"; Fraction of reads with deletion: "+str(frac_del))
                    if args.keepdel == 1:
                        consensus.append("-") #If there is a deletion in the mapping we print an '-'
                        continue
                    else:
                        continue #If there is a deletion in the mapping we continue without printing a nucleotide
                consensus.append(best)
            if len(consensus) == 0:
                continue
            all_ref.append('>'+ref+'_consensus\n'+''.join(consensus))
        return('\n'.join(all_ref))

if __name__ == '__main__':
    consensus = makeConsensus(args.infile, args.mindepth)
    with args.outfile as out_consensus:
        print(consensus, file=out_consensus)