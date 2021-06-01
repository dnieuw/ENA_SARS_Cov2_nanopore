#!/usr/bin/env python3
# Original author: David Newhouse
import sys, re, argparse, multiprocessing, itertools, textwrap
import pysam
from collections import Counter
from datetime import date

parser = argparse.ArgumentParser(description='Primer validation')

parser.add_argument('-b',
                    '--bam',
                    help="BAM file to turn into VCF",
                    type=str,
                    required = True)

parser.add_argument('-o',
                    '--out',
                    help='Output path for tab separated alignment result file',
                    type=str,
                   required = True)

parser.add_argument('-c',
                    '--cores',
                    help='Number of cores to use for processing',
                    default=1,
                    type=int,
                   required = False)

parser.add_argument('-d',
                    '--mindepth',
                    help='Minimal depth at which to not consider any alternative alleles',
                    default=10,
                    type=int,
                   required = False)

parser.add_argument('-af',
                    '--minAF',
                    help='Minimal allele frequency to output',
                    default=0.01,
                    type=float,
                   required = False)

parser.add_argument('-r',
                    '--reference',
                    help='Reference genome, must be the same as for BAM',
                    type=str,
                   required = True)

args = parser.parse_args()
insert_finder = re.compile("(.*)\+\d+(.*)")


def print_header(ref, outfile):
    file_date = date.today().strftime("%Y%m%d")
    columns = '\t'.join(["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"])
    print(textwrap.dedent(f"""\
        ##fileformat=VCFv4.0
        ##fileDate={file_date}
        ##source={' '.join(sys.argv)}
        ##reference={ref}
        ##contig={ref}
        ##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw Depth">
        ##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">
        ##INFO=<ID=DP4,Number=4,Type=Integer,Description="Counts for ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">
        ##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL">
        ##FILTER=<ID=minaf{args.minAF},Description="Allele frequency below indicated minimum">
        ##FILTER=<ID=mindp{args.mindepth},Description="Total depth below indicated minimum">
        #{columns}\
        """), file=outfile)


def work(start,stop):
    bamfile = pysam.AlignmentFile(args.bam, "rb")
    ref = bamfile.references[0]
    pileup = bamfile.pileup(ref, start=start, stop=stop, ignore_orphans=False, min_mapping_quality=0, min_base_quality=0)

    return([parse_column(p.reference_pos, p.get_query_sequences(add_indels=True),p.get_num_aligned()) for p in pileup])


def parse_column(ref_pos, allele_list, num_aln):
    if num_aln < args.mindepth:
        return(None)

    FWD_allele = Counter()
    REV_allele = Counter()
    COMB_allele = Counter()

    ref_allele = reference_seq[ref_pos]
    ref_FWD = 0
    ref_REV = 0

    def add_allele(ref, alt):
        if (ref.islower() | alt.islower()):
            REV_allele[(ref.upper(),alt.upper())] += 1
            COMB_allele[(ref.upper(),alt.upper())] += 1
        else:
            FWD_allele[(ref,alt)] += 1
            COMB_allele[(ref,alt)] += 1

    for var in allele_list:
        #Ignore positions that represent deletions
        if '*' in var:
            continue

        #Store the reference count as a special case
        if var.upper() == ref_allele:
            if var.isupper():
                ref_FWD += 1
            else:
                ref_REV += 1
        # - means next nucleotide is a deletion
        elif '-' in var:
            #Get the nucleotides that were deleted
            if var.islower():
                var=reference_seq[ref_pos:(ref_pos+len(var)-2)].lower()
            else:
                var=reference_seq[ref_pos:(ref_pos+len(var)-2)]
            add_allele(var, ref_allele)
        # + means next nucleotide is a insertion
        elif '+' in var:
            var = ''.join(insert_finder.match(var).groups())
            add_allele(ref_allele, var)
        else:
            add_allele(ref_allele, var)

    #Sort alleles by counts
    alleles = sorted(COMB_allele, key=COMB_allele.get, reverse=True)
    counts = [str(FWD_allele[k])+','+str(REV_allele[k]) for k in alleles]
    alt_freq = [(FWD_allele[k]+REV_allele[k])/num_aln for k in alleles]

    #Add ref allele in front or list
    alleles = [(ref_allele,ref_allele)]+alleles
    counts = [str(ref_FWD)+","+str(ref_REV)]+counts
    alt_freq = [(ref_FWD+ref_REV)/num_aln]+alt_freq

    return([str(num_aln), ref_pos, alleles, counts, alt_freq])


if __name__ == '__main__':
    bamfile = pysam.AlignmentFile(args.bam, "rb")
    ref = bamfile.references[0]
    reference_seq = pysam.Fastafile(args.reference).fetch(ref)

    #Make an array of start-stop intervals to parallelize processing
    genome_length = len(reference_seq)
    split = int(genome_length/args.cores)
    split_ranges = [[i*split+1,(i+1)*split] for i in range(args.cores)]
    #Adjust the last "stop" to be the genome length
    split_ranges[-1][1] = genome_length

    print("processing worklist")
    with multiprocessing.Pool(processes=args.cores) as p:
        resultlist = p.starmap(work, iter(split_ranges))

    print("writing file")
    with open(args.out,"w") as outfile:
        print_header(ref, outfile)
        for result in itertools.chain.from_iterable(resultlist):
            if result == None:
                continue
            num_aln, ref_pos, allele, counts, alt_freq = result
            for n,alt in enumerate(allele):
                #Do not print reference allele as a variant
                if alt[0]==alt[1]:
                    continue

                #Do not print minor alleles below 0.1%
                if alt_freq[n] < args.minAF:
                    continue

                #Print an INDEL flag if the alt is an indel
                indel = '' if (len(alt[0])+len(alt[1]))==2 else ';INDEL'

                #Print vcf lines, add 1 to the reference position to make the coordinates 1 based
                print(ref, ref_pos+1, ".", alt[0], alt[1], ".", "PASS", "DP="+num_aln+";AF="+"{:.6f}".format(alt_freq[n])+";DP4="+counts[0]+","+counts[n]+indel, sep='\t', file=outfile)

