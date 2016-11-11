#!/usr/bin/env python3
import pysam
import screed
import khmer_utils as ku
from sys import argv, stderr, stdout


def printread(read, file=None):
    if read.quality:
        print("@", read.name, ' ', read.description, sep='', file=file)
        print(read.sequence, file=file)
        print("+", file=file)
        print(read.quality, file=file)
    else:
        print(">", read.name, ' ', read.description, sep='', file=file)
        print(read.sequence, file=file)


def notag(readid):
    readid, _, _ = readid.partition(' ')
    if readid[-2] == '/' and readid[-1] in "12":
        return readid[:-2]
    return readid


def process_file(bamfile, regions, readfile, outstream=stdout):
    samf = pysam.AlignmentFile(bamfile, 'rb')
    to_extract = set()
    for region in regions:
        for read in samf.fetch(region):
            name = read.query_name
            if not read.is_secondary and not read.is_supplementary:
                to_extract.add(name)
    print("Collected", len(to_extract), "read names to extract", file=stderr)
    with screed.open(readfile) as fq:
        for n, is_pair, read1, read2 in ku.broken_paired_reader(fq):
            if is_pair:
                if notag(read1.name) in to_extract:
                    ku.write_record_pair(read1, read2, stdout)
            else:
                read = read1 if read1 else read2
                if notag(read.name) in to_extract:
                    ku.write_record(read, stdout)
    print("Done", file=stderr)


def cliparse():
    import argparse as ap
    psr = ap.ArgumentParser(
        "extract reads mapped to region of reference, including pairs",
    )
    psr.add_argument('-r', '--region', action='append',
                     help='Region to extract')
    psr.add_argument('fastq',
                     help='Interleaved fastq of reads')
    psr.add_argument('bam',
                     help='Alignment file')
    return psr.parse_args()


if __name__ == "__main__":
    args = cliparse()
    process_file(args.bam, args.region, args.fastq)
