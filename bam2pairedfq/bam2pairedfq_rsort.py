#!/usr/bin/env python3
import pysam
from sys import argv, stderr, stdout


def tophred(quals, offset=33):
    return "".join([chr(x+offset) for x in quals])


def printread(read, file=None):
    pair = ""
    if read.is_paired:
        pair = "/{}".format(1 if read.is_read1 else 2)
    tag = ""
    if not read.is_unmapped:
        tag = ' {}:{}-{}'.format(read.reference_name, read.reference_start + 1,
                                 read.reference_end)
    print("@", read.query_name, pair, tag, sep='', file=file)
    print(read.query_sequence, file=file)
    print("+", file=file)
    print(tophred(read.query_qualities), file=file)


def process_file(bamfiles, regions, outstream=stdout):
    for bamfile in bamfiles:
        samf = pysam.AlignmentFile(bamfile, 'rb')
        tids = {}
        for reg in regions:
            _, tid, start, end = samf.parse_region(reg)
            tids[tid] = (start, end)
        looking_for = None
        for read in samf:
            name = read.query_name
            if name == looking_for and not read.is_secondary and not read.is_supplementary:
                printread(read)
                looking_for = None
                continue
            tid = read.reference_id
            next_tid = read.next_reference_id
            if tid in tids or next_tid in tids:
                printread(read)
                looking_for = name


def cliparse():
    import argparse as ap
    psr = ap.ArgumentParser(
        "extract reads mapped to region of reference, including pairs",
    )
    psr.add_argument('-r', '--region', action='append',
                     help='Region to extract')
    psr.add_argument('bam', nargs='+',
                     help='Alignment file')
    return psr.parse_args()


if __name__ == "__main__":
    args = cliparse()
    process_file(args.bam, args.region)
