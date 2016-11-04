#!/usr/bin/env python3
from collections import namedtuple
import pysam
from sys import argv, stderr, stdout


Record = namedtuple('Record', ['name', 'seq', 'pair', 'rname', 'pos'])


def printrec(rec, file=None):
    print(">{}/{}\n{}".format(rec.name, rec.pair, rec.seq), file=file)


def printpair(a, b, file=None):
    if a.pair == '1':
        printrec(a, file)
        printrec(b, file)
    else:
        printrec(b, file)
        printrec(a, file)


def read2rec(read):
    name = read.query_name
    readno = '1' if read.is_read1 else '2'
    return Record(name, read.seq, readno, read.next_reference_name,
                  read.next_reference_start)


def process_file(bamfiles, regions, outstream=stdout):
    tofetch = {}
    for bamfile in bamfiles:
        samf = pysam.AlignmentFile(bamfile, 'rb')
        for region in regions:
            for read in samf.fetch(region):
                if read.is_unmapped or read.is_secondary:
                    continue
                this = read2rec(read)
                if not read.is_paired:
                    printrec(this)
                if this.name in tofetch:
                    other = tofetch.pop(this.name)
                    printpair(this, other, file=outstream)
                else:
                    tofetch[this.name] = this

        for name, have in tofetch.items():
            for read in samf.fetch(have.rname, have.pos - 10, have.pos + 10):
                need = read2rec(read)
                if need.name == have.name and need.pair != have.pair:
                    printpair(have, need, file=outstream)


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
