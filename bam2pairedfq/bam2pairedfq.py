#!/usr/bin/env python3
# bam2pairedfastq -- extract paired (or single) reads if one maps to some region
# Copyright (c) 2016 Kevin Murray <kdmfoss@gmail.com>
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
# OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# Uses khmer's utils module for a broken paired reader.

import pysam
import screed
import khmer_utils as ku
from sys import argv, stderr, stdout


def notag(readid):
    '''Read name to "basename" of read, as appears in sam/bam'''
    readid, _, _ = readid.partition(' ')
    if readid[-2] == '/' and readid[-1] in "12":
        return readid[:-2]
    return readid


def process_file(bamfile, regions, readfile, outstream=stdout):
    print("Extracting reads from ", readfile, "that align to", regions, "in",
          bamfile,  file=stderr)
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
