from collections import namedtuple
import pysam
from sys import argv, stderr, stdout


Record = namedtuple('Record', ['seq', 'pair', 'rname', 'pos'])


def printpair(name, a, b, file=None):
    first = a
    second = b
    if a.pair == '2':
        first = b
        second = a
    print(">{}/{}\n{}".format(name, first.pair, first.seq))
    print(">{}/{}\n{}".format(name, second.pair, second.seq))


def process_file(filename, region, outstream=stdout):
    samf = pysam.AlignmentFile(filename, 'rb')
    tofetch = {}
    for read in samf.fetch(region):
        name = read.query_name
        if not read.is_paired:
            print(">{}/{}\n{}\n".format(name, 1, read.seq), file=outstream)
        readno = '1' if read.is_read1 else '2'
        this = Record(read.seq, readno, read.next_reference_name,
                      read.next_reference_start)
        if name in tofetch:
            other = tofetch.pop(name)
            printpair(name, this, other, file=outstream)
        else:
            tofetch[name] = this

    for name, rec in tofetch.items():
        for read in samf.fetch(rec.rname, rec.pos - 10, rec.pos + 10):
            if read.query_name == name:
                readno = '1' if read.is_read1 else '2'
                other = Record(read.seq, readno, read.next_reference_name,
                               read.next_reference_start)
                printpair(name, rec, other, file=outstream)

if __name__ == "__main__":
    process_file(argv[1], argv[2])
