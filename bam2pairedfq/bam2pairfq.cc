/*
 * Copyright © 2016 Kevin Murray <kdmfoss@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <vector>
#include <string>
#include <iostream>
#include <unordered_set>
#include <getopt.h>


#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>

using namespace seqan;
using namespace std;

int region2tid(BamFileIn &bam, const CharString &region)
{
    int tid = 0;
    if (!getIdByName(tid, contigNamesCache(context(bam)), region)) {
        return -1;
    }
    return tid;
}

bool extract_readids(std::string &bamfile, std::unordered_set<string> &regions,
                     std::unordered_set<string> &readids)
{
    cerr << "Determining which reads in '" << bamfile << "' align to selected regions\n";

    BamFileIn bam_in;
    if (!open(bam_in, bamfile.c_str())) {
        cerr << "ERROR: Could not open " << bamfile << endl;
        return false;
    }

    BamIndex<Bai> bam_idx;
    if (!open(bam_idx, (bamfile + string(".bai")).c_str())) {
        cerr << "ERROR: Could not open bam index " << bamfile << endl;
        return false;
    }

    BamHeader header;
    readHeader(header, bam_in);

    // Check all regions are valid
    for (auto &region : regions) {
        if (region2tid(bam_in, region) < 0) {
            cerr << "ERROR: Unknown region " << region << endl;
            return false;
        }
    }

    for (auto &region : regions) {
        BamAlignmentRecord record;

        int tid = region2tid(bam_in, region);

        bool have_alns = false;
        if (!jumpToRegion(bam_in, have_alns, tid, 0, 1, bam_idx)) {
            cerr << "ERROR: Failed to seek to region " << region << endl;
            return false;
        }
        size_t n_alns = 0;
        size_t n_readids = readids.size();
        if (!have_alns) goto doneregion;

        while (!atEnd(bam_in)) {
            readRecord(record, bam_in);

            // If we are on the next reference or at the end already then we stop.
            if (record.rID != tid) break;

            n_alns++;
            readids.insert(toCString(record.qName));
        }
doneregion:
        cerr << "  Finished region '" << region << "'\n";
        cerr << "  Processed " << n_alns << " alignments\n";
        cerr << "  Added " << readids.size() - n_readids << " new read IDs to extract\n";
    }
    return true;
}

static inline bool readid_matches(string readid,
                                  const std::unordered_set<string> &readids)
{
    size_t pos = readid.find(" ");
    if (pos != string::npos) readid.erase(pos);
    pos = readid.rfind("/");
    if (pos == readid.size() - 2) readid.erase(pos);
    return readids.find(readid) != readids.end();
}

bool print_matching_reads(const string &fastqfile, const std::unordered_set<string> &readids)
{
    SeqFileIn reads_in(fastqfile.c_str());
    SeqFileOut reads_out(cout, Fastq());

    cerr << "Extracting reads with matching IDs from '" << fastqfile << "'\n";

    string id, seq, qual;
    size_t n_reads = 0, n_match = 0;
    while (!atEnd(reads_in)) {
        n_reads++;
        readRecord(id, seq, qual, reads_in);
        if (readid_matches(id, readids)) {
            n_match++;
            writeRecord(reads_out, id, seq, qual);
        }
    }
    cerr << "  Finished extracting reads\n";
    cerr << "  Processed " << n_reads << " reads\n";
    cerr << "  Output " << n_match << " matching reads\n";
}

static const string usage = 
    "USAGE:\n  bam2pairfq -r region FASTQ BAM\n\n"
    "Extracts reads (pairs) that align to `region`. One can supply -r\n"
    "many times with different regions.";  // no \n, done below

int
main(int argc, char *argv[])
{
    int ret = EXIT_FAILURE;
    int c = 0;
    std::unordered_set<string> regions;
    while ((c = getopt(argc, argv, "r:h")) > 0) {
        switch (c) {
        case 'r':
            regions.insert(optarg);
            break;
        case 'h':
            ret = EXIT_SUCCESS;
        default:
            cerr << usage << endl;
            return ret;
            break;
        }
    }
    if (optind + 2 > argc) {
        cerr << usage << endl;
        return ret;
    }

    string fastqfile = argv[optind];
    string bamfile = argv[optind + 1];

    std::unordered_set<string> readids;

    if (!extract_readids(bamfile, regions, readids)) return ret;
    if (!print_matching_reads(fastqfile, readids)) return ret;

    ret = EXIT_SUCCESS;
    return ret;
}
