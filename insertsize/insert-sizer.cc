/*
* ============================================================================
 *
 *       Filename:  insert-sizer.cc
 *    Description:  Calcualte insert size distribution from a BAM
 *        License:  GPLv3+
 *         Author:  Kevin Murray, spam@kdmurray.id.au
 *
 * ============================================================================
 */

#include <iostream>
#include <map>

#include <bamtools/api/BamReader.h>
#include <bamtools/api/BamAlignment.h>

class InsertSizeCalculator
{
public:
    InsertSizeCalculator();

    bool open(const std::string &bamfile);

    bool run();

    bool report();

protected:

    void _add_alignment(BamTools::BamAlignment &aln);

    BamTools::BamReader     _bam;
    std::map<size_t, size_t>
                            _histogram;
    std::map<int32_t, std::map<size_t, std::map<size_t, size_t>>>
                            _locus_sizes;
    size_t                  _num_alignments;  // All alignments, good or bad.
    size_t                  _num_paired_alignments;  // Good, PE alignments.
    size_t                  _num_merged_alignments;  // Good, SE merged reads.
};

InsertSizeCalculator::
InsertSizeCalculator()
    : _num_alignments(0)
    , _num_paired_alignments(0)
    , _num_merged_alignments(0)
{}

bool
InsertSizeCalculator::
open(const std::string &bamfile)
{
    return _bam.Open(bamfile);
}


bool
InsertSizeCalculator::
run()
{
    BamTools::BamAlignment aln;

    while (_bam.GetNextAlignment(aln)) {
        _num_alignments++;
        _add_alignment(aln);
    }
    return true;
}

void
InsertSizeCalculator::
_add_alignment(BamTools::BamAlignment &aln)
{
    if (aln.IsFailedQC()) {
        return;
    }
    if (!aln.IsMapped()) {
        return;
    }
    if (!aln.IsPrimaryAlignment()) {
        return;
    }

    if (aln.InsertSize > 0) {
        // InsertSize > 0 is true IFF read is paired, both pairs map and this
        // is the first of the pairs.
        if (aln.RefID != aln.MateRefID) {
            return;
        }
#if 0
        std::cerr << aln.Name << " " << aln.RefID << " " << aln.MateRefID << " ";
        std::cerr << aln.Position << " ";
        std::cerr << "pair -- P: " << aln.IsPaired() << " IS: " << aln.InsertSize << std::endl;
#endif

        _num_paired_alignments++;
        _histogram[aln.InsertSize] += 1;
        _locus_sizes[aln.RefID][aln.Position][aln.InsertSize] += 1;
        return;
    }
#if 0
    } else if (aln.IsFirstMate() && !aln.IsMateMapped()) {
        size_t frag_length = aln.AlignedBases.size();

//#if 0
        _names[aln.Name] += 1;
        std::cerr << aln.Name << " " << _names[aln.Name] << " ";
        std::cerr << aln.Position << " ";
        std::cerr << "merged -- IMM: " << aln.IsMateMapped() << " fl: " << frag_length << std::endl;
//#endif

        _histogram[frag_length] += 1;
        _num_merged_alignments++;
        _locus_sizes[aln.Position][frag_length] += 1;
        return;
    }
#endif
}

bool
InsertSizeCalculator::
report()
{
    using std::cout;
    using std::endl;
#if 0
    cout << "Processed " << _num_alignments << " alignments." << endl;
    cout << _num_paired_alignments << " were paired." << endl;
    cout << _num_merged_alignments << " were merged." << endl;
#endif

    cout << "RefID\tLocus\tSize\tNumTags\tModeSizeNumTags" << endl;
    for (auto &refpair: _locus_sizes) {
        int32_t             ref = refpair.first;
        for (auto &pair: refpair.second) {
            size_t              locus = pair.first;
            std::map<size_t, size_t>
                                counts = pair.second;
            size_t              locus_count = 0;
            size_t              mode_size = 0;
            size_t              max_count = 0;
            for (auto &countpair: counts) {
                size_t size = countpair.first;
                size_t count = countpair.second;

                locus_count += count;
                if (max_count < count) {
                    mode_size = size;
                    max_count = count;
                }

            }
            if (max_count > 1) {
                cout << ref << "\t" << locus << "\t" << mode_size << "\t" <<
                        locus_count  << "\t" << max_count << endl;
            }
        }
    }
    return true;
}


int
main (int argc, char *argv[])
{
    InsertSizeCalculator isc;

    if (argc != 2) {
        std::cerr << "USAGE: " << argv[0] << " <bamfile>" << std::endl;
        if (argc == 1) {
            // We had no args, which is our version of --help. So it's not a
            // failure.
            return EXIT_SUCCESS;
        }
        return EXIT_FAILURE;
    }
    if (!isc.open(argv[1])) {
        std::cerr << "Error opening BAM '" << argv[1] << "'" << std::endl;
        return EXIT_FAILURE;
    }
    if (!isc.run()) {
        std::cerr << "Error calcuating insert sizes." << std::endl;
        return EXIT_FAILURE;
    }
    isc.report();
    return EXIT_SUCCESS;
}
