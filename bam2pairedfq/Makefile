CXXFLAGS += -std=gnu++14 -O3 -isystem ext/seqan-library-2.1.1/include/ \
			-fopenmp -DSEQAN_HAS_ZLIB=1 -DSEQAN_HAS_OPENMP=1
LDFLAGS += -static -static-libgcc -static-libstdc++
LIBS += -lz -lbz2

bin/bam2pairfq: bam2pairfq.cc
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $< $(LIBS)

.PHONY: clean
clean:
	rm -rf bin
