CXXFLAGS += -std=gnu++14 -O3 -DSEQAN_HAS_ZLIB=1 -DSEQAN_HAS_BZIP2=1 -fsanitize=address
#LDFLAGS += -static -static-libstdc++ -static-libgcc
LIBS += -lz -lbz2 -pthread

bam2pairfq: bam2pairfq.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $< $(LIBS)

.PHONY: clean
clean:
	rm -f bam2pairfq
