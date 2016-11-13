CXXFLAGS += -std=gnu++14 -O3 $(shell pkg-config --cflags seqan-2) \
			-fopenmp -DSEQAN_HAS_ZLIB=1 -DSEQAN_HAS_OPENMP=1
LDFLAGS += -static -static-libgcc -static-libstdc++ 					\
		   $(shell pkg-config --static --libs-only-L seqan-2)
LIBS += $(shell pkg-config --static --libs-only-l --libs-only-other seqan-2)

bin/bam2pairfq: bam2pairfq.cc
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $< $(LIBS)

.PHONY: clean
clean:
	rm -rf bin
