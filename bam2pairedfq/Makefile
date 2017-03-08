CXXFLAGS += -std=gnu++14 -O3 -isystem ext/seqan-library-2.1.1/include/ \
			-fopenmp -DSEQAN_HAS_ZLIB=1 -DSEQAN_HAS_OPENMP=1
LDFLAGS += -static -static-libgcc -static-libstdc++
LIBS += -lz -lbz2
PREFIX ?= /usr/local

PROG=bam2pairfq

.PHONY: all
all: $(PROG)

%: %.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $< $(LIBS)

.PHONY: clean
clean:
	rm -rf $(PROG)

.PHONY: install
install: $(PROG)
	install -d $(PREFIX)/bin
	install -m 755 $< $(PREFIX)/bin
