ROOTLIBS=$(shell $(ROOT_CONFIG_BIN) --libs)
ROOTCXXFLAGS=$(shell $(ROOT_CONFIG_BIN) --cflags)

ROOFITLIBS=-lRooFit -lRooFitCore

CXXFLAGS=-std=c++14 -Wall -Wextra -Wpedantic -Wshadow -O3

CXX=$(GCC_BASE_DIR)g++

all: run_chicpolgen

run_chicpolgen: run_chicpolgen.cc chicpolgen.C smearing.h
	$(CXX) $(ROOTLIBS) $(ROOTCXXFLAGS) $(CXXFLAGS) -o $@ $< $(ROOFITLIBS)

clean:
	rm run_chicpolgen
