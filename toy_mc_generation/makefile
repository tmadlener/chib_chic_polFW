ROOTLIBS=$(shell $(ROOT_CONFIG_BIN) --libs)
ROOTCXXFLAGS=$(shell $(ROOT_CONFIG_BIN) --cflags)

ROOFITLIBS=-lRooFit -lRooFitCore

CXXFLAGS=-std=c++14 #-Wall -Wextra -Wpedantic -Wshadow -O3

CXX=$(GCC_BASE_DIR)g++
INCS=-I./interface/

all: run_chicpolgen

run_chicpolgen: run_chicpolgen.cc chicpolgen.C ./interface/smearing.h ./interface/efficiencies.h ./interface/range.h ./interface/select.h
	$(CXX) $(INCS) $(ROOTLIBS) $(ROOTCXXFLAGS) $(CXXFLAGS) -o $@ $< $(ROOFITLIBS)

clean:
	rm run_chicpolgen


newgen: chibpolgeneration.o ChibPolGenerator.o
	$(CXX) $^ -o $@ $(INCS) $(ROOTLIBS) $(ROOTCXXFLAGS) $(CXXFLAGS) $(ROOFITLIBS)
chibpolgeneration.o: chibpolgeneration.cc ChiPolGenerator.h
	$(CXX) $< -o $@ $(INCS) $(ROOTLIBS) $(ROOTCXXFLAGS) $(CXXFLAGS) $(ROOFITLIBS) -c
ChibPolGenerator.o: ChiPolGenerator.cc ChiPolGenerator.h
	$(CXX) $< -o $@ $(INCS) $(ROOTLIBS) $(ROOTCXXFLAGS) $(CXXFLAGS) $(ROOFITLIBS) -c
