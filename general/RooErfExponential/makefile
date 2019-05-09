ROOTLIBS=$(shell $(ROOT_CONFIG_BIN) --libs)
ROOTCXXFLAGS=$(shell $(ROOT_CONFIG_BIN) --cflags)
ROOTINCS=-I$(shell $(ROOT_CONFIG_BIN) --incdir)

CXXFLAGS=-std=c++14 -Wall -Wextra -Wpedantic -Wshadow -O3
ROOTFITLIBS=-lRooFit -lRooFitCore

CXX=$(GCC_BASE_DIR)g++

all: dictRooErfExponential.cxx RooErfExponential.so

dictRooErfExponential.cxx: RooErfExponential.h LinkDef.h
	rootcling -f $@ -c $^

RooErfExponential.so: dictRooErfExponential.cxx RooErfExponential.cxx
	$(CXX) $(ROOTINCS) -fPIC -shared $(CXXFLAGS) $(ROOTCXXFLAGS) -o $@ $(ROOTLIBS) $(ROOTFITLIBS) $^

clean:
	rm dictRooErfExponential* RooErfExponential.so
