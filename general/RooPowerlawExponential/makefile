ROOTLIBS=$(shell $(ROOT_CONFIG_BIN) --libs)
ROOTCXXFLAGS=$(shell $(ROOT_CONFIG_BIN) --cflags)
ROOTINCS=-I$(shell $(ROOT_CONFIG_BIN) --incdir)

CXXFLAGS=-std=c++14 -Wall -Wextra -Wpedantic -Wshadow -O3
ROOTFITLIBS=-lRooFit -lRooFitCore

CXX=$(GCC_BASE_DIR)g++

all: dictRooPowerlawExponential.cxx RooPowerlawExponential.so

dictRooPowerlawExponential.cxx: RooPowerlawExponential.h LinkDef.h
	rootcling -f $@ -c $^

RooPowerlawExponential.so: dictRooPowerlawExponential.cxx RooPowerlawExponential.cxx
	$(CXX) $(ROOTINCS) -fPIC -shared $(CXXFLAGS) $(ROOTCXXFLAGS) -o $@ $(ROOTLIBS) $(ROOTFITLIBS) $^

clean:
	rm dictRooPowerlawExponential* RooPowerlawExponential.so
