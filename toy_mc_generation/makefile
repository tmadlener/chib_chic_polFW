ROOTLIBS=$(shell $(ROOT_CONFIG_BIN) --libs)
ROOTCXXFLAGS=$(shell $(ROOT_CONFIG_BIN) --cflags)

ROOFITLIBS=-lRooFit -lRooFitCore

CXXFLAGS=-std=c++14 #-Wall -Wextra -Wpedantic -Wshadow -O3

CXX=$(GCC_BASE_DIR)g++
GENINCDIR := ../general/interface
MCINCDIR := ./interface
INCS=-I$(MCINCDIR) -I$(GENINCDIR)

DEPENDENCIES := $(MCINCDIR)/smearing.h $(MCINCDIR)/efficiencies.h $(MCINCDIR)/range.h $(MCINCDIR)/select.h
DEPENDENCIES += $(GENINCDIR)/calcAngles_defs.h $(GENINCDIR)/calcAngles.h

defaulttarget:
	@echo "No default target implemented."
	@echo "Available targets:"
	@echo "   all"
	@echo "   run_chicpolgen [clean]"
	@echo "   chibgen [cleanchibgen]"

all: run_chicpolgen

run_chicpolgen: run_chicpolgen.cc chicpolgen.C $(DEPENDENCIES) $(GENINCDIR)/ArgParser.h
	$(CXX) $(INCS) $(ROOTLIBS) $(ROOTCXXFLAGS) $(CXXFLAGS) -o $@ $< $(ROOFITLIBS)

clean:
	rm run_chicpolgen

chibgen: generate_chibpol.o ChiPolarizationGenerator.o 
	$(CXX) $^ -o $@ $(INCS) $(ROOTLIBS) $(ROOTCXXFLAGS) $(CXXFLAGS) $(ROOFITLIBS)
generate_chibpol.o: generate_chibpol.cc ChiPolarizationGenerator.h  $(GENINCDIR)/ArgParser.h $(GENINCDIR)/utils.h
	$(CXX) $< -o $@ $(INCS) $(ROOTLIBS) $(ROOTCXXFLAGS) $(CXXFLAGS) $(ROOFITLIBS) -c
ChiPolarizationGenerator.o: ChiPolarizationGenerator.cc ChiPolarizationGenerator.h $(DEPENDENCIES) $(GENINCDIR)/utils.h
	$(CXX) $< -o $@ $(INCS) $(ROOTLIBS) $(ROOTCXXFLAGS) $(CXXFLAGS) $(ROOFITLIBS) -c

cleanchibgen:
	rm -f chibgen generate_chibpol.o ChiPolarizationGenerator.o
