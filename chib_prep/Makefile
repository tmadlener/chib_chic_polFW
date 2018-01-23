CC = g++

BUILDDIR := build
SRCDIR := src
INCDIR := interface
TARGET := bin
GENSRCDIR := ../general/src
GENINCDIR := ../general/interface

ROOTFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)

ROOFITLIBS :=  -lRooFit -lRooFitCore

#CFLAGSINC := $(ROOTFLAGS) -I $(INCDIR) -ggdb
CFLAGSINC := $(ROOTFLAGS) -I $(INCDIR) -I $(GENINCDIR)

LIB := $(ROOTLIBS)

defaulttarget:
	@echo "No default target implemented."
	@echo "Available targets:"
	@echo "   preselection [cleanpreselection]"
	@echo "   fitting [cleanfitting]"


#####################
# CHIB PRESELECTION 

preselection : $(BUILDDIR)/preselection.o $(BUILDDIR)/paralooper.o $(BUILDDIR)/treeproc.o
	$(CC) $^ -o $(TARGET)/$@ $(CFLAGSINC) $(LIB)
$(BUILDDIR)/preselection.o: $(SRCDIR)/chib_preselection.cc $(INCDIR)/chib_preselection.h  $(GENINCDIR)/calcAngles.h  $(GENINCDIR)/misc_utils.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c
cleanpreselection:
	rm -f $(BUILDDIR)/preselection.o $(BUILDDIR)/paralooper.o $(BUILDDIR)/treeproc.o $(TARGET)/preselection
#
#####################


#####################
# CHIB FITTING 

fitting : $(BUILDDIR)/fitting.o $(BUILDDIR)/fitter.o $(BUILDDIR)/fitanalyser.o $(BUILDDIR)/utils.o
	$(CC) $^ -o $(TARGET)/$@ $(CFLAGSINC) $(LIB) $(ROOFITLIBS)
$(BUILDDIR)/fitting.o: $(SRCDIR)/chib_fitting.cc
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c
cleanfitting:
	rm -f $(BUILDDIR)/fitting.o $(BUILDDIR)/fitter.o $(BUILDDIR)/fitanalyser.o $(TARGET)/fitting
#
#####################





#
# PARALLEL TREE LOOPER	
#
$(BUILDDIR)/paralooper.o: $(GENSRCDIR)/ParallelTreeLooper.cc $(GENINCDIR)/ParallelTreeLooper.h $(GENINCDIR)/utils.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c

#
# TREE PROCESSOR
#
$(BUILDDIR)/treeproc.o: $(GENSRCDIR)/TreeProcessor.cc $(GENINCDIR)/TreeProcessor.h $(GENINCDIR)/ParallelTreeLooper.h $(GENINCDIR)/utils.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c 
		
#
# FITTER
#
$(BUILDDIR)/fitter.o: $(GENSRCDIR)/Fitter.cc $(GENINCDIR)/Fitter.h $(GENINCDIR)/utils.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c


#
# FITANALYSER
#
$(BUILDDIR)/fitanalyser.o: $(GENSRCDIR)/FitAnalyser.cc $(GENINCDIR)/FitAnalyser.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c

#
# UTILS
#
$(BUILDDIR)/utils.o: $(GENSRCDIR)/utils.cc $(GENINCDIR)/utils.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c
