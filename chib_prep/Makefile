CC = g++

BUILDDIR := build
SRCDIR := src
INCDIR := interface
TARGET := bin
GENSRCDIR := ../general/src
GENINCDIR := ../general/interface

ROOTFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)

ROOFITLIBS :=  -lRooFit -lRooFitCore -lRooStats

#CFLAGSINC := $(ROOTFLAGS) -I $(INCDIR) -ggdb
CFLAGSINC := $(ROOTFLAGS) -I $(INCDIR) -I $(GENINCDIR)

LIB := $(ROOTLIBS)

defaulttarget:
	@echo "No default target implemented."
	@echo "Available targets:"
	@echo "   preselection [cleanpreselection]"
	@echo "   fitting [cleanfitting]"
	@echo "   output [cleanoutput]"

all: _dir_setup preselection fitting

clean: _dir_cleanup

#####################
# CHIB PRESELECTION

preselection : $(BUILDDIR)/preselection.o $(BUILDDIR)/paralooper.o $(BUILDDIR)/treeproc.o
	$(CC) $^ -o $(TARGET)/$@ $(CFLAGSINC) $(LIB)
$(BUILDDIR)/preselection.o: $(SRCDIR)/chib_preselection.cc $(INCDIR)/chib_preselection.h  $(GENINCDIR)/calcAngles.h  $(GENINCDIR)/misc_utils.h $(looper_header)
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c
cleanpreselection:
	rm -f $(BUILDDIR)/preselection.o $(BUILDDIR)/paralooper.o $(BUILDDIR)/treeproc.o $(TARGET)/preselection
#
#####################


#####################
# CHIB FITTING

fitting : $(BUILDDIR)/fitting.o $(BUILDDIR)/fitter.o $(BUILDDIR)/fitanalyser.o $(BUILDDIR)/utils.o $(BUILDDIR)/chiorganizer.o
	$(CC) $^ -o $(TARGET)/$@ $(CFLAGSINC) $(LIB) $(ROOFITLIBS)
$(BUILDDIR)/fitting.o: $(SRCDIR)/chib_fitting.cc $(fit_header)
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c
cleanfitting:
	rm -f $(BUILDDIR)/fitting.o $(BUILDDIR)/fitter.o $(BUILDDIR)/fitanalyser.o $(TARGET)/fitting
#
#####################

#####################
# CHIB OUTPUT

output : $(BUILDDIR)/output.o $(BUILDDIR)/treemerger.o $(BUILDDIR)/paralooper.o $(BUILDDIR)/utils.o
	$(CC) $^ -o $(TARGET)/$@ $(CFLAGSINC) $(LIB)
$(BUILDDIR)/output.o: $(SRCDIR)/chib_output.cc $(output_header)
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c
cleanoutput:
	rm -f $(BUILDDIR)/output.o $(BUILDDIR)/treemerger.o $(BUILDDIR)/utils.o $(BUILDDIR)/paralooper.o $(TARGET)/output
#
#####################




_looper_header := ParallelTreeLooper.h TreeProcessor.h utils.h
looper_header := $(patsubst %,$(GENINCDIR)/%,$(_looper_header))

_output_header := ParallelTreeLooper.h TreeMerger.h utils.h
output_header := $(patsubst %,$(GENINCDIR)/%,$(_output_header))

_fit_header := FitAnalyser.h Fitter.h utils.h
fit_header := $(patsubst %,$(GENINCDIR)/%,$(_fit_header)) $(INCDIR)/ChiOrganizer.h

all_header := $(fit_header) $(looper_header) $(output_header)

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
# TREE Merger
#
$(BUILDDIR)/treemerger.o: $(GENSRCDIR)/TreeMerger.cc $(GENINCDIR)/TreeMerger.h $(GENINCDIR)/ParallelTreeLooper.h $(GENINCDIR)/utils.h
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

#
# CHIORGANIZER
#
$(BUILDDIR)/chiorganizer.o: $(SRCDIR)/ChiOrganizer.cc $(INCDIR)/ChiOrganizer.h $(GENINCDIR)/utils.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c

#
# DIRECTORY SETUP AND CLEANUP
#
_dir_setup:
	mkdir -p $(TARGET) $(BUILDDIR)
_dir_cleanup:
	rm -rf $(TARGET) $(BUILDDIR)
