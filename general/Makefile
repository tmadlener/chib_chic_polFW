CC = g++

BUILDDIR := build
SRCDIR := src
INCDIR := interface
TARGET := bin
TESTDIR := testing
LIBTARGET := lib

ROOTFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)

CFLAGSINC := $(ROOTFLAGS) -I $(INCDIR) -ggdb
#CFLAGSINC := $(ROOTFLAGS) -I $(INCDIR)

LIB := $(ROOTLIBS)

defaulttarget:
	@echo "No default target implemented"
	@echo "Available targets:"
	@echo "   paralooper [cleanparalooper]"
	@echo "   treeproc [cleantreeproc]"
	@echo "   treemerger [cleantreemerger]"
	@echo "   fitter [cleanfitter]"

################################
# PARALLEL TREE LOOPER EXAMPLE #
################################
paralooper : $(BUILDDIR)/testparalooper.o $(BUILDDIR)/paralooper.o
	$(CC) $^ -o $(TARGET)/$@ $(CFLAGSINC) $(LIB)
$(BUILDDIR)/testparalooper.o: $(TESTDIR)/paralleltreelooper_example.cc $(INCDIR)/ParallelTreeLooper.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c
cleanparalooper:
	rm -f $(BUILDDIR)/paralooper.o $(BUILDDIR)/testparalooper.o $(TARGET)/paralooper
	
##########################
# TREE PROCESSOR EXAMPLE #
##########################
treeproc : $(BUILDDIR)/testtreeproc.o $(BUILDDIR)/paralooper.o $(BUILDDIR)/treeproc.o
	$(CC) $^ -o $(TARGET)/$@ $(CFLAGSINC) $(LIB)
$(BUILDDIR)/testtreeproc.o: $(TESTDIR)/treeprocessor_example.cc $(INCDIR)/TreeProcessor.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c
cleantreeproc:
	rm -f $(BUILDDIR)/testtreeproc.o $(BUILDDIR)/paralooper.o $(BUILDDIR)/treeproc.o $(TARGET)/treeproc


##########################
# TREE MERGER EXAMPLE #
##########################
treemerger : $(BUILDDIR)/testtreemerger.o $(BUILDDIR)/paralooper.o $(BUILDDIR)/treemerger.o
	$(CC) $^ -o $(TARGET)/$@ $(CFLAGSINC) $(LIB)
$(BUILDDIR)/testtreemerger.o: $(TESTDIR)/treemerger_example.cc $(INCDIR)/TreeMerger.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c
cleantreemerger:
	rm -f $(BUILDDIR)/testtreemerger.o $(BUILDDIR)/paralooper.o $(BUILDDIR)/treemerger.o $(TARGET)/treemerger


##################
# FITTER EXAMPLE #
##################
fitter : $(BUILDDIR)/fitterexample.o $(BUILDDIR)/fitter.o $(BUILDDIR)/utils.o $(BUILDDIR)/fitanalyser.o $(BUILDDIR)/ExpGaussExp.o $(BUILDDIR)/dict.o
	$(CC) $^ -o $(TARGET)/$@ $(CFLAGSINC) $(LIB)  -lRooFit -lRooFitCore -lRooStats
$(BUILDDIR)/fitterexample.o: $(TESTDIR)/fitter_example.cc $(INCDIR)/Fitter.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c
cleanfitter: cleanexpgaussexp
	rm -f $(BUILDDIR)/fitterexample.o $(BUILDDIR)/fitter.o $(TARGET)/fitter



#
# PARALLEL TREE LOOPER	
#
$(BUILDDIR)/paralooper.o: $(SRCDIR)/ParallelTreeLooper.cc $(INCDIR)/ParallelTreeLooper.h $(INCDIR)/utils.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c
	
#
# TREE PROCESSOR
#
$(BUILDDIR)/treeproc.o: $(SRCDIR)/TreeProcessor.cc $(INCDIR)/TreeProcessor.h $(INCDIR)/ParallelTreeLooper.h $(INCDIR)/utils.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c

#
# TREE MERGER
#
$(BUILDDIR)/treemerger.o: $(SRCDIR)/TreeMerger.cc $(INCDIR)/TreeMerger.h $(INCDIR)/ParallelTreeLooper.h $(INCDIR)/utils.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c
	
#
# FITTER
#
$(BUILDDIR)/fitter.o: $(SRCDIR)/Fitter.cc $(INCDIR)/Fitter.h $(INCDIR)/utils.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c


#
# FITANALYSER
#
$(BUILDDIR)/fitanalyser.o: $(SRCDIR)/FitAnalyser.cc $(INCDIR)/FitAnalyser.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c

#
# UTILS
#
$(BUILDDIR)/utils.o: $(SRCDIR)/utils.cc $(INCDIR)/utils.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c

#
# EXPGAUSSEXP ROOFIT FUNCTION
#
dictbld=$(BUILDDIR)/dict
dictpcm=dict_rdict.pcm
$(BUILDDIR)/ExpGaussExp.o : $(SRCDIR)/ExpGaussExp.cxx $(INCDIR)/ExpGaussExp.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c
$(dictbld).o: $(INCDIR)/ExpGaussExp.h  $(INCDIR)/linkdef.h
	rm -f $(dictbld).cxx
	rootcling $(dictbld).cxx -c $^
	$(CC) $(dictbld).cxx -o $@ $(CFLAGSINC) $(LIB) -c -I.
	cp $(BUILDDIR)/$(dictpcm) $(TARGET)
cleanexpgaussexp:
	rm -f $(BUILDDIR)/$(dictpcm) $(TARGET)/$(dictpcm) $(BUILDDIR)/ExpGaussExp.o $(dictbld).o

	
