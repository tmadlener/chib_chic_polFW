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

CUSTOMFUNCTIONS := $(BUILDDIR)/ExpGaussExp.o $(BUILDDIR)/DoubleSidedCB.o

defaulttarget:
	@echo "No default target implemented"
	@echo "Available targets:"
	@echo "   paralooper [cleanparalooper]"
	@echo "   treeproc [cleantreeproc]"
	@echo "   treemerger [cleantreemerger]"
	@echo "   fitter [cleanfitter]"
	@echo "   sharedfunctions [cleansharedfunctions]"

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
fitter : $(BUILDDIR)/fitterexample.o $(BUILDDIR)/fitter.o $(BUILDDIR)/utils.o $(BUILDDIR)/fitanalyser.o $(CUSTOMFUNCTIONS) $(BUILDDIR)/dict.o
	$(CC) $^ -o $(TARGET)/$@ $(CFLAGSINC) $(LIB)  -lRooFit -lRooFitCore -lRooStats
$(BUILDDIR)/fitterexample.o: $(TESTDIR)/fitter_example.cc $(INCDIR)/Fitter.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c
cleanfitter: cleancustomfuncs
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
# CUSTOM ROOFIT FUNCTIONS
#

# Add new custom functins to CUSTOMFUNCTIONS variable

dictbld=$(BUILDDIR)/dict
dictpcm=dict_rdict.pcm
$(BUILDDIR)/ExpGaussExp.o : $(SRCDIR)/ExpGaussExp.cxx $(INCDIR)/ExpGaussExp.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c
$(BUILDDIR)/DoubleSidedCB.o : $(SRCDIR)/DoubleSidedCB.cxx $(INCDIR)/DoubleSidedCB.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c
$(dictbld).o: $(INCDIR)/ExpGaussExp.h $(INCDIR)/DoubleSidedCB.h  $(INCDIR)/linkdef.h
	rm -f $(dictbld).cxx
	rootcling $(dictbld).cxx -c $^
	$(CC) $(dictbld).cxx -o $@ $(CFLAGSINC) $(LIB) -c -I.
	cp $(BUILDDIR)/$(dictpcm) $(TARGET)
cleancustomfuncs:
	rm -f $(BUILDDIR)/$(dictpcm) $(TARGET)/$(dictpcm) $(CUSTOMFUNCTIONS) $(dictbld).o

sharedfunctions: $(TARGET)/CustomRooFitFunctions.so
dictbldshared=$(BUILDDIR)/dictshared
dictpcmshared=dictshared_rdict.pcm
$(TARGET)/CustomRooFitFunctions.so : $(BUILDDIR)/ExpGaussExp_shared.o  $(BUILDDIR)/DoubleSidedCB_shared.o $(dictbldshared).o
	$(CC) -shared $^ $(CFLAGSINC) $(LIB) -o $@ -lRooFit -lRooFitCore -lRooStats
$(BUILDDIR)/ExpGaussExp_shared.o : $(SRCDIR)/ExpGaussExp.cxx $(INCDIR)/ExpGaussExp.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c -fPIC
$(BUILDDIR)/DoubleSidedCB_shared.o : $(SRCDIR)/DoubleSidedCB.cxx $(INCDIR)/DoubleSidedCB.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c -fPIC
$(dictbldshared).o: $(INCDIR)/ExpGaussExp.h $(INCDIR)/DoubleSidedCB.h $(INCDIR)/linkdef.h
	rm -f $(dictbldshared).cxx
	rootcling $(dictbldshared).cxx -c -fPIC $^
	$(CC) $(dictbldshared).cxx -o $@ $(CFLAGSINC) $(LIB) -c -I. -fPIC
	cp $(BUILDDIR)/$(dictpcmshared) $(TARGET)
cleansharedfunctions:
	rm -f $(BUILDDIR)/$(dictpcmshared) $(TARGET)/$(dictpcmshared) $(TARGET)/CustomRooFitFunctions.so \
	$(BUILDDIR)/DoubleSidedCB_shared.o $(BUILDDIR)/ExpGaussExp_shared.o $(dictbldshared).o
