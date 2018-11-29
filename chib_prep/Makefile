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

dictbld=$(BUILDDIR)/dict
CUSTOMFUNCTIONS := $(dictbld).o $(BUILDDIR)/ExpGaussExp.o $(BUILDDIR)/DoubleSidedCB.o

defaulttarget:
	@echo "No default target implemented."
	@echo "Available targets:"
	@echo "   preselection [cleanpreselection]"
	@echo "   fitting [cleanfitting]"
	@echo "   output [cleanoutput]"
	@echo "   helper [cleanhelper]"

all: _dir_setup preselection fitting

clean: _dir_cleanup

#####################
# CHIB PRESELECTION

preselection : $(BUILDDIR)/preselection.o $(BUILDDIR)/paralooper.o $(BUILDDIR)/treeproc.o $(BUILDDIR)/utils.o
	$(CC) $^ -o $(TARGET)/$@ $(CFLAGSINC) $(LIB)
$(BUILDDIR)/preselection.o: $(SRCDIR)/chib_preselection.cc $(INCDIR)/chib_preselection.h  $(GENINCDIR)/calcAngles.h  $(GENINCDIR)/misc_utils.h $(looper_header)
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c
cleanpreselection: cleanutils  cleantreeproc cleanparalooper
	rm -f $(BUILDDIR)/preselection.o $(TARGET)/preselection
#
#####################


#####################
# CHIB FITTING

fitting : $(BUILDDIR)/fitting.o $(BUILDDIR)/fitter.o $(BUILDDIR)/fitanalyser.o $(BUILDDIR)/utils.o $(BUILDDIR)/chiorganizer.o $(CUSTOMFUNCTIONS)
	$(CC) $^ -o $(TARGET)/$@ $(CFLAGSINC) $(LIB) $(ROOFITLIBS)
$(BUILDDIR)/fitting.o: $(SRCDIR)/chib_fitting.cc $(fit_header)
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c
cleanfitting: cleanchiorganizer cleanutils cleanfitanalyser cleanfitter cleancustomfuncs
	rm -f $(BUILDDIR)/fitting.o $(TARGET)/fitting
#
#####################

#####################
# CHIB OUTPUT

output : $(BUILDDIR)/output.o $(BUILDDIR)/treemerger.o $(BUILDDIR)/paralooper.o $(BUILDDIR)/utils.o $(BUILDDIR)/chiorganizer.o
	$(CC) $^ -o $(TARGET)/$@ $(CFLAGSINC) $(LIB)
$(BUILDDIR)/output.o: $(SRCDIR)/chib_output.cc $(output_header)
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c
cleanoutput: cleanchiorganizer cleanutils cleantreemerger  cleanparalooper
	rm -f $(BUILDDIR)/output.o $(TARGET)/output
#
#####################


#####################
# HELPER TOOLS

helper : getworkspacefile fitresults fitresultplot

getworkspacefile: $(BUILDDIR)/getworkspacefile.o  $(BUILDDIR)/chiorganizer.o $(BUILDDIR)/utils.o
	$(CC) $^ -o $(TARGET)/$@ $(CFLAGSINC) $(LIB)
$(BUILDDIR)/getworkspacefile.o: $(SRCDIR)/getworkspacefile.cc $(INCDIR)/ChiOrganizer.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c
cleangetworkspacefile:
	rm -f $(BUILDDIR)/getworkspacefile.o $(TARGET)/getworkspacefile
	
fitresults: $(BUILDDIR)/fitresults.o  $(BUILDDIR)/chiorganizer.o $(BUILDDIR)/fitanalyser.o $(BUILDDIR)/fitter.o $(BUILDDIR)/utils.o $(CUSTOMFUNCTIONS)
	$(CC) $^ -o $(TARGET)/$@ $(CFLAGSINC) $(LIB) $(ROOFITLIBS)
$(BUILDDIR)/fitresults.o: $(SRCDIR)/chib_fitresults.cc $(INCDIR)/ChiOrganizer.h $(GENINCDIR)/FitAnalyser.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c
cleanfitresults:
	rm -f $(BUILDDIR)/fitresults.o $(TARGET)/fitresults

fitresultplot: $(BUILDDIR)/fitresult_plot.o $(BUILDDIR)/chiorganizer.o $(BUILDDIR)/fitanalyser.o $(BUILDDIR)/fitter.o $(BUILDDIR)/utils.o $(CUSTOMFUNCTIONS)
	$(CC) $^ -o $(TARGET)/$@ $(CFLAGSINC) $(LIB) $(ROOFITLIBS)	
$(BUILDDIR)/fitresult_plot.o: $(SRCDIR)/fitresult_plot.cc $(INCDIR)/ChiOrganizer.h $(GENINCDIR)/FitAnalyser.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c
cleanfitresultplot:
	rm -f $(BUILDDIR)/fitresult_plot.o $(TARGET)/fitresultplot


cleanhelper: cleanchiorganizer cleanfitanalyser cleanutils cleanfitresults cleangetworkspacefile cleancustomfuncs cleanfitresultplot
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
cleanparalooper:
	rm -f $(BUILDDIR)/paralooper.o

#
# TREE PROCESSOR
#
$(BUILDDIR)/treeproc.o: $(GENSRCDIR)/TreeProcessor.cc $(GENINCDIR)/TreeProcessor.h $(GENINCDIR)/ParallelTreeLooper.h $(GENINCDIR)/utils.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c
cleantreeproc:
	rm -f $(BUILDDIR)/treeproc.o

#
# TREE Merger
#
$(BUILDDIR)/treemerger.o: $(GENSRCDIR)/TreeMerger.cc $(GENINCDIR)/TreeMerger.h $(GENINCDIR)/ParallelTreeLooper.h $(GENINCDIR)/utils.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c
cleantreemerger:
	rm -f $(BUILDDIR)/treemerger.o

#
# FITTER
#
$(BUILDDIR)/fitter.o: $(GENSRCDIR)/Fitter.cc $(GENINCDIR)/Fitter.h $(GENINCDIR)/utils.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c
cleanfitter:
	rm -f $(BUILDDIR)/fitter.o


#
# FITANALYSER
#
$(BUILDDIR)/fitanalyser.o: $(GENSRCDIR)/FitAnalyser.cc $(GENINCDIR)/FitAnalyser.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c
cleanfitanalyser:
	rm -f $(BUILDDIR)/fitanalyser.o

#
# UTILS
#
$(BUILDDIR)/utils.o: $(GENSRCDIR)/utils.cc $(GENINCDIR)/utils.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c
cleanutils:
	rm -f $(BUILDDIR)/utils.o

#
# CHIORGANIZER
#
$(BUILDDIR)/chiorganizer.o: $(SRCDIR)/ChiOrganizer.cc $(INCDIR)/ChiOrganizer.h $(GENINCDIR)/utils.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c
cleanchiorganizer:
	rm -f $(BUILDDIR)/chiorganizer.o

	
	
#
# CUSTOM ROOFIT FUNCTIONS
#

# Add new custom functins to CUSTOMFUNCTIONS variable

dictpcm=dict_rdict.pcm
$(BUILDDIR)/ExpGaussExp.o : $(GENSRCDIR)/ExpGaussExp.cxx $(GENINCDIR)/ExpGaussExp.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c
$(BUILDDIR)/DoubleSidedCB.o : $(GENSRCDIR)/DoubleSidedCB.cxx $(GENINCDIR)/DoubleSidedCB.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c
$(dictbld).o: $(GENINCDIR)/ExpGaussExp.h $(GENINCDIR)/DoubleSidedCB.h  $(GENINCDIR)/linkdef.h
	rm -f $(dictbld).cxx
	rootcling $(dictbld).cxx -c $^
	$(CC) $(dictbld).cxx -o $@ $(CFLAGSINC) $(LIB) -c -I.
	cp $(BUILDDIR)/$(dictpcm) $(TARGET)
cleancustomfuncs:
	rm -f $(BUILDDIR)/$(dictpcm) $(TARGET)/$(dictpcm) $(CUSTOMFUNCTIONS) $(dictbld).o


#
# DIRECTORY SETUP AND CLEANUP
#
_dir_setup:
	mkdir -p $(TARGET) $(BUILDDIR)
_dir_cleanup:
	rm -rf $(TARGET) $(BUILDDIR)



cleanall: cleanchiorganizer cleanutils cleanfitanalyser cleanfitter cleantreemerger cleantreeproc cleanparalooper
