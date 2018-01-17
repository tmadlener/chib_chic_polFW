CC = g++

BUILDDIR := build
SRCDIR := src
INCDIR := interface
TARGET := bin
GENSRCDIR := ../general/src
GENINCDIR := ../general/interface

ROOTFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)

#CFLAGSINC := $(ROOTFLAGS) -I $(INCDIR) -ggdb
CFLAGSINC := $(ROOTFLAGS) -I $(INCDIR) -I $(GENINCDIR)

LIB := $(ROOTLIBS)

defaulttarget:
	@echo "No default target implemented."
	@echo "Available targets:"
	@echo "   preselection [cleanpreselection]"

#####################
# CHIB PRESELECTION #
#####################
preselection : $(BUILDDIR)/preselection.o $(BUILDDIR)/paralooper.o $(BUILDDIR)/treeproc.o
	$(CC) $^ -o $(TARGET)/$@ $(CFLAGSINC) $(LIB)
$(BUILDDIR)/preselection.o: $(SRCDIR)/chib_preselection.cc $(INCDIR)/chib_preselection.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c
cleanpreselection:
	rm -f $(BUILDDIR)/preselection.o $(BUILDDIR)/paralooper.o $(BUILDDIR)/treeproc.o $(TARGET)/preselection



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