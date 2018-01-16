CC = g++

BUILDDIR := build
SRCDIR := src
INCDIR := interface
TARGET := bin
TESTDIR := testing
LIBTARGET := lib

ROOTFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)

#CFLAGSINC := $(ROOTFLAGS) -I $(INCDIR) -ggdb
CFLAGSINC := $(ROOTFLAGS) -I $(INCDIR)

LIB := $(ROOTLIBS)

#############################
# PARALLEL TREE LOOPER TEST #
#############################
paralooper : $(BUILDDIR)/testparalooper.o $(BUILDDIR)/paralooper.o
	$(CC) $^ -o $(TARGET)/$@ $(CFLAGSINC) $(LIB)
$(BUILDDIR)/testparalooper.o: $(TESTDIR)/paralleltreelooper_example.cc $(INCDIR)/ParallelTreeLooper.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c	
$(BUILDDIR)/paralooper.o: $(SRCDIR)/ParallelTreeLooper.cc $(INCDIR)/ParallelTreeLooper.h $(INCDIR)/utils.h
	$(CC) $< -o $@ $(CFLAGSINC) $(LIB) -c
cleanparalooper:
	rm -f $(BUILDDIR)/paralooper.o $(BUILDDIR)/testparalooper.o paralooper