ROOTLIBS=$(shell $(ROOT_CONFIG_BIN) --libs)
ROOTCXXFLAGS=$(shell $(ROOT_CONFIG_BIN) --cflags)

CXXFLAGS=-std=c++14 -Wall -Wextra -Wpedantic -Wshadow -O3 -g # -Wconversion # can't use conversion warnings, because there are a lot in ROOT
ROOFITLIBS=-lRooFit -lRooFitCore -lMinuit -lFoam

CXX=$(GCC_BASE_DIR)g++

INCS=-I./interface -I../ -I.

BIN=./bin/
SRC=./src/

SRCS=$(wildcard $(SRC)*.cc)

PROG=$(SRCS:$(SRC)%.cc=%)
LIST=$(addprefix $(BIN), $(PROG))

DEPS=$(LIST:%=%.d)

-include $(DEPS)

all: setup $(LIST)

setup:
	mkdir -p $(BIN)

$(BIN)%: $(SRC)%.cc
	$(CXX) $(INCS) $(ROOTINCS) $(CXXFLAGS) $(ROOTCXXFLAGS) -MMD -MP $< -o $@ $(ROOTLIBS) $(ROOFITLIBS)
