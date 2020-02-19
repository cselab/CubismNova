SHELL := /bin/bash

CC = mpic++
LD = $(CC)

ROOT = ../../LegacyCubism

dump ?= 0

bsx ?= 32
bsy ?= $(bsx)
bsz ?= $(bsx)
ap ?= float
omp ?= 0
align ?= 16

CPPFLAGS = -std=c++11 -I$(ROOT)/include -g -O3
CPPFLAGS += -DCUBISM_ALIGNMENT=$(align) -D_BLOCKSIZEX_=$(bsx) -D_BLOCKSIZEY_=$(bsy) -D_BLOCKSIZEZ_=$(bsz)
LIBS = -lm -lz

OBJ += legacy.cpp \
	   $(ROOT)/src/ArgumentParser.cpp \
	   $(ROOT)/src/Histogram.cpp

ifeq "$(ap)" "float"
	CPPFLAGS += -D_FLOAT_PRECISION_
endif

ifeq "$(omp)" "1"
	CPPFLAGS += -fopenmp
endif

ifeq "$(dump)" "1"
	CPPFLAGS += -D_DUMP_ -I$(HDF5_ROOT)/include -DCUBISM_USE_HDF
	LIBS += -L$(HDF5_ROOT)/lib -lhdf5
endif

LIBS += -ldl

.DEFAULT_GOAL := legacy

all: legacy

legacy: $(OBJ)
	$(CC) $(CPPFLAGS) $(extra) $^ -o $@ $(LIBS)

legacy-profile: $(OBJ)
	$(CC) $(CPPFLAGS) -pg $(extra) $^ -o $@ $(LIBS)
	./legacy-profile
	gprof ./legacy-profile > profile

%.o: %.cpp
	$(CC) $(CPPFLAGS) -c $^ -o $@

clean:
	rm -f *.o legacy legacy-profile Cubism/src/*.o
	rm -f *.h5 *.xmf

cleandeep: clean
	find . -iname "*~" -exec rm -f {} \;
