SHELL := /bin/bash

CC = mpic++
LD = $(CC)

ROOT = ${HOME}/code/amrex

dump ?= 0
pg ?= 0
ap ?= float

omp ?= 0
align ?= 16

CPPFLAGS = -std=c++11 -g -O3 -DNDEBUG
LIBS =

OBJ += amrex.cpp

target = amrex

ifeq "$(ap)" "float"
ifeq "$(pg)" "0"
	target = amrex
	# CPPFLAGS += -I$(ROOT)/gnu3D_SP/include
	# LIBS += $(ROOT)/gnu3D_SP/lib/libamrex.a
	CPPFLAGS += -I$(ROOT)/gnu3D_SP_NOMPI/include
	LIBS += $(ROOT)/gnu3D_SP_NOMPI/lib/libamrex.a
else
	target = amrex-profile
	CPPFLAGS += -I$(ROOT)/gnu3D_gprof_SP/include
	LIBS += $(ROOT)/gnu3D_gprof_SP/lib/libamrex.a
endif
else
ifeq "$(pg)" "0"
	target = amrex
	CPPFLAGS += -I$(ROOT)/gnu3D_DP/include
	LIBS += $(ROOT)/gnu3D_DP/lib/libamrex.a
else
	target = amrex-profile
	CPPFLAGS += -I$(ROOT)/gnu3D_gprof_DP/include
	LIBS += $(ROOT)/gnu3D_gprof_DP/lib/libamrex.a
endif
endif

ifeq "$(dump)" "1"
	CPPFLAGS += -D_DUMP_
endif

LIBS += -lgfortran -pthread -lm -lz -ldl

.DEFAULT_GOAL := amrex

all: amrex

amrex: $(OBJ)
	$(CC) $(CPPFLAGS) $(extra) $^ -o $@ $(LIBS)

amrex-profile: $(OBJ)
	$(CC) $(CPPFLAGS) -pg $(extra) $^ -o $@ $(LIBS)
	./amrex-profile
	gprof ./amrex-profile > profile

%.o: %.cpp
	$(CC) $(CPPFLAGS) -c $^ -o $@

clean:
	rm -f *.o amrex amrex-profile
	rm -rf plt* Backtrace.*

cleandeep: clean
	find . -iname "*~" -exec rm -f {} \;
