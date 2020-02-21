SHELL := /bin/bash

CC = mpic++
LD ?= $(CC)

ROOT = ../../../CubismNova

dump ?= 0
debug ?= 0

CPPFLAGS = -DCUBISM_ALIGNMENT=16
ifeq "$(debug)" "1"
	CPPFLAGS += -std=c++11 -I$(ROOT)/include -g
else
	CPPFLAGS += -std=c++11 -I$(ROOT)/include -g -O3 -DNDEBUG
endif

OBJ += nova.cpp

ifeq "$(dump)" "1"
	CPPFLAGS += -D_DUMP_ -DCUBISM_USE_HDF
	LIBS += -L$(ROOT)/lib -lCubismIO
endif

.DEFAULT_GOAL := nova

all: nova

nova: nova.cpp
	$(CC) $(CPPFLAGS) $(extra) $^ -o $@ $(LIBS)

nova-profile: nova.cpp
	$(CC) $(CPPFLAGS) -pg $(extra) $^ -o $@ $(LIBS)
	./nova-profile
	gprof ./nova-profile > profile

%.o: %.cpp
	$(CC) $(CPPFLAGS) -c $^ -o $@

clean:
	rm -f *.o nova nova-profile
	rm -f *.h5 *.xmf

cleandeep: clean
	find . -iname "*~" -exec rm -f {} \;
