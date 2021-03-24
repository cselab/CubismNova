SHELL := /bin/bash

CC = mpic++
LD = $(CC)

ROOT = ${HOME}/code/aphros/install

dump ?= 0
pg ?= 0
ap ?= float

omp ?= 0
align ?= 16

CPPFLAGS = -std=c++14 -g -O3 -DNDEBUG
LIBS =

CPPFLAGS += -I$(ROOT)/include -I$(ROOT)/include/hydro
LIBS += $(ROOT)/lib/aphros_approx_eb.so \
		$(ROOT)/lib/aphros_approx.so \
		$(ROOT)/lib/aphros_color.so \
		$(ROOT)/lib/aphros_convdiffe.so \
		$(ROOT)/lib/aphros_convdiffi.so \
		$(ROOT)/lib/aphros_convdiffvg.so \
		$(ROOT)/lib/aphros_cubismnc.so \
		$(ROOT)/lib/aphros_cubism.so \
		$(ROOT)/lib/aphros_curv.so \
		$(ROOT)/lib/aphros_distrbasic.so \
		$(ROOT)/lib/aphros_distr.so \
		$(ROOT)/lib/aphros_distrsolver.so \
		$(ROOT)/lib/aphros_dumper.so \
		$(ROOT)/lib/aphros_embed.so \
		$(ROOT)/lib/aphros_events.so \
		$(ROOT)/lib/aphros_gitgen.so \
		$(ROOT)/lib/aphros_git.so \
		$(ROOT)/lib/aphros_histogram.so \
		$(ROOT)/lib/aphros_hypre.so \
		$(ROOT)/lib/aphros_hypresub.so \
		$(ROOT)/lib/aphros_init.so \
		$(ROOT)/lib/aphros_local.so \
		$(ROOT)/lib/aphros_normal.so \
		$(ROOT)/lib/aphros_overlap.so \
		$(ROOT)/lib/aphros_parser.so \
		$(ROOT)/lib/aphros_partstrmeshm.so \
		$(ROOT)/lib/aphros_posthook_default.so \
		$(ROOT)/lib/aphros_proj_eb.so \
		$(ROOT)/lib/aphros_proj.so \
		$(ROOT)/lib/aphros_report.so \
		$(ROOT)/lib/aphros_simple.so \
		$(ROOT)/lib/aphros_solver.so \
		$(ROOT)/lib/aphros_subcomm.so \
		$(ROOT)/lib/aphros_suspender.so \
		$(ROOT)/lib/aphros_sysinfo.so \
		$(ROOT)/lib/aphros_timer.so \
		$(ROOT)/lib/aphros_tvd.so \
		$(ROOT)/lib/aphros_utilconvdiff.so \
		$(ROOT)/lib/aphros_utilfluid.so \
		$(ROOT)/lib/aphros_utilhydro.so \
		$(ROOT)/lib/aphros_utilvof.so \
		$(ROOT)/lib/aphros_vars.so \
		$(ROOT)/lib/aphros_vofm.so \
		$(ROOT)/lib/aphros_vof.so \
		$(ROOT)/lib/aphros_young.so \
		$(ROOT)/lib/libHYPRE.so \
		$(ROOT)/lib/libh5.so \
		$(ROOT)/lib/libmarch.so \
		$(ROOT)/lib/liboverlap.so \
		$(ROOT)/lib/libvof.so \
		$(ROOT)/lib/libvofi.so \
		$(ROOT)/lib/libvtk.a \
		$(ROOT)/lib/libtable.a \


OBJ += aphros.cpp


ifeq "$(ap)" "float"
	CPPFLAGS += -D_SINGLE_PRECISION_
endif

ifeq "$(dump)" "1"
	CPPFLAGS += -D_DUMP_
endif

LIBS += -lm -lz -ldl

.DEFAULT_GOAL := aphros

all: aphros

aphros: $(OBJ)
	$(CC) $(CPPFLAGS) $(extra) $^ -o $@ $(LIBS)

aphros-profile: $(OBJ)
	$(CC) $(CPPFLAGS) -pg $(extra) $^ -o $@ $(LIBS)
	./aphros-profile
	gprof ./aphros-profile > profile

%.o: %.cpp
	$(CC) $(CPPFLAGS) -c $^ -o $@

clean:
	rm -f *.o aphros aphros-profile
	rm -rf *.xmf *.h5

cleandeep: clean
	find . -iname "*~" -exec rm -f {} \;
