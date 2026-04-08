PROJECT=int2         # for historical reasons we always like to
                     # call the executable int2, but you could set it to
                     # something more descriptive

# Auto-pick HOST and CC_WRAP based on the OS so this build works
# unchanged on macOS (Homebrew gcc) and Linux (system gcc).
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  HOST ?= macos-gfortran-openmp-quad
  # Use the newest Homebrew gcc available; fall back to gcc-14.
  UNAME_M := $(shell uname -m)
  ifeq ($(UNAME_M),arm64)
    HOMEBREW_PREFIX := /opt/homebrew
  else
    HOMEBREW_PREFIX := /usr/local
  endif
  HOMEBREW_GCCS := $(wildcard $(HOMEBREW_PREFIX)/bin/gcc-1[0-9]*)
  CC_WRAP ?= $(if $(HOMEBREW_GCCS),$(lastword $(HOMEBREW_GCCS)),gcc-14)
else
  HOST ?= linux-gfortran-openmp-quad
  CC_WRAP ?= gcc
endif

include make.inc

# C compiler for omp_wrappers (quad timer fix)
CWFLAGS  = -O3 -fopenmp -DF2C_FLOAT128

.PHONY: $(PROJECT) clean list

.f.$(OBJSUF):
	$(FC) $(FFLAGS) $<

.f.$(MODSUF):
	$(FC) $(FFLAGS) $<

.SUFFIXES: $(MODSUF) .$(OBJSUF) .f .c

# SOURCE FILE LIST
#
vpath %.f .:../lib:../geo

# omp_wrappers.c is shipped with f2c and installed to $(HOME)/lib by
# the f2c top-level Makefile (`make install`).  Tell make to look there
# so we don't have to keep a stale copy in this directory.
vpath omp_wrappers.c $(HOME)/lib

FMODS =

FSRCS =  test7.f chunkmatc.f inter2dn.f hank103.f \
           cadavect.f hqsuppquad.f \
           cgmres6-rel.f \
           cqrsolve.f \
	   legeexps.f prini.f

ifeq ($(WITH_SECOND),1)
FSRCS += second-r8.f
endif

#
# object files list
MODS    =  $(FMODS:.f=.$(MODSUF))
OBJS    =  $(FMODS:.f=.$(OBJSUF)) $(FSRCS:.f=.$(OBJSUF))
OMP_WRAP = omp_wrappers.o
#
$(PROJECT):   $(MODS)   $(OBJS) $(OMP_WRAP)
	rm -f $(PROJECT)
	$(FLINK) $(OBJS) $(OMP_WRAP)
	./$(PROJECT)
#
omp_wrappers.o: omp_wrappers.c
	$(CC_WRAP) $(CWFLAGS) -c $<

clean:
	rm -f $(OBJS) $(OMP_WRAP)
#
list: $(FSRCS)
	echo $^
#
pack: $(FSRCS)
	cat $^ > _tmp_.pack.f
#
