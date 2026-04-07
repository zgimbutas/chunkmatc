#
#  test7-f2c-omp.make -- build test7 via f2c + gcc with OpenMP
#
#  Usage:
#    make -f test7-f2c-omp.make
#    make -f test7-f2c-omp.make run
#    make -f test7-f2c-omp.make clean
#

PROJECT = int2

F2C     = $(HOME)/linux/repositories/f2c/src/f2c
F2CFLAGS = -fomp
F2C_INC = $(HOME)/linux/repositories/f2c
F2C_LIBDIR = $(HOME)/lib
LIBF2C  = $(HOME)/repositories/f2c/libf2c-x86_64/libf2c.a

# Portable defaults: pick a working gcc and MPFR/GMP path for the host OS.
# - macOS: auto-detect highest Homebrew gcc-N (Apple clang lacks __float128).
#          Use /opt/homebrew on arm64, /usr/local on Intel -- mixing them
#          picks an x86_64 gcc-N on Apple Silicon and breaks linking against
#          arm64 libf2c/mpfr/gmp.
# - Linux: use system gcc; rely on default library search for mpfr/gmp.
# Any of these can be overridden on the make command line.
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  UNAME_M := $(shell uname -m)
  ifeq ($(UNAME_M),arm64)
    HOMEBREW_PREFIX ?= /opt/homebrew
  else
    HOMEBREW_PREFIX ?= /usr/local
  endif
  HOMEBREW_GCCS := $(shell ls $(HOMEBREW_PREFIX)/bin/gcc-[0-9]* 2>/dev/null | sort -V)
  CC_F2C   ?= $(if $(HOMEBREW_GCCS),$(lastword $(HOMEBREW_GCCS)),gcc-14)
  MPFR_LIB ?= -L$(HOMEBREW_PREFIX)/lib -lmpfr
  GMP_LIB  ?= -L$(HOMEBREW_PREFIX)/lib -lgmp
else
  CC_F2C   ?= gcc
  MPFR_LIB ?= -lmpfr
  GMP_LIB  ?= -lgmp
endif

# gcc 15+ defaults to C23 where () means (void), which breaks f2c's
# K&R-style extern declarations.  Force gnu17 in that case.
GCC_MAJOR := $(shell $(CC_F2C) -dumpversion 2>/dev/null | cut -d. -f1)
ifeq ($(shell [ "$(GCC_MAJOR)" -ge 15 ] 2>/dev/null && echo y),y)
  GCC_STD := -std=gnu17
endif

CFLAGS  = -O3 -march=native -fopenmp $(GCC_STD) -I$(F2C_INC)
LDLIBS   = $(LIBF2C) $(MPFR_LIB) $(GMP_LIB) -lquadmath -lm -lgomp

SRCDIR   = .
LIBDIR   = $(HOME)/develop/lib
BUILDDIR = _f2c_build_omp

FSRCS = test7.f chunkmatc.f inter2dn.f hank103.f cadavect.f \
        hqsuppquad.f cgmres6-rel.f cqrsolve.f legeexps.f prini.f \
        $(LIBDIR)/second-r8.f

EXTRA_C = omp_wrappers.c

CSRCS = $(addprefix $(BUILDDIR)/, $(notdir $(FSRCS:.f=.c)))
OBJS  = $(CSRCS:.c=.o) $(addprefix $(BUILDDIR)/, $(EXTRA_C:.c=.o))

.PHONY: all clean run
.PRECIOUS: $(BUILDDIR)/%.c

all: $(BUILDDIR)/$(PROJECT)
	@echo "=== Built $(BUILDDIR)/$(PROJECT) ==="

run: $(BUILDDIR)/$(PROJECT)
	cd $(BUILDDIR) && OMP_NUM_THREADS=4 ./$(PROJECT)

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

$(BUILDDIR)/%.c: $(SRCDIR)/%.f | $(BUILDDIR)
	@cp $< $(BUILDDIR)/$(notdir $<)
	@cd $(BUILDDIR) && $(F2C) $(F2CFLAGS) $(notdir $<) 2>/dev/null; rm -f $(notdir $<)

$(BUILDDIR)/%.c: $(LIBDIR)/%.f | $(BUILDDIR)
	@cp $< $(BUILDDIR)/$(notdir $<)
	@cd $(BUILDDIR) && $(F2C) $(F2CFLAGS) $(notdir $<) 2>/dev/null; rm -f $(notdir $<)

$(BUILDDIR)/%.o: $(BUILDDIR)/%.c
	$(CC_F2C) $(CFLAGS) -c $< -o $@

$(BUILDDIR)/%.o: $(F2C_LIBDIR)/%.c | $(BUILDDIR)
	$(CC_F2C) $(CFLAGS) -c $< -o $@

$(BUILDDIR)/$(PROJECT): $(OBJS)
	$(CC_F2C) $(CFLAGS) -o $@ $(OBJS) $(LDLIBS)
	@echo "  LINK $@"

clean:
	rm -rf $(BUILDDIR)
