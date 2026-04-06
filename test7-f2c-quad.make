#
#  test7-f2c-quad.make -- build test7 via f2c + gcc with quad precision
#
#  Usage:
#    make -f test7-f2c-quad.make
#    make -f test7-f2c-quad.make run
#    make -f test7-f2c-quad.make clean
#

PROJECT = int2

F2C     = $(HOME)/linux/repositories/f2c/src/f2c
F2CFLAGS = -freal-8-real-16
F2C_INC = $(HOME)/linux/repositories/f2c
LIBF2C  = $(HOME)/repositories/f2c/libf2c-x86_64/libf2c.a

CC_F2C  = gcc-14
CFLAGS  = -O3 -march=native -DF2C_FLOAT128 -I$(F2C_INC)

MPFR_LIB = -L/opt/homebrew/Cellar/mpfr/4.2.2/lib -lmpfr
GMP_LIB  = -L/opt/homebrew/Cellar/gmp/6.3.0/lib -lgmp
LDLIBS   = $(LIBF2C) $(MPFR_LIB) $(GMP_LIB) -lquadmath -lm

SRCDIR   = .
LIBDIR   = $(HOME)/develop/lib
BUILDDIR = _f2c_build_quad

FSRCS = test7.f chunkmatc.f inter2dn.f hank103.f cadavect.f \
        hqsuppquad.f cgmres6-rel.f cqrsolve.f legeexps.f prini.f \
        $(LIBDIR)/second-r8.f

CSRCS = $(addprefix $(BUILDDIR)/, $(notdir $(FSRCS:.f=.c)))
OBJS  = $(CSRCS:.c=.o)

.PHONY: all clean run
.PRECIOUS: $(BUILDDIR)/%.c

all: $(BUILDDIR)/$(PROJECT)
	@echo "=== Built $(BUILDDIR)/$(PROJECT) ==="

run: $(BUILDDIR)/$(PROJECT)
	cd $(BUILDDIR) && ./$(PROJECT)

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

$(BUILDDIR)/$(PROJECT): $(OBJS)
	$(CC_F2C) $(CFLAGS) -o $@ $(OBJS) $(LDLIBS)
	@echo "  LINK $@"

clean:
	rm -rf $(BUILDDIR)
