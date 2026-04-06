PROJECT=int2         # for historical reasons we always like to
                     # call the executable int2, but you could set it to
                     # something more descriptive

###HOST=linux-gfortran-openmp-quad
HOST=macos-gfortran-openmp-quad

include make.inc

# C compiler for omp_wrappers (quad timer fix)
CC_WRAP  = gcc-14
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
