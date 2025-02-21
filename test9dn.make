PROJECT=int2         # for historical reasons we always like to
                     # call the executable int2, but you could set it to
                     # something more descriptive

###HOST=linux-gfortran
HOST=linux-gfortran-openmp
###HOST=macos-gfortran
###HOST=macos-gfortran-openmp

include make.inc

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

FSRCS =  test9dn.f chunkmatc.f inter2dn.f hank103.f \
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
#
$(PROJECT):   $(MODS)   $(OBJS)
	rm -f $(PROJECT)
	$(FLINK) $(OBJS)
	./$(PROJECT)
#
clean: 
	rm -f $(OBJS)
# 
list: $(FSRCS)
	echo $^
#
pack: $(FSRCS)
	cat $^ > _tmp_.pack.f
#
