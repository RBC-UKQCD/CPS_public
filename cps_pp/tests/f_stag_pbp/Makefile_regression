#
#  Simple makefile for regression testing for the generation
#  of unquenched gauge configuratios.
#
# Options:
#     make     : create the executable 
#
#


#
#  include the platform specific part of the makefile
#
include ../Makefile_common

#
#  targets
#

all:regression

regression:
	$(COMPILER) -c $(CFLAGS) $(INCLIST) -o main.o main.C
	$(COMPILER) -o regression_cps.x main.o $(LIBLIST) $(LDFLAGS)  



clean:
	rm *.dat *.o  regression_cps.x





