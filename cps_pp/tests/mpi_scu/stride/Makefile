#include<config.h>
CPS_START_NAMESPACE
# Makefile for the MPI implementation of the UKQCD@QCDOC code.
#
# CVS v. $Id: Makefile,v 1.1.1.1 2003-06-22 13:34:49 mcneile Exp $
# A.N. Jackson: andrewj@epcc.ed.ac.uk 
#

# name of final library archive
LIBRARY = ../../lib/libcommsMPI.a

# source files 
SOURCES = stridetest.C \

# rules for compilation
#
# first you need to specify what suffices might appear
# add more separated by spaces (eg ' .cc' for c++ sources)
.SUFFIXES: .c .f90 .cxx .cc .C 

# Compiler details and flags and so on:
#
# The C compiler 
CC = /opt/SUNWhpc/bin/tmcc
# Flags for the C compiler
CFLAGS = -g -DVERBOSE
# Include paths for headers and libraries
CINCS = -I. -I/opt/SUNWhpc/include -I../../
CLINCS = -L. -L../../lib -L/opt/SUNWhpc/lib
CLIBS = -lmpi

# The C++ compiler
CPPC = /opt/SUNWhpc/bin/tmCC
CPPFLAGS = $(CFLAGS)
CPPINCS = $(CINCS)

# The F90 compiler
#FC = f90
#FFLAGS = -c

# now you need to specify how to get .o files from .whatever files
#
# for c sources. You could add other flags here
.c.o:
	$(CC) -c $(CFLAGS) $(CINCS) $(CLINCS) $?

# for fortran 90 source.
.f90.o:
	$(FC) -c $(FFLAGS) $?

# for C++ sources, left commented out for now:
.cc.o:
	$(CPPC) -c $(CPPINCS) $(CPPFLAGS) $?

.C.o:
	$(CPPC) -c $(CPPINCS) $(CPPFLAGS) $?


# This line links the objectfiles to make the executables
# You shouldn't need to touch this unless you want a different linker or
# specify some spangly options
#
stridetest: $(LIBRARY) $(SOURCES)
	$(CPPC) -o stridetest $(CPPFLAGS) $(SOURCES) $(CPPINCS) -L $(LIBRARY) $(CLIBS)

clean: 
	rm -f *.o *~ comms.log.* core

CPS_END_NAMESPACE
