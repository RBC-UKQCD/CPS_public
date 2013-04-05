# $Id: Makefile,v 1.2 2013-04-05 17:46:30 chulwoo Exp $
# Makefile.users.  Generated from Makefile.users.in by configure.
#   This is the makefile for all the commands
#   that are common for all testing applications.
#----------------------------------------------------------------------


SRCDIR = /Users/tblum/src/cps++/CPS508-Lanczos/obj/..
BUILDDIR = /Users/tblum/src/cps++/CPS508-Lanczos/obj
QOS = 
QOSLIB = ${QOS}/quser/gcc-lib-user///
CC = mpicc -m64
CXX = mpicxx -m64
AS  = as
LDFLAGS =  -L/Users/tblum/src/cps++/v5_0_8-EM/build/qmp-2.3.1//lib -lqmp -L/usr/local/qio-2.3.8//lib -L/Users/tblum/src/cps++/v5_0_8-EM/build/qmp-2.3.1//lib -lqio -llime -lqmp
QOS = 
QOSLIB = ${QOS}/quser/gcc-lib-user///

AS  = as
LDFLAGS =   -L/Volumes/ssd2/izubuchi/src/scidac/qmp-2.3.1/inst/lib  -lqmp \
            -L/Volumes/ssd2/izubuchi/src/scidac/qio-2.3.8/inst/lib -lqio \
            -llime 

me = $(notdir $(PWD))
BIN = cps-lanczos.x


#
# include list for the Columbia code
#
INCLIST = -I${BUILDDIR} -I${SRCDIR}/include  -I/Users/tblum/src/cps++/v5_0_8-EM/build/qmp-2.3.1//include    -I/usr/local/qio-2.3.8//include

CFLAGS= -g -O2 -O2
CXXFLAGS=  -O2

#ICCOPT=  -O3 -xS  -fno-alias  -align -inline-max-total-size=10000 -ip  -i-static -fasm-blocks
ICCOPT=  -O3 -xS  -fno-alias  -align -inline-max-total-size=10000 -ip  -i-static -fasm-blocks
CC = mpicc -cc=icc $(ICCOPT)
CXX = mpicxx -CC=icpc $(ICCOPT)
CC_PATH = $(CXX)

CFLAGS=  -g -ggdb -O2
CXXFLAGS=  -g -ggdb -O2
ASFLAGS= 
DFLAGS +=  -DUSE_QMP -DUSE_QIO -DUSE_SSE



#
# Libraries for the Columbia code
#
# (These are for the scalar version of the code)
#
#

.PHONY: cps clean


LIBLIST =\
  $(BUILDDIR)/cps.a \

#
#  targets
#


all:$(BIN)

.SUFFIXES:
.SUFFIXES:  .o .C .S .c

CSRC :=$(wildcard *.c)
CCSRC :=$(wildcard *.C)
SSRC :=$(wildcard *.S)

COBJ=$(CSRC:.c=.o)
CCOBJ=$(CCSRC:.C=.o)
SOBJ=$(SSRC:.S=.o)

OBJS_SRC = $(SOBJ) $(CCOBJ) $(COBJ)
OBJS := $(notdir $(OBJS_SRC))

$(BIN):  $(OBJS) $(LIBLIST)
	@echo OBJS = $(OBJS)
	$(CXX) $(OBJS) $(LIBLIST) $(LDFLAGS) -o $(BIN)

.c.o:
	$(CC) -o $@ $(CFLAGS) $(DFLAGS) -c $(INCLIST) $<
.C.o:
	$(CXX) -o $@ $(CXXFLAGS) $(DFLAGS) -c $(INCLIST) $<
.S.o:
	$(AS) -o $@ $(ASFLAGS) -c $(INCLIST) $<

cps:
	$(MAKE) -C $(BUILDDIR)

clean:
	rm -f *.dat *.o  $(BIN)
	rm -f ../regressions/*$(me).dat
	rm -f ../regressions/*$(me).checklog
