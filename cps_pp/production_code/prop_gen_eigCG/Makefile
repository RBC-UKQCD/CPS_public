# $Id: Makefile,v 1.1 2011-11-12 03:59:56 chulwoo Exp $
# Makefile.users.  Generated from Makefile.users.in by configure.
#   This is the makefile for all the commands
#   that are common for all testing applications.
#----------------------------------------------------------------------


SRCDIR = /gpfs/home1/chulwoo/CPS/v5_0_15/cps++
BUILDDIR = /gpfs/home1/chulwoo/CPS/v5_0_15/fen
QOS = 
QOSLIB = ${QOS}/quser/gcc-lib-user///
CC = /bgl/BlueLight/ppcfloor/bglsys/bin/mpicc
CXX = /bgl/BlueLight/ppcfloor/bglsys/bin/mpicxx
AS  = /bgl/BlueLight/ppcfloor/blrts-gnu/bin//powerpc-bgl-blrts-gnu-as
LDFLAGS = -I/bgl/BlueLight/ppcfloor/bglsys/include/ -L/bgl/BlueLight/ppcfloor/bglsys/lib/ -lrts.rts -L/gpfs/home1/chulwoo/CPS/gmp/fen//lib -lgmp -L/gpfs/scratch1/chulwoo/SciDAC//scidac/install/qmp-2.3.1-BGL/lib -lqmp -L/gpfs/scratch1/chulwoo/SciDAC/scidac/install/qio-2.3.8-BGL//lib -L/gpfs/scratch1/chulwoo/SciDAC//scidac/install/qmp-2.3.1-BGL/lib -lqio -llime -lqmp

me = $(notdir $(PWD))
BIN = BGL.x

#VPATH :=$(SRCDIR)/tests/$(me)



#
# include list for the Columbia code
#
INCLIST = -I${BUILDDIR} -I${SRCDIR}/include  -I/bgl/BlueLight/ppcfloor/bglsys/include/ -I/gpfs/home1/chulwoo/CPS/gmp/fen//include -I/gpfs/scratch1/chulwoo/SciDAC//scidac/install/qmp-2.3.1-BGL/include    -I/gpfs/scratch1/chulwoo/SciDAC/scidac/install/qio-2.3.8-BGL//include

CFLAGS= -g -O2 -O -Wall
CXXFLAGS=  -O -Wall
ASFLAGS= 
DFLAGS +=  -DGMP -DUSE_QMP -DUSE_QIO

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
