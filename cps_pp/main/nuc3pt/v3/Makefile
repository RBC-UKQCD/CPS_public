# $Id: Makefile,v 1.2 2011-02-26 00:33:24 chulwoo Exp $
# Makefile.users.  Generated from Makefile.users.in by configure.
#   This is the makefile for all the commands
#   that are common for all testing applications.
#----------------------------------------------------------------------


SRCDIR = /home/chulwoo/CPS/v5_0_8-RICC/cps++
BUILDDIR = /home/chulwoo/CPS/v5_0_8-RICC/Ds
QOS = 
QOSLIB = ${QOS}/quser/gcc-lib-user///
CC = mpicc
CXX = mpiCC
AS  = as
LDFLAGS =  -Lyes/lib -lgmp -L/home/chulwoo/scidac/qmp-2.3.1-Izubuchi-Ds/lib -lqmp -L/home/chulwoo/scidac/qio-2.3.8-Ds/lib -L/home/chulwoo/scidac/qmp-2.3.1-Izubuchi-Ds/lib -lqio -llime -lqmp -fopenmp -L/home/chulwoo/CPS/v5_0_8-RICC/cps++/aux/GotoBLAS2 -lgoto2 -lgfortran

me = $(notdir $(PWD))
BIN = NOARCH.x

#VPATH :=$(SRCDIR)/tests/$(me)



#
# include list for the Columbia code
#
INCLIST = -I${BUILDDIR} -I${SRCDIR}/include  -Iyes/include -I/home/chulwoo/scidac/qmp-2.3.1-Izubuchi-Ds/include    -I/home/chulwoo/scidac/qio-2.3.8-Ds/include

CFLAGS= -O3  -march=amdfam10 -mfpmath=sse -ffast-math -funroll-loops -fprefetch-loop-arrays -fomit-frame-pointer -O3 -msse -msse2 -msse3
CXXFLAGS= -O3  -march=amdfam10 -mfpmath=sse -ffast-math -funroll-loops -fprefetch-loop-arrays -fomit-frame-pointer -O3 -msse -msse2 -msse3
ASFLAGS= 
DFLAGS +=  -DUSE_SSE -DGMP -DUSE_QMP -DUSE_QIO

#
# Libraries for the Columbia code
#
# (These are for the scalar version of the code)
#
#

.PHONY: cps clean


LIBLIST =\
  $(BUILDDIR)/cps.a \
  $(BUILDDIR)/d_op_wilson_sse.a \
  $(BUILDDIR)/d_op_dwf.a \
  $(BUILDDIR)/d_op_dwf_sse.a \
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
