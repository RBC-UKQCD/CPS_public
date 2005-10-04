# $Id: Makefile,v 1.2 2005-10-04 15:29:14 chulwoo Exp $
# tests/Makefile_common.  Generated from Makefile_common.in by configure.
#   This is the makefile for all the commands
#   that are common for all testing applications.
#----------------------------------------------------------------------


SRCDIR = /home/Q01/mike/v4_9_2_mike/cps++
BUILDDIR = /home/Q01/mike/v4_9_2_mike/qos-2.6.0g
QOS = /qcdoc/sfw/qos/v2.6.0/aix5.2g/
QOSLIB = ${QOS}/quser/gcc-lib-user///
CC = env QGCC_EXEC_PREFIX=${QOSLIB} /qcdoc/sfw/packages/v1.5/aix5.2/local/bin/powerpc-gnu-elf-gcc
CXX = env QGCC_EXEC_PREFIX=${QOSLIB} /qcdoc/sfw/packages/v1.5/aix5.2/local/bin/powerpc-gnu-elf-g++
AS  = env QGCC_EXEC_PREFIX=${QOSLIB} /qcdoc/sfw/packages/v1.5/aix5.2/local/bin/powerpc-gnu-elf-as
LDFLAGS =  -L/qcdoc/sfw/packages/gmp/qos-2.5.9a/lib -lgmp

me = $(notdir $(PWD))
BIN = QCDOC.x




#
# include list for the Columbia code
#
INCLIST = -I${BUILDDIR} -I${SRCDIR}/include  -I/qcdoc/sfw/packages/gmp/qos-2.5.9a/include

CFLAGS= -D_GCC_ -O3
CXXFLAGS= -D_GCC_ -O3
ASFLAGS= 
DFLAGS +=  -DGMP

#
# Libraries for the Columbia code
#
# (These are for the scalar version of the code)
#
#

.PHONY: cps clean


LIBLIST = $(BUILDDIR)/cps.a 
#LIBLIST =  $(BUILDDIR)/cps.xla $(BUILDDIR)/cps.a 
	 

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
	rm -f $(SRCDIR)/tests/$(me)/*.o
	rm -f ../regressions/*$(me).dat
	rm -f ../regressions/*$(me).checklog

