#include<config.h>
CPS_START_NAMESPACE
##--------------------------------------------------------------------
##  CVS keywords
##
##  $Author: mcneile $
##  $Date: 2003-06-22 13:34:47 $
##  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/various/overlay/Makefile,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
##  $Id: Makefile,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
##  $Name: not supported by cvs2svn $
##  $Locker:  $
##  $Log: not supported by cvs2svn $
##  Revision 1.2  2001/06/19 18:12:31  anj
##  Serious ANSIfication.  Plus, degenerate double64.h files removed.
##  Next version will contain the new nga/include/double64.h.  Also,
##  Makefile.gnutests has been modified to work properly, propagating the
##  choice of C++ compiler and flags all the way down the directory tree.
##  The mpi_scu code has been added under phys/nga, and partially
##  plumbed in.
##
##  Everything has newer dates, due to the way in which this first alteration was handled.
##
##  Anj.
##
##  Revision 1.2  2001/05/25 06:16:04  cvs
##  Added CVS keywords to phys_v4_0_0_preCVS
##
##  $RCSfile: Makefile,v $
##  $Revision: 1.1.1.1 $
##  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/various/overlay/Makefile,v $
##  $State: Exp $
##
##--------------------------------------------------------------------
# The runtime library is:
RTS_LIB = tcio30bs.olb


#----------------------------------------------------------
# No editing below this line is needed
#----------------------------------------------------------
MAKEFILE	= Makefile

TCPP	= tcpp
TI    	= asm30
TA	= tasm

TAFLAGS		= -v30 -mb
TIFLAGS		= -q -mb -v30 -lxs -dSCU_ON=1 -dCBUF=0 -dROUND=1 -dSTAND_ALONE=0

CPPFLAGS = 
TCPPFLAGS = -mb -km

INCLUDES 	= $(INTERFACE_INCLUDE_FILE) $(INCLUDES_LOCAL)

SRCS_CPP:sh	=  echo ; ls | gawk '/\.C$/{print}';
SRCS_C:sh	=  echo ; ls | gawk '/\.c$/{print}';
SRCS_ASM:sh	=  echo ; ls | gawk '/\.asm$/{print}';
SRCS_TASM:sh	=  echo ; ls | gawk '/\.tasm$/{print}';

TOF_OBJS_CPP	= $(SRCS_CPP:.C=.tof)
TOF_OBJS_C	= $(SRCS_C:.c=.tof)
TOF_OBJS_ASM	= $(SRCS_ASM:.asm=.tof)
TOF_OBJS_TASM	= $(SRCS_TASM:.tasm=.tof)
TOF_OBJS	= $(TOF_OBJS_CPP) $(TOF_OBJS_C) $(TOF_OBJS_ASM) $(TOF_OBJS_TASM)

TOF_OUT		= qcdsp

.SUFFIXES:

.SUFFIXES: .C .tof .asm .tasm $(SUFFIXES)

all : $(TOF_OUT).out

force : clean all

$(TOF_OUT).out : $(TOF_OBJS)  $(INCLUDES) $(MAKEFILE) link.lcf
	tcpp  -e $(TOF_OUT).outtof $(TCPPFLAGS) $(TOF_OBJS) $(RTS_LIB)
	t2c $(TOF_OUT).outtof $(TOF_OUT).out

$(TOF_OBJS) : $(INCLUDES) $(MAKEFILE)


.C.tof	:
	$(TCPP) -c $(TCPPFLAGS) $<

.asm.tof:
	$(TI) $(TIFLAGS) $<
	c2t $*.obj $*.tof

.tasm.tof:
	$(TA) $(TAFLAGS) $<

clean:
	rm -f *.o *.tof *.out *.outtof *.map *.ctl  *.lst *.obj core 

CPS_END_NAMESPACE
