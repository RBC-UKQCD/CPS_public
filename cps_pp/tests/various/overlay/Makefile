#include<config.h>
CPS_START_NAMESPACE
##--------------------------------------------------------------------
##  CVS keywords
##
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
