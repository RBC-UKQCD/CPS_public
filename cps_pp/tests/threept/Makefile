#include<config.h>
CPS_START_NAMESPACE
##--------------------------------------------------------------------
##  CVS keywords
##
##  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/threept/Makefile,v $
##  $State: Exp $
##
##--------------------------------------------------------------------
# The interface header file is:
INTERFACE_INCLUDE_FILE = ../../util/include/*.h ../../alg/include/*.h ../../task/include/*.h ../../nga/include/*.h

# The interface library is:
LIB = ../../lib/alg.lib \
	../../lib/util.lib \
	../../lib/d_op_stag_cpp.lib \
	../../lib/d_op_stag_types.lib \
	../../lib/d_op_clover.lib \
	../../lib/lapack.lib \
	../../lib/d_op_dwf.lib \
	../../lib/d_op_wilson_scl.lib \
	../../lib/d_op_wilson_types.lib \
	../../lib/d_op_base.lib \
	../../lib/glb_sim.lib \
	../../lib/scu_sim.lib \
	../../lib/mem.lib \
	../../lib/cbuf_sim.lib

TOF_LIB = ../../lib/alg.olb \
	../../lib/util.olb \
	../../lib/d_op_stag_opt.olb \
	../../lib/d_op_stag_types.olb \
	../../lib/d_op_clover.olb \
	../../lib/lapack.olb \
	../../lib/d_op_dwf_nos.olb \
	../../lib/d_op_wilson_opt_nos.olb \
	../../lib/d_op_base.olb \
	../../lib/d_op_wilson_types.olb \
	../../lib/glb_cpp.olb \
	../../lib/scu.olb \
	../../lib/mem.olb \
	../../lib/cbuf.olb


# The runtime library is:
RTS_LIB = tcio30bs.olb

# The linker control file is:
LCF = ../../mem/include/link_p2v.lcf

# If there are local *.h files uncoment the line below
#INCLUDES_LOCAL:sh	= ls *.h



#----------------------------------------------------------
# No editing below this line is needed
#----------------------------------------------------------
MAKEFILE	= Makefile

CPP	= CC
TCPP	= tcpp

CPPFLAGS = 
TCPPFLAGS = -mb -km

INCLUDES 	= $(INTERFACE_INCLUDE_FILE) $(INCLUDES_LOCAL)

SRCS:sh		= ls *.C
OBJS		= $(SRCS:.C=.o)
TOF_OBJS	= $(SRCS:.C=.tof)

OUT		= sun
TOF_OUT		= qcdsp

.SUFFIXES:

.SUFFIXES: .C .tof $(SUFFIXES)

all : $(OUT).out  $(TOF_OUT).out

force : clean all

$(OUT).out : $(OBJS) $(LIB) $(INCLUDES) $(MAKEFILE)
	$(CPP) $(CPPFLAGS) $(OBJS) $(LIB) -lm -o $@

$(TOF_OUT).out : $(TOF_OBJS)  $(INCLUDES) $(MAKEFILE) $(TOF_LIB) $(LCF)
	- tcpp  -e $(TOF_OUT).outtof $(TCPPFLAGS) $(TOF_OBJS) $(TOF_LIB) $(RTS_LIB) -fl $(LCF)
	t2c $(TOF_OUT).outtof $(TOF_OUT).out

$(OBJS) : $(INCLUDES) $(MAKEFILE)
$(TOF_OBJS) : $(INCLUDES) $(MAKEFILE)

.C.o	:
	$(CPP) -c $(CPPFLAGS) $<

.C.tof	:
	$(TCPP) -c $(TCPPFLAGS) $<

clean:
	rm -f *.o *.tof *.out *.outtof *.map *.ctl core 

CPS_END_NAMESPACE
