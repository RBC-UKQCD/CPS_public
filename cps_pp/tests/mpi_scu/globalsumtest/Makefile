#include<config.h>
CPS_START_NAMESPACE
##--------------------------------------------------------------------
##  CVS keywords
##
##  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/mpi_scu/globalsumtest/Makefile,v $
##  $State: Exp $
##
##--------------------------------------------------------------------
MAKEFILE = Makefile

PWD:sh = pwd

#  LIST contains the directories to search for other makefiles
#  The directory lib is last in the LIST

LIST1:sh = ls | awk '{if ($1 != "lib") print $0}
LIST = $(LIST1) lib

.SILENT :	

all	:=  ARG = all
force	:=  ARG = force
clean	:=  ARG = clean

all	: recurse
force	: recurse
clean	: recurse

recurse:
	for i in $(LIST) ; \
	  do if [ -f $$i/Makefile ] ; \
	      then cd $$i ; \
	      echo Running make $(ARG) in `pwd` ; \
	      make $(ARG) ; \
	      cd $(PWD) ; \
	    fi ; \
	done

CPS_END_NAMESPACE
