#include<config.h>
CPS_START_NAMESPACE
##--------------------------------------------------------------------
##  CVS keywords
##
##  $Author: chulwoo $
##  $Date: 2004-06-04 21:14:16 $
##  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/Makefile,v 1.2 2004-06-04 21:14:16 chulwoo Exp $
##  $Id: Makefile,v 1.2 2004-06-04 21:14:16 chulwoo Exp $
##  $Name: not supported by cvs2svn $
##  $Locker:  $
##  $Revision: 1.2 $
##  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/Makefile,v $
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
