#include<config.h>
CPS_START_NAMESPACE
##--------------------------------------------------------------------
##  CVS keywords
##
##  $Author: mcneile $
##  $Date: 2003-06-22 13:34:49 $
##  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/doc/Makefile,v 1.1.1.1 2003-06-22 13:34:49 mcneile Exp $
##  $Id: Makefile,v 1.1.1.1 2003-06-22 13:34:49 mcneile Exp $
##  $Name: not supported by cvs2svn $
##  $Locker:  $
##  $RCSfile: Makefile,v $
##  $Revision: 1.1.1.1 $
##  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/doc/Makefile,v $
##  $State: Exp $
##
##--------------------------------------------------------------------
MAKEFILE = Makefile

PWD:sh = pwd

#  LIST contains the directories to search for other makefiles
#  The directory lib is last in the LIST

LIST1:sh = ls 
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
