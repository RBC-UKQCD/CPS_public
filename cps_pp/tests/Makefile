#include<config.h>
CPS_START_NAMESPACE
##--------------------------------------------------------------------
##  CVS keywords
##
##  $Author: mcneile $
##  $Date: 2003-06-22 13:34:47 $
##  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/tests/Makefile,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
##  $Id: Makefile,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
##  $Name: not supported by cvs2svn $
##  $Locker:  $
##  $Log: not supported by cvs2svn $
##  Revision 1.2  2001/06/19 18:12:20  anj
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
