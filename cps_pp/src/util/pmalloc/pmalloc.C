#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:46 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/pmalloc/pmalloc.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Id: pmalloc.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.2  2001/06/19 18:13:34  anj
//  Serious ANSIfication.  Plus, degenerate double64.h files removed.
//  Next version will contain the new nga/include/double64.h.  Also,
//  Makefile.gnutests has been modified to work properly, propagating the
//  choice of C++ compiler and flags all the way down the directory tree.
//  The mpi_scu code has been added under phys/nga, and partially
//  plumbed in.
//
//  Everything has newer dates, due to the way in which this first alteration was handled.
//
//  Anj.
//
//  Revision 1.2  2001/05/25 06:16:10  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: pmalloc.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/pmalloc/pmalloc.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
CPS_END_NAMESPACE
#include<util/pmalloc.h>
#include<util/verbose.h>
#include <stdlib.h>
CPS_START_NAMESPACE

void* pmalloc(int request){
  void* ptr;
  
  ptr =  malloc(request);
  VRB.Pmalloc("","pmalloc(i)","", ptr, request);
  return ptr;
}

void pfree(void* p){
  VRB.Pfree("","pfree(v*)","",p);
  free((char*) p);
}



void pclear(void){};




CPS_END_NAMESPACE
