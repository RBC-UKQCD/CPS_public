#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Dynamic memory management routines.	

  $Id: pmalloc.C,v 1.2 2003-07-24 16:53:54 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-07-24 16:53:54 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/pmalloc/pmalloc.C,v 1.2 2003-07-24 16:53:54 zs Exp $
//  $Id: pmalloc.C,v 1.2 2003-07-24 16:53:54 zs Exp $
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
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/pmalloc/pmalloc.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
CPS_END_NAMESPACE
#include <util/pmalloc.h>
#include <util/verbose.h>
#include <stdlib.h>
CPS_START_NAMESPACE

//! Allocate memory
/*!
  \param request The amount of memory (in bytes) to allocate
  \return A pointer to the allocated memory
*/
void* pmalloc(int request){
  void* ptr;
  
  ptr =  malloc(request);
  VRB.Pmalloc("","pmalloc(i)","", ptr, request);
  return ptr;
}

//! Free allocate memory
/*!
  \param p Pointer to the memory to be freed.
*/
void pfree(void* p){
  VRB.Pfree("","pfree(v*)","",p);
  free((char*) p);
}

void pclear(void){};





CPS_END_NAMESPACE
