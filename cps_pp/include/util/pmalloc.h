#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Declaration of dynamic memory management routines.	

  $Id: pmalloc.h,v 1.2 2003-07-24 16:53:53 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-07-24 16:53:53 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/pmalloc.h,v 1.2 2003-07-24 16:53:53 zs Exp $
// $Id: pmalloc.h,v 1.2 2003-07-24 16:53:53 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.2  2001/06/19 18:13:18  anj
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
//  Revision 1.2  2001/05/25 06:16:09  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: pmalloc.h,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/util/pmalloc.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
#ifndef _pmalloc_h
#define _pmalloc_h             //!< Prevent multiple inclusion


void* pmalloc(int request);

void pfree(void* p);

void pclear(void);


#endif /* !_pmalloc_h */

CPS_END_NAMESPACE
