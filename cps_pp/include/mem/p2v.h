#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-07-24 16:53:53 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/mem/p2v.h,v 1.2 2003-07-24 16:53:53 zs Exp $
//  $Id: p2v.h,v 1.2 2003-07-24 16:53:53 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.2  2001/06/19 18:11:41  anj
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
//  Revision 1.2  2001/05/25 06:16:01  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: p2v.h,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/mem/p2v.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------//
//
// p2v.h
//
// Physical to virtual memory copy routines.
//
//------------------------------------------------------------------//


#ifndef INCLUDED_P2V_H
#define INCLUDED_P2V_H


//------------------------------------------------------------------//
// Copies into CRAM the Vector utilities assembly code
//------------------------------------------------------------------//
void p2vVector();


//------------------------------------------------------------------//
// Copies into CRAM the Wilson library 
//------------------------------------------------------------------//
void p2vWilsonLib();


//------------------------------------------------------------------//
// Copies into CRAM the Staggered dirac_serial assembly code
//------------------------------------------------------------------//
void p2vStagDs();


//------------------------------------------------------------------//
// Copies into CRAM the Gauge heat bath (ghb) assembly code
//------------------------------------------------------------------//
void p2vGhb();


//------------------------------------------------------------------//
// Copies into CRAM the Matrix Multiplication assembly code
//------------------------------------------------------------------//
void p2vCloverLib();

#endif


CPS_END_NAMESPACE
