#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-07-24 16:53:53 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/alg/myenum.h,v 1.2 2003-07-24 16:53:53 zs Exp $
//  $Id: myenum.h,v 1.2 2003-07-24 16:53:53 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.5  2002/03/11 22:25:43  anj
//  This should now be the correct, fully merged code from our two versions. Anj
//
//  Revision 1.2.2.1  2002/03/08 16:35:01  anj
//  Checking in the Columbia code branch on tag Columbia4_1_1_test-branch, to be
//  merged with the UKQCD head branch shortly.  Anj
//
//  Revision 1.2  2001/06/19 18:11:30  anj
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
//  Revision 1.2  2001/05/25 06:16:00  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: myenum.h,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/alg/myenum.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
#ifndef INCLUDED_MY_ENUM_H
#define INCLUDED_MY_ENUM_H

enum { VECT_LEN = 6, MATRIX_SIZE = 18 };
enum HadronType { SMESON, SMOMMESON, SNUCLEON, SNONLOCAL};

#endif

CPS_END_NAMESPACE
