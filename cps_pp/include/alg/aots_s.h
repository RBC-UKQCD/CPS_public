#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:52 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/alg/aots_s.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Id: aots_s.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.2  2001/06/19 18:11:29  anj
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
//  $RCSfile: aots_s.h,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/alg/aots_s.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
// aots_s.h
#ifndef INCLUDED_AOTS_S_H
#define INCLUDED_AOTS_S_H

class Aots {

  static char cname[];

  int start;
  int end;
  int step;
  int current;

public :

  //-----------------------------------------------------------------
  // CTOR
  //-----------------------------------------------------------------
  Aots(int s, int e, int stride);

  //-----------------------------------------------------------------
  // DTOR
  //-----------------------------------------------------------------
  ~Aots() {}

  //-----------------------------------------------------------------
  // advancing operator: prefix 
  //-----------------------------------------------------------------
  Aots & operator++() { current += step; return *this; }

  //-----------------------------------------------------------------
  // conversion operator, combined with operator ++, used in the idiom:
  //  	 for(Aots it(); it; ++it){ ... }
  //-----------------------------------------------------------------
  operator const void *() const { return current <= end ? this : 0; }

  //-----------------------------------------------------------------
  // if you don't like the above, here is the alternatives
  //-----------------------------------------------------------------
  void advance() { current += step; }
  int finish() const { return current > end; }

  //-----------------------------------------------------------------
  // status of this Aots
  //-----------------------------------------------------------------
  int begin() const { return start == current; } 
  int last() const { return (current+step) > end; }
  int slice() const { return current; }

  //-----------------------------------------------------------------
  // number of slices of AOTS: minimum = 1
  //-----------------------------------------------------------------
  int numSlices() { return (current - start)/step + 1; }
};

#endif
CPS_END_NAMESPACE
