#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-06-04 21:13:57 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/alg/aots_s.h,v 1.3 2004-06-04 21:13:57 chulwoo Exp $
//  $Id: aots_s.h,v 1.3 2004-06-04 21:13:57 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.3 $
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
