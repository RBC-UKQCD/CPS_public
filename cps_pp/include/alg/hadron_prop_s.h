#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2005-03-09 20:35:45 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/alg/hadron_prop_s.h,v 1.7 2005-03-09 20:35:45 chulwoo Exp $
//  $Id: hadron_prop_s.h,v 1.7 2005-03-09 20:35:45 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: hadron_prop_s.h,v $
//  $Revision: 1.7 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/alg/hadron_prop_s.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
// hadron_prop_s.h

#ifndef INCLUDED_HADRON_PROP_S_H
#define INCLUDED_HADRON_PROP_S_H

CPS_END_NAMESPACE
#include <util/vector.h>
#include <util/gjp.h>
#include <alg/myenum.h>
CPS_START_NAMESPACE

/*
//```````````````````````````````````````````````````````````````
// CLASS INHERITANCE:
// 
//		    ----------------
//		   |	HadronPropS  |
//     		    ----------------
//		   ^	  ^	   ^
//		  ||       \\	    \\
//	-----------  -------------  ------------
//	|MesonPropS| |NucleonPropS| |NLocalPropS|
//	-----------  -------------  ------------
//,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
*/
//
//
//---------------------------------------------------------------
// FORWARD CLASS DECLARATION
//---------------------------------------------------------------
class MesonPropS;
class MomMesonPropS;
class NucleonPropS;
class NLocalPropS;
class Lattice;

//---------------------------------------------------------------
// ABSTRACT BASE CLASS 
//---------------------------------------------------------------
class HadronPropS {

  static char cname[];

  Lattice &lat;

  int dir;
	// direction along which the hadrons propagate 

  int slice;
	// the source location

  int stride;
	// stride = 1 if sum over the entire hyperplane
	// stride = 2 if sum over the EVEN coordinates sites

  int n_props;
	// number of kinds of propagator: 
	// 4 -- mesons 
	// 1 -- nucleon
	// 4 -- nlocal

  Complex *prop;
	// ptr to propagator for the current configuration on this node

  enum { CR = 2 };
	// Now prop is stored as complex number; for real part CR = 1

  void collect_prop(HadronType type, Float *sum_buf, 
		    const Float *float_p, int unit, int dir, int t);

protected:

  virtual void localVal(Complex *currt_p, int *s) = 0;

  virtual void getNeighbors(int t);

  int X_OFFSET(const int *x);

public:

  //-------------------------------------------------------------
  // CTOR
  //-------------------------------------------------------------
  HadronPropS(Lattice &lattice, int num, int direction, int srcslice, int incr); 

  //-------------------------------------------------------------
  // DTOR
  //-------------------------------------------------------------
  virtual ~HadronPropS();

  //-------------------------------------------------------------
  // Compute Hadron Propagator 
  //-------------------------------------------------------------
  void getHadronPropS();

  //-------------------------------------------------------------
  // down load hadron propagator and do AOTS
  //-------------------------------------------------------------
  void download_prop(HadronType type, Float *buf);

  //-------------------------------------------------------------
  // ACCESSORS
  //-------------------------------------------------------------
  const Float *propPtr() const { return (Float *)prop; }
  
  // return the boundary condition in the hadron propagation dir
  BndCndType bcd();

  //-------------------------------------------------------------
  // return the number of Float's in propagators on one node
  //-------------------------------------------------------------
  int propLen() const; 

  //-------------------------------------------------------------
  // return the number of Float's in propagators on whole lattice
  //-------------------------------------------------------------
  int propLenTotal() const; 

  //-------------------------------------------------------------
  // return the number of Float's in propagators on whole lattice
  // after folding
  //-------------------------------------------------------------
  int propLenFold() const; 

  friend class NucleonPropS;
  friend class MesonPropS;
  friend class MomMesonPropS;
  friend class NLocalPropS;
};

#endif

CPS_END_NAMESPACE
