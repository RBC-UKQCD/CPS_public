#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:57:36 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/alg/nucl_prop_s.h,v 1.4 2004-08-18 11:57:36 zs Exp $
//  $Id: nucl_prop_s.h,v 1.4 2004-08-18 11:57:36 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.4 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/alg/nucl_prop_s.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
// nucleon_prop_s.h

#ifndef INCLUDED_NUCLEON_PROP_S_H
#define INCLUDED_NUCLEON_PROP_S_H

CPS_END_NAMESPACE
#include <alg/hadron_prop_s.h>
#include <util/rcomplex.h>
#include <alg/s_spect_arg.h>
CPS_START_NAMESPACE

//-----------------------------------------------------
// FORWARD DECLARATION
//-----------------------------------------------------
class QuarkPropS;

//=====================================================
// DERIVED CLASS for nucleon propagator measurements 
//=====================================================

class NucleonPropS : public HadronPropS {

  static char cname[];

  Float **qp0;
  Float **qp1;
  Float **qp2;

private:

  //-----------------------------------------------------
  // compute determinant of 3x3 matrix of quark propagator
  //-----------------------------------------------------
  Complex det(Float *G[], int offset);
 	// det of M = (G[0]+offset, G[1]+offset, G[2]+offset)

  Complex det(Complex *v0, Complex *v1, Complex *v2);
	// det of M = (v1, v2, v3)

  //-----------------------------------------------------
  // if Quarks are degenerate
  //-----------------------------------------------------
  int isDegenerateQuarks() 
  { return qp0==qp1 && qp1==qp2; }

protected:

  void localVal(Complex *currt_p, int *s);
 
public:

  //-----------------------------------------------------
  // CTOR
  //-----------------------------------------------------
  NucleonPropS(Lattice& lattice, StagNucleonArg &arg);

  //-----------------------------------------------------
  // DTOR
  //-----------------------------------------------------
  ~NucleonPropS();

};

#endif

CPS_END_NAMESPACE
