#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004/08/18 11:57:36 $
//  $Header: /space/cvs/cps/cps++/include/alg/nlocal_prop_s.h,v 1.4 2004/08/18 11:57:36 zs Exp $
//  $Id: nlocal_prop_s.h,v 1.4 2004/08/18 11:57:36 zs Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $Revision: 1.4 $
//  $Source: /space/cvs/cps/cps++/include/alg/nlocal_prop_s.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
// nlocal_prop_s.h

#ifndef INCLUDED_NLOCAL_PROP_S_H
#define INCLUDED_NLOCAL_PROP_S_H

CPS_END_NAMESPACE
#include <alg/hadron_prop_s.h>
#include <alg/s_spect_arg.h>
CPS_START_NAMESPACE

//==============================================================
// DERIVED CLASS for Nonlocal operator measurements :
// including delta in 8 and  8' rep and nucleon in 16 rep
// Note: in this version Degenerate quarks are assumed.
//==============================================================

class NLocalPropS : public HadronPropS {

static char cname[];

Float **qp0;
Float **qp1;
Float **qp2;

Float *buffer0[3];
Float *buffer1[3];
Float *buffer2[3];
	// buffers for the quark propagators of off-node sites
	// 3 is the number of colors, propagators 
	// of each color stored seperately
	// because CG solves for propagator color by color

private:

  inline int isDegQuarks()
  { return (qp0 == qp1) && (qp1 == qp2); }

  void transfer(Float **qp, Float **buffer, int t);

  Complex determinant(Complex *v0, Complex *v1, Complex *v2);
	// det of M = (v1, v2, v3)

  Complex element(Float **p0, Float **p1, Float **p2);
	// E_abc * E_ijk * G_ai(x1) * G_bj(x2) * G_ck(x3)

  int map(int *x);

protected:

  void getNeighbors(int t);

  void localVal(Complex *currp, int *s); 
 
public:

// CTOR
   NLocalPropS(Lattice& lattice, StagNonLocalArg& narg);
	
// DTOR
   ~NLocalPropS();
};

#endif

CPS_END_NAMESPACE
