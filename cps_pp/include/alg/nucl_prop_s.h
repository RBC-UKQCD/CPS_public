#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:52 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/alg/nucl_prop_s.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Id: nucl_prop_s.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.2  2001/06/19 18:11:31  anj
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
//  $RCSfile: nucl_prop_s.h,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/alg/nucl_prop_s.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
// nucleon_prop_s.h

#ifndef INCLUDED_NUCLEON_PROP_S_H
#define INCLUDED_NUCLEON_PROP_S_H

CPS_END_NAMESPACE
#include<alg/hadron_prop_s.h>
#include<util/rcomplex.h>
#include<alg/s_spect_arg.h>
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
