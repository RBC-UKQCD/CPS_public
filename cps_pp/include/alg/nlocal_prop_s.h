#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:52 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/alg/nlocal_prop_s.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Id: nlocal_prop_s.h,v 1.1.1.1 2003-06-22 13:34:52 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
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
//  $RCSfile: nlocal_prop_s.h,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/alg/nlocal_prop_s.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
// nlocal_prop_s.h

#ifndef INCLUDED_NLOCAL_PROP_S_H
#define INCLUDED_NLOCAL_PROP_S_H

CPS_END_NAMESPACE
#include<alg/hadron_prop_s.h>
#include<alg/s_spect_arg.h>
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
