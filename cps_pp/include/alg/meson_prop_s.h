#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-07-24 16:53:53 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/alg/meson_prop_s.h,v 1.2 2003-07-24 16:53:53 zs Exp $
//  $Id: meson_prop_s.h,v 1.2 2003-07-24 16:53:53 zs Exp $
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
//  $RCSfile: meson_prop_s.h,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/alg/meson_prop_s.h,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
// meson_prop_s.h

#ifndef INCLUDED_MESON_PROP_S_H
#define INCLUDED_MESON_PROP_S_H

CPS_END_NAMESPACE
#include <alg/hadron_prop_s.h>
#include <util/rcomplex.h>
#include <alg/s_spect_arg.h>
CPS_START_NAMESPACE

//-------------------------------------------------------------------
// DERIVED CLASS for meson propagators measurements 
//-------------------------------------------------------------------
class MesonPropS : public HadronPropS {

   static char cname[];

   Float **qp0;
   Float **qp1;

   //----------------------------------------------------------------
   // Tr(G1^Dagger * G2)
   //----------------------------------------------------------------
   Complex traceG1DagG2(Float *G1[], Float *G2[], int offset);

   //----------------------------------------------------------------
   // phase factors in staggered fermion meson operators
   // signPS = 1
   //----------------------------------------------------------------
   int signSC(int x[]);
   int signPV(int x[]);
   int signVT(int x[]);

protected:

   //----------------------------------------------------------------
   // currt_p is the ptr to the current memmory to write the 
   // 4 meson (4 complex numbers) propagators of the next slice
   //----------------------------------------------------------------
   void localVal(Complex *currt_p, int *s);

public:

   //----------------------------------------------------------------
   // CTOR & DTOR
   //----------------------------------------------------------------
   MesonPropS(Lattice &lattice, StagMesonArg& arg);

   ~MesonPropS();
};

#endif

CPS_END_NAMESPACE
