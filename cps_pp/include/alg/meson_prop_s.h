#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:57:36 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/include/alg/meson_prop_s.h,v 1.4 2004-08-18 11:57:36 zs Exp $
//  $Id: meson_prop_s.h,v 1.4 2004-08-18 11:57:36 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.4 $
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
