#include<config.h>
CPS_START_NAMESPACE
// mom_meson_p_s.h

#ifndef INCLUDED_MOM_MESON_PROP_S_H
#define INCLUDED_MOM_MESON_PROP_S_H

CPS_END_NAMESPACE
#include <alg/hadron_prop_s.h>
#include <util/rcomplex.h>
#include <alg/s_spect_arg.h>
#include <util/mom.h>
CPS_START_NAMESPACE

//-------------------------------------------------------------------
// DERIVED CLASS for meson propagators measurements 
//-------------------------------------------------------------------
class MomMesonPropS : public HadronPropS {

   static char cname[];

   Float **qp0;
   Float **qp1;

   // Mom & mom;
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
   // (4 x number of momenta) mesons ==
   // (4 x no_of_mom) complex numbers propagators of the next slice
   //----------------------------------------------------------------
   void localVal(Complex *currt_p, int *s);

public:

   //----------------------------------------------------------------
   // CTOR & DTOR
   //----------------------------------------------------------------
   MomMesonPropS(Lattice &lattice, StagMomMesonArg& arg);

   ~MomMesonPropS();
};

#endif

CPS_END_NAMESPACE
