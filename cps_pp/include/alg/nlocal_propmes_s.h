#ifndef INCLUDED_NLSMESON_PROP_S_H
#define INCLUDED_NLSMESON_PROP_S_H
#include<config.h>
#include <alg/hadron_prop_s.h>
#include <alg/s_spect_arg.h>
CPS_START_NAMESPACE
//==============================================================
// DERIVED CLASS for Nonlocal operator measurements :
// This is for 15 meson operator in staggered fermion.
// Note: in this version Degenerate quarks are assumed.
//==============================================================
class NLSMesonPropS : public HadronPropS { 
static char cname[];
Float **qp0[8];
Float **qp00;
Float **qp01;
Float **qp02;
Float **qp03;
Float **qp04;
Float **qp05;
Float **qp06;
Float **qp07;
Float *buffer0[8][3];
	// buffers for the quark propagators of off-node sites
	// 3 is the number of colors, propagators 
	// of each color stored seperately
	// because CG solves for propagator color by color
private:

  inline int isDegQuarks()
    {return 1;}
  //((qp0 == qp1) && (qp1 == qp2)); }

  void transfer(Float **qp, Float **buffer, int t);

  int x_map(int *x);

  int IsOffNode(int *y);

  int epsilon(int *y);

  int zeta(int *y, int l);

  int eta(int *y, int l);

  Complex traceG1DagG2(Float *G1[], Float *G2[], int offset1,int offset2);

protected:

  void getNeighbors(int t);

  void localVal(Complex *currp, int *s); 
 
public:

// CTOR
   NLSMesonPropS(Lattice& lattice, NLStagMesonArg& narg);
	
// DTOR
   ~NLSMesonPropS();
};


CPS_END_NAMESPACE
#endif
