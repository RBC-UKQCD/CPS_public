#include<config.h>
CPS_START_NAMESPACE
// mom_meson_p_s.C
// the main changes compared to meson_prop_s.v1.C 
// concern the addition of momenta in localVal()
//

CPS_END_NAMESPACE
#include <util/qcdio.h>  // for printf  ( really?  where? )
CPS_START_NAMESPACE

CPS_END_NAMESPACE
#include <math.h>
#include <util/mom.h>  
#include <alg/mom_meson_p_s.h>
#include <alg/quark_prop_s.h>
CPS_START_NAMESPACE


char MomMesonPropS::cname[] = "MomMesonPropS";

//--------------------------------------------------------------------
// private functions 
//--------------------------------------------------------------------

Complex MomMesonPropS::traceG1DagG2(Float *G1[], Float *G2[], int offset)
{
  Complex res;
  for(int color = 0; color < 3; ++color) {
    Complex *a = (Complex *) (G1[color]+offset);
    Complex *b = (Complex *) (G2[color]+offset);
    for (int i = 0; i < 3; ++i) {
      res += conj(a[i]) * b[i];
    }
  }
  return res;
}

//--------------------------------------------------------------------
// when dir = T, the following are normal hadron operator. 
// phase factors unclear when dir != T
//--------------------------------------------------------------------

int MomMesonPropS::signSC(int x[])         // SC
{ return ((x[(dir+1)%4]+x[(dir+2)%4]+x[(dir+3)%4]) & 1) == 0 ? 1 : -1; }

int MomMesonPropS::signVT(int x[])         // VT
{
  int res = (x[(dir+1)%4] & 1) == 0 ? 1 : -1;
  res += (x[(dir+2)%4] & 1) == 0 ? 1 : -1;
  res += (x[(dir+3)%4] & 1) == 0 ? 1 : -1;
  return res;
}

int MomMesonPropS::signPV(int x[])         // PV
{
  int i = (dir+1)%4;
  int j = (dir+2)%4;
  int k = (dir+3)%4;

  int res = ((x[i]+x[j]) & 1) == 0 ? 1 : -1;
  res += ((x[j]+x[k]) & 1) == 0 ? 1 : -1;
  res += ((x[i]+x[k]) & 1) == 0 ? 1 : -1;
  return res;
}

//--------------------------------------------------------------------
// protected function
//--------------------------------------------------------------------
void MomMesonPropS::localVal(Complex *currt_p, int *s)
{

  int no_of_mom = n_props/4;
  int offset    = X_OFFSET(s); 

  // take trace and multiply with momentum for each state
  Complex trace = traceG1DagG2(qp0, qp1, offset);
  

  for (int imom=0; imom < no_of_mom; imom++){

    *currt_p++ += trace * MOM.fact(imom,s);             // pseudo scalar pion
    *currt_p++ += trace * signVT(s)*MOM.fact(imom,s);   // vector        rho
    *currt_p++ += trace * signPV(s)*MOM.fact(imom,s);   // pseudo vector rho_2
    *currt_p++ += trace * signSC(s)*MOM.fact(imom,s);   // scalar        pion_2
  }

  return; 
}

MomMesonPropS::MomMesonPropS(Lattice &lattice, StagMomMesonArg& arg)
: HadronPropS(lattice, 4*arg.no_of_momenta, arg.dir, QuarkPropSMng::srcSlice(arg.qid0), 1),
  qp0(QuarkPropSMng::prop(arg.qid0)),
  qp1(QuarkPropSMng::prop(arg.qid1)){}
 
MomMesonPropS::~MomMesonPropS(){}





CPS_END_NAMESPACE
