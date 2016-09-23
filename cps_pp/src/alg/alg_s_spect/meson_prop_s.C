#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004/08/18 11:57:39 $
//  $Header: /space/cvs/cps/cps++/src/alg/alg_s_spect/meson_prop_s.C,v 1.6 2004/08/18 11:57:39 zs Exp $
//  $Id: meson_prop_s.C,v 1.6 2004/08/18 11:57:39 zs Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: meson_prop_s.C,v $
//  $Revision: 1.6 $
//  $Source: /space/cvs/cps/cps++/src/alg/alg_s_spect/meson_prop_s.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
// meson_prop_s.C

CPS_END_NAMESPACE
#include <math.h>
#include <alg/meson_prop_s.h>
#include <alg/quark_prop_s.h>
CPS_START_NAMESPACE

char MesonPropS::cname[] = "MesonPropS";

//--------------------------------------------------------------------
// private functions 
//--------------------------------------------------------------------

Complex MesonPropS::traceG1DagG2(Float *G1[], Float *G2[], int offset)
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

int MesonPropS::signSC(int x[])         // SC
{ return ((x[(dir+1)%4]+x[(dir+2)%4]+x[(dir+3)%4]) & 1) == 0 ? 1 : -1; }

int MesonPropS::signVT(int x[])         // VT
{
  int res = (x[(dir+1)%4] & 1) == 0 ? 1 : -1;
  res += (x[(dir+2)%4] & 1) == 0 ? 1 : -1;
  res += (x[(dir+3)%4] & 1) == 0 ? 1 : -1;
  return res;
}

int MesonPropS::signPV(int x[])         // PV
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
void MesonPropS::localVal(Complex *currt_p, int *s)
{
  int offset = X_OFFSET(s); 

  Complex trace = traceG1DagG2(qp0, qp1, offset);

  *currt_p++ += trace;
  *currt_p++ += trace * signVT(s);
  *currt_p++ += trace * signPV(s);
  *currt_p   += trace * signSC(s);

  return; 
}

MesonPropS::MesonPropS(Lattice &lattice, StagMesonArg& arg)
: HadronPropS(lattice, 4, arg.dir, QuarkPropSMng::srcSlice(arg.qid0), 1),
  qp0(QuarkPropSMng::prop(arg.qid0)),
  qp1(QuarkPropSMng::prop(arg.qid1)) {}
 
MesonPropS::~MesonPropS(){}


CPS_END_NAMESPACE
