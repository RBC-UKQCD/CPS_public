#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004/08/18 11:57:40 $
//  $Header: /space/cvs/cps/cps++/src/alg/alg_s_spect/nucl_prop_s.C,v 1.6 2004/08/18 11:57:40 zs Exp $
//  $Id: nucl_prop_s.C,v 1.6 2004/08/18 11:57:40 zs Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $RCSfile: nucl_prop_s.C,v $
//  $Revision: 1.6 $
//  $Source: /space/cvs/cps/cps++/src/alg/alg_s_spect/nucl_prop_s.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
// nucl_prop_s.C

CPS_END_NAMESPACE
#include <alg/nucl_prop_s.h>
#include <alg/quark_prop_s.h>
#include <math.h>

CPS_START_NAMESPACE

char NucleonPropS::cname[] = "NucleonPropS";

//----------------------------------------------------------------
// get the determinant of 3x3 matrix G(x,t; x',t') in color space
//----------------------------------------------------------------
Complex NucleonPropS::det(Float *G[], int offset)
{
  Complex res(0,0);
  Complex *a = (Complex *) (G[0]+offset);	// ptr to 1st column
  Complex *b = (Complex *) (G[1]+offset);	// ptr to 2nd column
  Complex *c = (Complex *) (G[2]+offset);	// ptr to 3rd column

  res = a[0]*b[1]*c[2] + a[2]*b[0]*c[1] + a[1]*b[2]*c[0]
      - (a[0]*b[2]*c[1] + a[1]*b[0]*c[2] + a[2]*b[1]*c[0]);

  return res;
}

//----------------------------------------------------------------
// get determinant for 3x3 matrix with 3 columns: v0, v1, v2
//----------------------------------------------------------------
Complex NucleonPropS::det(Complex *v0, Complex *v1, Complex *v2)
{
  return v0[0] * v1[1] * v2[2] + v0[2] * v1[0] * v2[1] +	
	 v0[1] * v1[2] * v2[0] - v0[1] * v1[0] * v2[2] -	
	 v0[2] * v1[1] * v2[0] - v0[0] * v1[2] * v2[1]; 	
}

//----------------------------------------------------------------
// 
//----------------------------------------------------------------
void NucleonPropS::localVal(Complex *currt_p, int *s)
{
  int offset = X_OFFSET(s); 

  if( isDegenerateQuarks() )
  { 
    *currt_p += det(qp0, offset) * 6;
    return;
  }

  //-------------------------------
  // NON-DEGENERATE Quarks
  //-------------------------------

  *currt_p += det((Complex *)(qp0[0]+offset), 
	    	  (Complex *)(qp1[1]+offset), 
		  (Complex *)(qp2[2]+offset)) +  

  	      det((Complex *)(qp0[2]+offset), 
		  (Complex *)(qp1[0]+offset), 
		  (Complex *)(qp2[1]+offset)) +  

	      det((Complex *)(qp0[1]+offset), 
		  (Complex *)(qp1[2]+offset), 
		  (Complex *)(qp2[0]+offset)) -  

  	      det((Complex *)(qp0[1]+offset), 
		  (Complex *)(qp1[0]+offset), 
		  (Complex *)(qp2[2]+offset)) -  

  	      det((Complex *)(qp0[2]+offset), 
		  (Complex *)(qp1[1]+offset), 
		  (Complex *)(qp2[0]+offset)) -  

  	      det((Complex *)(qp0[0]+offset), 
		  (Complex *)(qp1[2]+offset), 
		  (Complex *)(qp2[1]+offset));   

  return ;
}

NucleonPropS::NucleonPropS(Lattice& lattice, StagNucleonArg &arg)
: HadronPropS(lattice, 1, arg.dir, QuarkPropSMng::srcSlice(arg.qid0), 2),
  qp0(QuarkPropSMng::prop(arg.qid0)),
  qp1(QuarkPropSMng::prop(arg.qid1)),
  qp2(QuarkPropSMng::prop(arg.qid2)) {}
 
NucleonPropS::~NucleonPropS(){}


CPS_END_NAMESPACE
