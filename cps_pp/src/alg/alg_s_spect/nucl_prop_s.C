#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:45 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_s_spect/nucl_prop_s.C,v 1.1.1.1 2003-06-22 13:34:45 mcneile Exp $
//  $Id: nucl_prop_s.C,v 1.1.1.1 2003-06-22 13:34:45 mcneile Exp $
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
//  $RCSfile: nucl_prop_s.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/alg/alg_s_spect/nucl_prop_s.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
// nucl_prop_s.C

CPS_END_NAMESPACE
#include<alg/nucl_prop_s.h>
#include<alg/quark_prop_s.h>
#include<math.h>
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
