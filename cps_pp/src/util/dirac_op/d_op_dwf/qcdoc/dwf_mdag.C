#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-07-15 22:23:05 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/qcdoc/dwf_mdag.C,v 1.5 2004-07-15 22:23:05 chulwoo Exp $
//  $Id: dwf_mdag.C,v 1.5 2004-07-15 22:23:05 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: dwf_mdag.C,v $
//  $Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/qcdoc/dwf_mdag.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
// dwf_mdag.C
//
// dwf_mdag is the dagger of the fermion matrix. 
// The in, out fields are defined on the checkerboard lattice
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include<util/dwf.h>
#include<util/gjp.h>
#include<util/vector.h>
#include<util/verbose.h>
#include<util/error.h>
#include<util/dirac_op.h>
//#include "dwf_internal.h"
CPS_START_NAMESPACE


void  dwf_mdag(Vector *out, 
	       Matrix *gauge_field, 
	       Vector *in, 
	       Float mass, 
	       Dwf *dwf_lib_arg)
{

//------------------------------------------------------------------
// Initializations
//------------------------------------------------------------------
  int f_size = 24 * dwf_lib_arg->vol_4d * dwf_lib_arg->ls / 2;
#if 0
  printf("ls=%d\n",dwf_lib_arg->ls);
#endif
  Float minus_kappa_sq = - dwf_lib_arg->dwf_kappa * dwf_lib_arg->dwf_kappa;
  Vector *frm_tmp2 = (Vector *) dwf_lib_arg->frm_tmp2;

//------------------------------------------------------------------
// Apply Dslash E <- O
//------------------------------------------------------------------
#if 0
{IFloat *tmp = (IFloat *)in;
 for(int i =0;i<24;i++) printf("in[%d]=%e\n",i,*tmp++);
}
#endif
  dwf_dslash(frm_tmp2, gauge_field, in, mass, 1, 1, dwf_lib_arg);
#if 0
{IFloat *tmp = (IFloat *)frm_tmp2;
 for(int i =0;i<24;i++) printf("frm_tmp2[%d]=%e\n",i,*tmp++);
}
#endif


//------------------------------------------------------------------
// Apply Dslash O <- E
//------------------------------------------------------------------
  dwf_dslash(out, gauge_field, frm_tmp2, mass, 0, 1, dwf_lib_arg);
#if 0
{IFloat *tmp = (IFloat *)out;
printf("out[0]=%e\n",*tmp);
}
#endif

//------------------------------------------------------------------
// out = in - dwf_kappa_sq * out
//------------------------------------------------------------------
  //out->FTimesV1PlusV2(minus_kappa_sq, out, in, f_size); 
  vaxpy3(out,&minus_kappa_sq,out,in,f_size/6);
  DiracOp::CGflops+=2*f_size;
#if 0
{IFloat *tmp = (IFloat *)out;
printf("out[0]=%e\n",*tmp);}
#endif
}
CPS_END_NAMESPACE
