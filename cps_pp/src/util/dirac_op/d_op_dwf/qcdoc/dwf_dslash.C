#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:57:50 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/qcdoc/dwf_dslash.C,v 1.6 2004-08-18 11:57:50 zs Exp $
//  $Id: dwf_dslash.C,v 1.6 2004-08-18 11:57:50 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: dwf_dslash.C,v $
//  $Revision: 1.6 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/qcdoc/dwf_dslash.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// dwf_dslash.C
//
// dwf_dslash is the derivative part of the fermion matrix. 
// The in, out fields are defined on the full lattice
// dag = 0/1 <--> Dslash/Dslash^dagger is calculated.
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include<util/dwf.h>
#include<util/gjp.h>
#include<util/vector.h>
#include<util/dirac_op.h>
#include<util/time.h>
#include<util/verbose.h>
#include<util/error.h>
#include<stdio.h>
CPS_START_NAMESPACE

#undef PROFILE
void dwf_dslash(Vector *out, 
		Matrix *gauge_field, 
		Vector *in, 
		Float mass,
		int cb, 
		int dag, 
		Dwf *dwf_lib_arg)
{
//------------------------------------------------------------------
// Apply 4-dimensional Dslash
//------------------------------------------------------------------
#ifdef PROFILE
  DiracOp::CGflops=0;
  Float time = -dclock();
#endif
//  printf("out=%p gauge_field=%p in=%p\n",out,gauge_field,in);
  dwf_dslash_4(out, gauge_field, in, cb, dag, dwf_lib_arg);

#ifdef PROFILE
  time += dclock();
  print_flops("","dwf_dslash_4()",DiracOp::CGflops,time);
#endif
  int temp_size = 49152;



//  printf("dslash : %e %e\n",out->NormSqNode(temp_size),in->NormSqNode(temp_size));


//------------------------------------------------------------------
// Apply 5th-direction Dslash
//------------------------------------------------------------------
#ifdef PROFILE
  time = -dclock();
  DiracOp::CGflops=0;
#endif
  dwf_dslash_5_plus(out, in, mass, dag, dwf_lib_arg);


#ifdef PROFILE
  time += dclock();
  print_flops("","dwf_dslash_5_plus()",DiracOp::CGflops,time);
#endif
}

CPS_END_NAMESPACE
