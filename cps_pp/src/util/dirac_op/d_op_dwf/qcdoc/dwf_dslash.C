#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-07-01 21:26:51 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/qcdoc/dwf_dslash.C,v 1.3 2004-07-01 21:26:51 chulwoo Exp $
//  $Id: dwf_dslash.C,v 1.3 2004-07-01 21:26:51 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: dwf_dslash.C,v $
//  $Revision: 1.3 $
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
#include<util/verbose.h>
#include<util/error.h>
#include<stdio.h>
CPS_START_NAMESPACE


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
  dwf_dslash_4(out, gauge_field, in, cb, dag, dwf_lib_arg);
  int temp_size = 49152;


//  printf("dslash : %e %e\n",out->NormSqNode(temp_size),in->NormSqNode(temp_size));


//------------------------------------------------------------------
// Apply 5th-direction Dslash
//------------------------------------------------------------------
  dwf_dslash_5_plus(out, in, mass, dag, dwf_lib_arg);
//  printf("dslash 5 plus : %e %e\n",out->NormSqNode(temp_size),in->NormSqNode(temp_size));

}

CPS_END_NAMESPACE
