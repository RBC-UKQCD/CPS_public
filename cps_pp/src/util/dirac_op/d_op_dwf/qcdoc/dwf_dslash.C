#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-07-01 21:22:36 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_dwf/qcdoc/dwf_dslash.C,v 1.2 2004-07-01 21:22:36 chulwoo Exp $
//  $Id: dwf_dslash.C,v 1.2 2004-07-01 21:22:36 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.1.4.1  2004/06/09 04:31:54  chulwoo
//  *** empty log message ***
//
//  Revision 1.1.2.1  2004/05/20 14:36:31  pab
//  Files for three optimised dirac operators.
//  Added some patches to the ReadLattice too.
//
//  Revision 1.1.1.1  2003/06/22 13:34:46  mcneile
//  This is the cleaned up version of the Columbia Physics System.
//  The directory structure has been changed.
//  The include paths have been updated.
//
//
//  Revision 1.2  2001/06/19 18:12:38  anj
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
//  Revision 1.2  2001/05/25 06:16:05  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: dwf_dslash.C,v $
//  $Revision: 1.2 $
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
