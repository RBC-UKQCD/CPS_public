#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:46 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_stag_types/d_op_stag_types.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Id: d_op_stag_types.C,v 1.1.1.1 2003-06-22 13:34:46 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.2  2001/06/19 18:12:46  anj
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
//  Revision 1.2  2001/05/25 06:16:06  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: d_op_stag_types.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_stag_types/d_op_stag_types.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//------------------------------------------------------------------
//
// d_op_stag_types.C
//
// Is derived from DiracOp and is relevant to
// all DiracOp classes with Staggered type fermions 
// These classes are derived from DiracOpStagTypes
//
//------------------------------------------------------------------

CPS_END_NAMESPACE
#include<util/dirac_op.h>
#include<util/lattice.h>
#include<util/gjp.h>
#include<util/verbose.h>
#include<util/error.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------
DiracOpStagTypes::DiracOpStagTypes(Lattice & latt,
				   Vector *f_field_out,
				   Vector *f_field_in,
				   CgArg *arg,
				   CnvFrmType cnv_frm_flg) :
                                   DiracOp(latt, 
				           f_field_out,
				           f_field_in, 
				           arg,
				           cnv_frm_flg)
{
  cname = "DiracOpStagTypes";
  char *fname = "DiracOpStagTypes(L&,V*,V*,CgArg*,CnvFrmType)";
  VRB.Func(cname,fname);
}


//------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------
DiracOpStagTypes::~DiracOpStagTypes() {
  char *fname = "~DiracOpStagTypes()";
  VRB.Func(cname,fname);
}
CPS_END_NAMESPACE
