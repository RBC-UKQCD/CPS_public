#include <config.h>
CPS_START_NAMESPACE
/*! \file
  \brief  Definition of DiracOpStagTypes class constructor and destructor.

  $Id: d_op_stag_types.C,v 1.2 2003-07-24 16:53:54 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-07-24 16:53:54 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_stag_types/d_op_stag_types.C,v 1.2 2003-07-24 16:53:54 zs Exp $
//  $Id: d_op_stag_types.C,v 1.2 2003-07-24 16:53:54 zs Exp $
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
//  $Revision: 1.2 $
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
#include <util/dirac_op.h>
#include <util/lattice.h>
#include <util/gjp.h>
#include <util/verbose.h>
#include <util/error.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
/*!
  Only one instance of this class is allowed to be in existence at
  any time.
  \param latt The lattice on which this Dirac operator is defined
  \param f_field_out A (pointer to) a spin-colour field (optionally). 
  \param f_field_in A (pointer to) a spin-colour field (optionally).
  \param arg Parameters for the solver.
  \param cnv_frm_flag Whether the lattice fields should be converted to
  to a new storage order appropriate for the type of fermion action.
  If this is ::CNV_FRM_NO, then just the gauge field is converted.
  If this is ::CNV_FRM_YES, then the fields \a f_field_out and \a f_field_in
  are also converted: This assumes they are initially in the same order as
  the gauge field.
 */
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
/*!
  If the storage order of any fields was changed by the constructor
  then they are changed to the canonical order by the destructor.
*/
//------------------------------------------------------------------
DiracOpStagTypes::~DiracOpStagTypes() {
  char *fname = "~DiracOpStagTypes()";
  VRB.Func(cname,fname);
}

CPS_END_NAMESPACE
