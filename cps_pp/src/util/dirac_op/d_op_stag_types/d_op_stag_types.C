#include <config.h>
CPS_START_NAMESPACE
/*! \file
  \brief  Definition of DiracOpStagTypes class constructor and destructor.

  $Id: d_op_stag_types.C,v 1.4 2004-06-04 21:14:07 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-06-04 21:14:07 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_stag_types/d_op_stag_types.C,v 1.4 2004-06-04 21:14:07 chulwoo Exp $
//  $Id: d_op_stag_types.C,v 1.4 2004-06-04 21:14:07 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: d_op_stag_types.C,v $
//  $Revision: 1.4 $
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
