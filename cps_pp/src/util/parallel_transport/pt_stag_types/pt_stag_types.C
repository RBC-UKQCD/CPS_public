#include <config.h>
CPS_START_NAMESPACE
/*! \file
  \brief  Definition of ParTransStagTypes class constructor and destructor.

  $Id: pt_stag_types.C,v 1.6 2004-08-18 11:58:07 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:58:07 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_stag_types/pt_stag_types.C,v 1.6 2004-08-18 11:58:07 zs Exp $
//  $Id: pt_stag_types.C,v 1.6 2004-08-18 11:58:07 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: pt_stag_types.C,v $
//  $Revision: 1.6 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_stag_types/pt_stag_types.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/pt.h>
#include <util/lattice.h>
#include <util/verbose.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
/*!
  \param latt The lattice on which this operation is defined
 */
//------------------------------------------------------------------

static StrOrdType old_str_ord;
ParTransStagTypes::ParTransStagTypes(Lattice & latt) :
                                   ParTrans(latt)
{
  cname = "ParTransStagTypes";
  char *fname = "ParTransStagTypes(Lattice&)";
  VRB.Func(cname,fname);
  old_str_ord = lat.StrOrd();
  if (old_str_ord != STAG){
    lat.Convert(STAG);
  }
  pt_init(lat);
  pt_init_g();
}


//------------------------------------------------------------------
ParTransStagTypes::~ParTransStagTypes() {
  char *fname = "~ParTransStagTypes()";
  VRB.Func(cname,fname);
  if ( old_str_ord !=STAG){
    lat.Convert(old_str_ord);
  }
  pt_delete_g();
  pt_delete();
}

CPS_END_NAMESPACE
