#include <config.h>
CPS_START_NAMESPACE
/*! \file
  \brief  Definition of ParTransStagTypes class constructor and destructor.

  $Id: pt_wilson_t.C,v 1.2 2004-08-09 07:47:26 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-08-09 07:47:26 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_wilson_types/pt_wilson_t.C,v 1.2 2004-08-09 07:47:26 chulwoo Exp $
//  $Id: pt_wilson_t.C,v 1.2 2004-08-09 07:47:26 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: pt_wilson_t.C,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_wilson_types/pt_wilson_t.C,v $
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
ParTransWilsonTypes::ParTransWilsonTypes(Lattice & latt) :
                                   ParTrans(latt)
{
  cname = "ParTransStagTypes";
  char *fname = "ParTransStagTypes(Lattice&)";
  VRB.Func(cname,fname);
  if (lat.StrOrd() != WILSON && lat.StrOrd() != CANONICAL ){
    old_str_ord = lat.StrOrd();
    lat.Convert(CANONICAL);
  }
  pt_init(lat);
  pt_init_g();
}


//------------------------------------------------------------------
ParTransWilsonTypes::~ParTransWilsonTypes() {
  char *fname = "~ParTransStagTypes()";
  VRB.Func(cname,fname);
  pt_delete_g();
  pt_delete();
  if (old_str_ord != WILSON && old_str_ord != CANONICAL ){
    lat.Convert(old_str_ord);
  }
}

CPS_END_NAMESPACE
