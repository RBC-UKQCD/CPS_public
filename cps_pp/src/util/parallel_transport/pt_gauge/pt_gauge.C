#include <config.h>
CPS_START_NAMESPACE
/*! \file
  \brief  Definition of ParTransStagTypes class constructor and destructor.

  $Id: pt_gauge.C,v 1.3 2004-08-18 11:58:07 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:58:07 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_gauge/pt_gauge.C,v 1.3 2004-08-18 11:58:07 zs Exp $
//  $Id: pt_gauge.C,v 1.3 2004-08-18 11:58:07 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: pt_gauge.C,v $
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_gauge/pt_gauge.C,v $
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
ParTransGauge::ParTransGauge(Lattice & latt) :
                                   ParTrans(latt)
{
  cname = "ParTransGauge";
  char *fname = "ParTransGauge(Lattice&)";
  VRB.Func(cname,fname);
  if (lat.StrOrd() != WILSON && lat.StrOrd() != CANONICAL ){
    old_str_ord = lat.StrOrd();
    lat.Convert(CANONICAL);
  }
  pt_init(lat);
  pt_init_g();
}


//------------------------------------------------------------------
ParTransGauge::~ParTransGauge() {
  char *fname = "~ParTransGauge()";
  VRB.Func(cname,fname);
  pt_delete_g();
  pt_delete();
  if (old_str_ord != WILSON && old_str_ord != CANONICAL ){
    lat.Convert(old_str_ord);
  }
}

CPS_END_NAMESPACE
