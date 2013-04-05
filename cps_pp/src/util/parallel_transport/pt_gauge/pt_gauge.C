#include <config.h>
CPS_START_NAMESPACE
/*! \file
  \brief  Definition of ParTransStagTypes class constructor and destructor.

  $Id: pt_gauge.C,v 1.5 2013-04-05 17:46:31 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-04-05 17:46:31 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_gauge/pt_gauge.C,v 1.5 2013-04-05 17:46:31 chulwoo Exp $
//  $Id: pt_gauge.C,v 1.5 2013-04-05 17:46:31 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: pt_gauge.C,v $
//  $Revision: 1.5 $
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
  if (lat.StrOrd() != WILSON && lat.StrOrd() != CANONICAL &&
      lat.StrOrd() != DWF_4D_EOPREC  && lat.StrOrd() != DWF_4D_EOPREC_EE ){
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
  if (lat.StrOrd() != WILSON && lat.StrOrd() != CANONICAL &&
      lat.StrOrd() != DWF_4D_EOPREC && lat.StrOrd() != DWF_4D_EOPREC_EE  ){
    lat.Convert(old_str_ord);
  }
}

CPS_END_NAMESPACE
