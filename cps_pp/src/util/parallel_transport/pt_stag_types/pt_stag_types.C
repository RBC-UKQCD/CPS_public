#include <config.h>
CPS_START_NAMESPACE
/*! \file
  \brief  Definition of ParTransStagTypes class constructor and destructor.

  $Id: pt_stag_types.C,v 1.3 2004-05-10 15:26:55 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-05-10 15:26:55 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_stag_types/pt_stag_types.C,v 1.3 2004-05-10 15:26:55 zs Exp $
//  $Id: pt_stag_types.C,v 1.3 2004-05-10 15:26:55 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: pt_stag_types.C,v $
//  $Revision: 1.3 $
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

ParTransStagTypes::ParTransStagTypes(Lattice & latt) :
                                   ParTrans(latt)
{
  cname = "ParTransStagTypes";
  char *fname = "ParTransStagTypes(Lattice&)";
  VRB.Func(cname,fname);
}


//------------------------------------------------------------------
ParTransStagTypes::~ParTransStagTypes() {
  char *fname = "~ParTransStagTypes()";
  VRB.Func(cname,fname);
}

CPS_END_NAMESPACE
