#include <config.h>
CPS_START_NAMESPACE
/*! \file
  \brief  Routine used internally in the DiracOpWilson class.

  $Id: wilson_end.C,v 1.4 2004-08-18 11:57:52 zs Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:57:52 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/noarch/wilson_end.C,v 1.4 2004-08-18 11:57:52 zs Exp $
//  $Id: wilson_end.C,v 1.4 2004-08-18 11:57:52 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.4 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/dirac_op/d_op_wilson/noarch/wilson_end.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/wilson.h>
#include <util/smalloc.h>
#include <util/verbose.h>
CPS_START_NAMESPACE

  
void wilson_end( Wilson *wilson_p)
{
  char *cname = " ";
  char *fname = "wilson_end(Wilson*)";
  VRB.Func(cname,fname);

  VRB.Sfree(cname,fname, "af[1]", wilson_p->af[1]);
  sfree(wilson_p->af[1]);

  VRB.Sfree(cname,fname, "af[0]", wilson_p->af[0]);
  sfree(wilson_p->af[0]);

  VRB.Sfree(cname,fname, "ptr", wilson_p->ptr);
  sfree(wilson_p->ptr);

}

CPS_END_NAMESPACE
