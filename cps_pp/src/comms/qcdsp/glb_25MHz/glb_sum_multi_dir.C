#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:57:43 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp/glb_25MHz/glb_sum_multi_dir.C,v 1.4 2004-08-18 11:57:43 zs Exp $
//  $Id: glb_sum_multi_dir.C,v 1.4 2004-08-18 11:57:43 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.4 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp/glb_25MHz/glb_sum_multi_dir.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
/*
 *  glb_sum_dir.C
 */

CPS_END_NAMESPACE
#include<comms/glb.h>
#include<util/error.h>
CPS_START_NAMESPACE

void glb_sum_multi_dir(Float * float_p, int dir, int len)
{
  char *cname = " ";
  char *fname = "glb_sum_multi_dir(F*,i)";
  ERR.NotImplemented(cname, fname);
}


CPS_END_NAMESPACE
