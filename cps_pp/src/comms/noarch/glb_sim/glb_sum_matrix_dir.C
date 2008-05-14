#include<config.h>
#ifndef PARALLEL
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2008-05-14 21:20:52 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/noarch/glb_sim/glb_sum_matrix_dir.C,v 1.5 2008-05-14 21:20:52 chulwoo Exp $
//  $Id: glb_sum_matrix_dir.C,v 1.5 2008-05-14 21:20:52 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/noarch/glb_sim/glb_sum_matrix_dir.C,v $
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

void glb_sum_matrix_dir(Matrix * float_p, int dir)
{
  Matrix tmp;
  tmp = *float_p;
}


CPS_END_NAMESPACE
#endif
