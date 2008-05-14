#include<config.h>
#ifndef PARALLEL
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2008-05-14 21:20:52 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/noarch/glb_sim/glb_sum.C,v 1.5 2008-05-14 21:20:52 chulwoo Exp $
//  $Id: glb_sum.C,v 1.5 2008-05-14 21:20:52 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/noarch/glb_sim/glb_sum.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
/*
 *  glb_sum.C 
 *  Fake simulation routine
 */

CPS_END_NAMESPACE
#include<comms/glb.h>
CPS_START_NAMESPACE


void glb_sum(Float * float_p){
  Float tmp;
  tmp = *float_p;
}

CPS_END_NAMESPACE
#endif
