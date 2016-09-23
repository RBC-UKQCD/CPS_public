#include<config.h>
#ifndef USE_QMP
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2012/03/26 13:50:11 $
//  $Header: /space/cvs/cps/cps++/src/comms/noarch/glb_sim/glb_min_max_sim.C,v 1.2 2012/03/26 13:50:11 chulwoo Exp $
//  $Id: glb_min_max_sim.C,v 1.2 2012/03/26 13:50:11 chulwoo Exp $
//  $Name: v5_0_16_hantao_io_test_v7 $
//  $Locker:  $
//  $Revision: 1.2 $
//  $Source: /space/cvs/cps/cps++/src/comms/noarch/glb_sim/glb_min_max_sim.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
/*
 *  glb_min_max.C
 *  Fake simulation routine
 */


CPS_END_NAMESPACE
#include<comms/glb.h>
CPS_START_NAMESPACE


void glb_max(Float * float_p){
  Float tmp;
  tmp = *float_p;
}


void glb_min(Float * float_p){
  Float tmp;
  tmp = *float_p;
}

CPS_END_NAMESPACE
#endif
