#include<config.h>
#ifndef PARALLEL
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2008-05-14 21:20:52 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/noarch/scu_sim/get_data.C,v 1.5 2008-05-14 21:20:52 chulwoo Exp $
//  $Id: get_data.C,v 1.5 2008-05-14 21:20:52 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/noarch/scu_sim/get_data.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
CPS_END_NAMESPACE
#include<comms/scu.h>
CPS_START_NAMESPACE


//------------------------------------------------------------------
// mu = {0,1,2,3,4} corresponds to {x,y,z,t,s}
//------------------------------------------------------------------

void getPlusData(IFloat *rcv_buf, IFloat *send_buf, int len, int mu)
{
  for(int i = 0; i < len; ++i) {
    *rcv_buf++ = *send_buf++;
  }
}


void getMinusData(IFloat* rcv_buf, IFloat* send_buf, int len, int mu)
{
  for(int i = 0; i < len; ++i) {
    *rcv_buf++ = *send_buf++;
  }
}

void getMinus2Data(IFloat* rcv_buf, IFloat* send_buf, int len, int mu, int nu)
{
  for(int i = 0; i < len; ++i) {
    *rcv_buf++ = *send_buf++;
  }
}

void getMinus3Data(IFloat* rcv_buf, IFloat* send_buf, int len, int dir)
{
  for(int i = 0; i < len; ++i) {
    *rcv_buf++ = *send_buf++;
  }
}
CPS_END_NAMESPACE
#endif
