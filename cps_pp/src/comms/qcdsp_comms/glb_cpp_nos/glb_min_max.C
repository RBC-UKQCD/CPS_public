#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2003-10-31 14:15:33 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp_comms/glb_cpp_nos/glb_min_max.C,v 1.2 2003-10-31 14:15:33 zs Exp $
//  $Id: glb_min_max.C,v 1.2 2003-10-31 14:15:33 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: glb_min_max.C,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp_comms/glb_cpp_nos/glb_min_max.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
/*
 *  glb_min_max.C (C++ version)
 */


CPS_END_NAMESPACE
#include<comms/glb.h>
#include<comms/scu.h>
#include<util/gjp.h>
#include <comms/sysfunc.h>
CPS_START_NAMESPACE

#define max(A, B) ((A) > (B) ? (A) : (B))
#define min(A, B) ((A) < (B) ? (A) : (B))


static Float transmit_buf;
static Float receive_buf;
static Float gsum_buf;
static IFloat *send_buf = (IFloat *) &transmit_buf;
static IFloat *rcv_buf = (IFloat *) &receive_buf;


static volatile unsigned* dsp_scu_base0x10 = 
        (volatile unsigned* )(DSP_SCU_BASE + 0x10);


void glb_max(Float * float_p)
{
  int NP[4] = {GJP.Xnodes(), GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes()};

  gsum_buf = *float_p;

  for(int i = 0; i < 4; ++i) {

      transmit_buf = gsum_buf;

      for (int itmp = 1; itmp < NP[i]; itmp++) {
	bsm(send_buf, sizeof(IFloat), 0, 1, gjp_scu_wire_map[2*i], TRANSMIT);
	bsm(rcv_buf, sizeof(IFloat), 0, 1, gjp_scu_wire_map[2*i+1], RECEIVE);
	while( *(dsp_scu_base0x10+gjp_scu_wire_map[2*i+1]) ) ;
	while( *(dsp_scu_base0x10+gjp_scu_wire_map[2*i]) ) ;

        gsum_buf = max(gsum_buf, receive_buf) ;
        transmit_buf = receive_buf;
      }
  }
  *float_p = gsum_buf;
}


void glb_min(Float * float_p)
{
  int NP[4] = {GJP.Xnodes(), GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes()};

  gsum_buf = *float_p;

  for(int i = 0; i < 4; ++i) {

      transmit_buf = gsum_buf;

      for (int itmp = 1; itmp < NP[i]; itmp++) {
	bsm(send_buf, sizeof(IFloat), 0, 1, gjp_scu_wire_map[2*i], TRANSMIT);
	bsm(rcv_buf, sizeof(IFloat), 0, 1, gjp_scu_wire_map[2*i+1], RECEIVE);
	while( *(dsp_scu_base0x10+gjp_scu_wire_map[2*i+1]) ) ;
	while( *(dsp_scu_base0x10+gjp_scu_wire_map[2*i]) ) ;

        gsum_buf = min(gsum_buf, receive_buf) ;
        transmit_buf = receive_buf;
      }
  }
  *float_p = gsum_buf;
}

CPS_END_NAMESPACE
