#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-06-04 21:14:04 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp/glb_cpp_nos/glb_sum.C,v 1.3 2004-06-04 21:14:04 chulwoo Exp $
//  $Id: glb_sum.C,v 1.3 2004-06-04 21:14:04 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: glb_sum.C,v $
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp/glb_cpp_nos/glb_sum.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//--------------------------------------------------------------
// glb_sum
//
// Sum over all nodes 
// {GJP.Xnodes(), GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes()}
//--------------------------------------------------------------
CPS_END_NAMESPACE
#include<comms/glb.h>
#include<comms/scu.h>
#include<util/gjp.h>
#include<comms/double64.h>
#include <comms/sysfunc.h>
CPS_START_NAMESPACE



static Double64 transmit_buf;
static Double64 receive_buf;
static Double64 gsum_buf;
static IFloat *send_buf = (IFloat *) &transmit_buf;
static IFloat *rcv_buf = (IFloat *) &receive_buf;

static volatile unsigned* dsp_scu_base0x10 = 
        (volatile unsigned* )(DSP_SCU_BASE + 0x10);


void glb_sum(Float * float_p)
{
  int NP[4] = {GJP.Xnodes(), GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes()};

  gsum_buf = *float_p;

  for(int i = 0; i < 4; ++i) {

      transmit_buf = gsum_buf;

      for (int itmp = 1; itmp < NP[i]; itmp++) {
	bsm(send_buf, 2, 0, 1, gjp_scu_wire_map[2*i], TRANSMIT);
	bsm(rcv_buf, 2, 0, 1, gjp_scu_wire_map[2*i+1], RECEIVE);
	while( *(dsp_scu_base0x10+gjp_scu_wire_map[2*i+1]) ) ;
	while( *(dsp_scu_base0x10+gjp_scu_wire_map[2*i]) ) ;

        gsum_buf += receive_buf;
        transmit_buf = receive_buf;
      }
  }
  *float_p = gsum_buf;
}
CPS_END_NAMESPACE
