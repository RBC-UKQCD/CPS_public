#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:57:46 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp/glb_cpp_nos/glb_sum_five.C,v 1.4 2004-08-18 11:57:46 zs Exp $
//  $Id: glb_sum_five.C,v 1.4 2004-08-18 11:57:46 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: glb_sum_five.C,v $
//  $Revision: 1.4 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp/glb_cpp_nos/glb_sum_five.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//--------------------------------------------------------------
//  glb_sum_five
//
// Sum over all nodes of the "virtual" 5-dimensional volume.
// {GJP.Xnodes(), GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes(), GJP.Snodes}
// Relevant for spread-out DWF (GJP.s_nodes not 1) only.
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


void glb_sum_five(Float * float_p)
{
  int NP[5] = {GJP.Xnodes(), GJP.Ynodes(), GJP.Znodes(), GJP.Tnodes(), GJP.Snodes()};


  // Sum over the "virtual" 5-dimensional mesh
  //------------------------------------------------------------
  gsum_buf = *float_p;

  int i;
  for(i = 0; i < 5; ++i) {

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

  // Broadcast the result of node (0,0,0,0,0)
  //------------------------------------------------------------
  if(GJP.XnodeCoor() != 0 || 
     GJP.YnodeCoor() != 0 || 
     GJP.ZnodeCoor() != 0 || 
     GJP.TnodeCoor() != 0 || 
     GJP.SnodeCoor() != 0 ) {
    gsum_buf = 0;
  }

  for(i = 0; i < 5; ++i) {
    
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
