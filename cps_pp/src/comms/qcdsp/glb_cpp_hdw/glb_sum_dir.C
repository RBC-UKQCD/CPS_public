#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: zs $
//  $Date: 2004-08-18 11:57:46 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp/glb_cpp_hdw/glb_sum_dir.C,v 1.4 2004-08-18 11:57:46 zs Exp $
//  $Id: glb_sum_dir.C,v 1.4 2004-08-18 11:57:46 zs Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: glb_sum_dir.C,v $
//  $Revision: 1.4 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp/glb_cpp_hdw/glb_sum_dir.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//--------------------------------------------------------------
// glb_sum_dir
//
// Sum over all nodes along a direction
// (0,1,2,3,4) <-> (x,y,z,t,s)
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


void glb_sum_dir(Float * float_p, int dir)
{
  int NP[5] = {GJP.Xnodes(), 
	       GJP.Ynodes(), 
	       GJP.Znodes(), 
	       GJP.Tnodes(), 
	       GJP.Snodes()};

  int COOR[5] = {GJP.XnodeCoor(), 
		 GJP.YnodeCoor(), 
		 GJP.ZnodeCoor(), 
		 GJP.TnodeCoor(), 
		 GJP.SnodeCoor()}; 

  // Sum along dir
  //--------------------------------------------------------------
  gsum_buf = *float_p;

  transmit_buf = gsum_buf;

  int itmp;
  for (itmp = 1; itmp < NP[dir]; itmp++) {
    SCUDirArg send(&transmit_buf, gjp_scu_dir[2*dir], SCU_SEND, 2);
    SCUDirArg rcv(&receive_buf, gjp_scu_dir[2*dir+1], SCU_REC, 2);

    SCUTrans(&send);
    SCUTrans(&rcv);

    SCUTransComplete();

    gsum_buf += receive_buf;
    transmit_buf = receive_buf;
  }


  // Broadcast the result of node with dir coordinate == 0
  //--------------------------------------------------------------

  if(COOR[dir] != 0) {
    gsum_buf = 0;
  }
    
  transmit_buf = gsum_buf;
  
  for (itmp = 1; itmp < NP[dir]; itmp++) {
    SCUDirArg send(&transmit_buf, gjp_scu_dir[2*dir], SCU_SEND, 2);
    SCUDirArg rcv(&receive_buf, gjp_scu_dir[2*dir+1], SCU_REC, 2);
    
    SCUTrans(&send);
    SCUTrans(&rcv);
    
    SCUTransComplete();
    
    gsum_buf += receive_buf;
    transmit_buf = receive_buf;
  }


  *float_p = gsum_buf;
}


CPS_END_NAMESPACE
