#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-01-13 20:39:14 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp/glb_cpp_hdw_nos/glb_sum_matrix_dir.C,v 1.2 2004-01-13 20:39:14 chulwoo Exp $
//  $Id: glb_sum_matrix_dir.C,v 1.2 2004-01-13 20:39:14 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: glb_sum_matrix_dir.C,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp/glb_cpp_hdw_nos/glb_sum_matrix_dir.C,v $
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
#include<util/vector.h>
#include <comms/sysfunc.h>
CPS_START_NAMESPACE



static Matrix transmit_buf;
static Matrix receive_buf;
static Matrix gsum_buf;
static IFloat *send_buf ;
static IFloat *rcv_buf ;

static volatile unsigned* dsp_scu_base0x10 =
        (volatile unsigned* )(DSP_SCU_BASE + 0x10);



void glb_sum_matrix_dir(Matrix * matrix_p, int dir)
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
  gsum_buf = *matrix_p;
  transmit_buf = gsum_buf;

  send_buf = (IFloat *) &transmit_buf;
  rcv_buf = (IFloat *) &receive_buf;

  int blocksize=sizeof(Float)*sizeof(Matrix);
  int itmp;
  for (itmp = 1; itmp < NP[dir]; itmp++) {
    bsm(send_buf, blocksize, 0, 1, gjp_scu_wire_map[2*dir], TRANSMIT);
    bsm(rcv_buf, blocksize, 0, 1, gjp_scu_wire_map[2*dir+1], RECEIVE);
    while( *(dsp_scu_base0x10+gjp_scu_wire_map[2*dir+1]) ) ;
    while( *(dsp_scu_base0x10+gjp_scu_wire_map[2*dir]) ) ;

 
    gsum_buf += receive_buf;
    transmit_buf = receive_buf;
  }


  // Broadcast the result of node with dir coordinate == 0
  //--------------------------------------------------------------

  if(COOR[dir] != 0) {
    gsum_buf.ZeroMatrix();
  }
    
  transmit_buf = gsum_buf;
  
  for (itmp = 1; itmp < NP[dir]; itmp++) {
     bsm(send_buf, blocksize, 0, 1, gjp_scu_wire_map[2*dir], TRANSMIT);
    bsm(rcv_buf, blocksize, 0, 1, gjp_scu_wire_map[2*dir+1], RECEIVE);
    while( *(dsp_scu_base0x10+gjp_scu_wire_map[2*dir+1]) ) ;
    while( *(dsp_scu_base0x10+gjp_scu_wire_map[2*dir]) ) ;
 
    gsum_buf += receive_buf;
    transmit_buf = receive_buf;
  }


  *matrix_p = gsum_buf;
}


CPS_END_NAMESPACE
