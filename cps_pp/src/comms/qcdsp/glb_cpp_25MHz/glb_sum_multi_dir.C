#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-06-04 21:14:03 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp/glb_cpp_25MHz/glb_sum_multi_dir.C,v 1.3 2004-06-04 21:14:03 chulwoo Exp $
//  $Id: glb_sum_multi_dir.C,v 1.3 2004-06-04 21:14:03 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: glb_sum_multi_dir.C,v $
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp/glb_cpp_25MHz/glb_sum_multi_dir.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
//--------------------------------------------------------------
// glb_sum_dir
// sum over multiple IFloating point data with single precision
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

#define MAX_NUM_WORDS 40  //set to larger than size of double precision 3x3 matrix
static Double64 transmit_buf[MAX_NUM_WORDS];
static Double64 receive_buf[MAX_NUM_WORDS];
static Double64 gsum_buf[MAX_NUM_WORDS];


void glb_sum_multi_dir(Float * float_p, int dir, int len)
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

  int j;
  
  int blocksize=2*len*sizeof(Float);

  // Sum along dir
  //--------------------------------------------------------------
  for(j=0; j<len; j++){
    *(gsum_buf+j) = *(float_p+j);
    *(transmit_buf+j) = *(gsum_buf+j);
  }


  int itmp;
  for (itmp = 1; itmp < NP[dir]; itmp++) {
    SCUDirArg send(transmit_buf, gjp_scu_dir[2*dir], SCU_SEND, blocksize);
    SCUDirArg rcv(receive_buf, gjp_scu_dir[2*dir+1], SCU_REC, blocksize);

    SCUTrans(&send);
    SCUTrans(&rcv);

    SCUTransComplete();
    for(j=0;j<len;j++){
      *(gsum_buf+j) += *(receive_buf+j);
      *(transmit_buf+j) = *(receive_buf+j);
    }
  }


  // Broadcast the result of node with dir coordinate == 0
  //--------------------------------------------------------------

  if(COOR[dir] != 0) {
    for(j=0;j<len;j++){
      *(gsum_buf+j) = 0;
    }
  }
  
  for(j=0; j<len; j++){
    *(transmit_buf+j) = *(gsum_buf+j);
  }

  for (itmp = 1; itmp < NP[dir]; itmp++) {
    SCUDirArg send(transmit_buf, gjp_scu_dir[2*dir], SCU_SEND, blocksize);
    SCUDirArg rcv(receive_buf, gjp_scu_dir[2*dir+1], SCU_REC, blocksize);

    SCUTrans(&send);
    SCUTrans(&rcv);

    SCUTransComplete();
    for(j=0;j<len;j++){
      *(gsum_buf+j) += *(receive_buf+j);
      *(transmit_buf+j) = *(receive_buf+j);
    }
  }
 
  for(j=0; j<len; j++){ 
   *(float_p+j) = *(gsum_buf+j);
  }
 
}


CPS_END_NAMESPACE
