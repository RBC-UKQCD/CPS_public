#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2004-01-13 20:39:16 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp/scu_nos/get_data.C,v 1.2 2004-01-13 20:39:16 chulwoo Exp $
//  $Id: get_data.C,v 1.2 2004-01-13 20:39:16 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: get_data.C,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp/scu_nos/get_data.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
CPS_END_NAMESPACE
#include<comms/scu.h>
#include<comms/nga_reg.h>
#include<util/smalloc.h>
#include<util/gjp.h>
#include <comms/sysfunc.h>
CPS_START_NAMESPACE

static volatile unsigned* dsp_scu_base0x10 = 
	(volatile unsigned* )(DSP_SCU_BASE + 0x10);


void getPlusData(IFloat *rcv_buf, IFloat *send_buf, int len, int mu)
{
  if(gjp_local_axis[mu] == 0) {
    bsm(send_buf, len, 0, 1, gjp_scu_wire_map[2*mu+1], TRANSMIT);
    bsm(rcv_buf, len, 0, 1, gjp_scu_wire_map[2*mu], RECEIVE);
    while( *(dsp_scu_base0x10+gjp_scu_wire_map[2*mu]) ) ;
    while( *(dsp_scu_base0x10+gjp_scu_wire_map[2*mu+1]) ) ;
  }
  else {
    for(int i = 0; i < len; ++i) {
      *rcv_buf++ = *send_buf++;
    }
  }
}


void getMinusData(IFloat* rcv_buf, IFloat* send_buf, int len, int mu)
{
  if(gjp_local_axis[mu] == 0) {
    bsm(send_buf, len, 0, 1, gjp_scu_wire_map[2*mu], TRANSMIT);
    bsm(rcv_buf, len, 0, 1, gjp_scu_wire_map[2*mu+1], RECEIVE);
    while( *(dsp_scu_base0x10+gjp_scu_wire_map[2*mu+1]) ) ;
    while( *(dsp_scu_base0x10+gjp_scu_wire_map[2*mu]) ) ;
  }
  else {
    for(int i = 0; i < len; ++i) {
      *rcv_buf++ = *send_buf++;
    }
  }
}


//====================================================================
//*  SUI
//*  used in NLocalProp class
//====================================================================

const SCUDir pos_dir[] = { SCU_XP, SCU_YP, SCU_ZP, SCU_TP };
const SCUDir neg_dir[] = { SCU_XM, SCU_YM, SCU_ZM, SCU_TM };

//-------------------------------------------------------------------
//  get data from any ONE site on the 1-D side of 
//  the spatial cube: 2 times SCU transfer
//-------------------------------------------------------------------
void getMinus2Data(IFloat* rcv_buf, IFloat* send_buf, int len, int mu, int nu)
{
    IFloat *tmp_buf = (IFloat *)smalloc(len*sizeof(IFloat));

    SCUDirArg send1(send_buf, pos_dir[mu], SCU_SEND, len);
    SCUDirArg recv1(tmp_buf, neg_dir[mu], SCU_REC, len);
    SCUTrans(&send1);
    SCUTrans(&recv1);
    SCUTransComplete();

    SCUDirArg send2(tmp_buf, pos_dir[nu], SCU_SEND, len);
    SCUDirArg recv2(rcv_buf, neg_dir[nu], SCU_REC, len);
    SCUTrans(&send2);
    SCUTrans(&recv2);
    SCUTransComplete();

    sfree(tmp_buf);
}

//-------------------------------------------------------------------
//  get data from (-1, -1, -1): with dir being the normal direction
//  orthogonal to this hyperplane
//-------------------------------------------------------------------
void getMinus3Data(IFloat* rcv_buf, IFloat* send_buf, int len, int dir)
{
    IFloat *tmp_buf = (IFloat *)smalloc(len*sizeof(IFloat));

    int i = (dir+1)%4;
    int j = (dir+2)%4;
    int k = (dir+3)%4;

    //--------------------------------------------------------------
    // send_buf --> rcv_buf(as a temporary buffer) 
    //--------------------------------------------------------------
    SCUDirArg send1(send_buf, pos_dir[i], SCU_SEND, len);
    SCUDirArg recv1(rcv_buf, neg_dir[i], SCU_REC, len);
    SCUTrans(&send1);
    SCUTrans(&recv1);
    SCUTransComplete();

    //--------------------------------------------------------------
    // rcv_buf --> tmp_buf 
    //--------------------------------------------------------------
    SCUDirArg send2(rcv_buf, pos_dir[j], SCU_SEND, len);
    SCUDirArg recv2(tmp_buf, neg_dir[j], SCU_REC, len);
    SCUTrans(&send2);
    SCUTrans(&recv2);
    SCUTransComplete();

    //--------------------------------------------------------------
    // tmp_buf --> rcv_buf 
    //--------------------------------------------------------------
    SCUDirArg send3(tmp_buf, pos_dir[k], SCU_SEND, len);
    SCUDirArg recv3(rcv_buf, neg_dir[k], SCU_REC, len);
    SCUTrans(&send3);
    SCUTrans(&recv3);
    SCUTransComplete();

    sfree(tmp_buf);
}
















CPS_END_NAMESPACE
