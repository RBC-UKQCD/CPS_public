#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: mcneile $
//  $Date: 2003-06-22 13:34:47 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp_comms/scu/get_data.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Id: get_data.C,v 1.1.1.1 2003-06-22 13:34:47 mcneile Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.4  2001/08/16 10:50:03  anj
//  The float->Float changes in the previous version were unworkable on QCDSP.
//  To allow type-flexibility, all references to "float" have been
//  replaced with "IFloat".  This can be undone via a typedef for QCDSP
//  (where Float=rfloat), and on all other machines allows the use of
//  double or float in all cases (i.e. for both Float and IFloat).  The I
//  stands for Internal, as in "for internal use only". Anj
//
//  Revision 1.2  2001/06/19 18:12:16  anj
//  Serious ANSIfication.  Plus, degenerate double64.h files removed.
//  Next version will contain the new nga/include/double64.h.  Also,
//  Makefile.gnutests has been modified to work properly, propagating the
//  choice of C++ compiler and flags all the way down the directory tree.
//  The mpi_scu code has been added under phys/nga, and partially
//  plumbed in.
//
//  Everything has newer dates, due to the way in which this first alteration was handled.
//
//  Anj.
//
//  Revision 1.2  2001/05/25 06:16:03  cvs
//  Added CVS keywords to phys_v4_0_0_preCVS
//
//  $RCSfile: get_data.C,v $
//  $Revision: 1.1.1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdsp_comms/scu/get_data.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
CPS_END_NAMESPACE
#include<comms/scu.h>
#include<comms/nga_reg.h>
#include<util/smalloc.h>
#include<util/gjp.h>
#include <sysfunc.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
// mu = {0,1,2,3,4} corresponds to {x,y,z,t,s}
//------------------------------------------------------------------

void getPlusData(IFloat *rcv_buf, IFloat *send_buf, int len, int mu)
{
  if(gjp_local_axis[mu] == 0) {
    SCUDirArg send(send_buf, gjp_scu_dir[2*mu+1], SCU_SEND, len);
    SCUDirArg rcv(rcv_buf, gjp_scu_dir[2*mu], SCU_REC, len);
    SCUTrans(&send);
    SCUTrans(&rcv);
    SCUTransComplete();
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
    SCUDirArg send(send_buf, gjp_scu_dir[2*mu], SCU_SEND, len);
    SCUDirArg rcv(rcv_buf, gjp_scu_dir[2*mu+1], SCU_REC, len);
    SCUTrans(&send);
    SCUTrans(&rcv);
    SCUTransComplete();
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
