#include<config.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief  Definitions of communications routines

  $Id: get_data.C,v 1.5 2013-04-05 17:51:13 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-04-05 17:51:13 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/mpi/scu/get_data.C,v 1.5 2013-04-05 17:51:13 chulwoo Exp $
//  $Id: get_data.C,v 1.5 2013-04-05 17:51:13 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: get_data.C,v $
//  $Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/mpi/scu/get_data.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
CPS_END_NAMESPACE
#include<comms/scu.h>
//#include<comms/nga_reg.h>
#include<util/smalloc.h>
#include<util/gjp.h>
#include <comms/sysfunc_cps.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
// mu = {0,1,2,3,4} corresponds to {x,y,z,t,s}
//------------------------------------------------------------------

//------------------------------------------------------------------
/*!
  Gets a contiguous block of floating point data from the neighbouring node
  in a positive direction.

  This function also handles the case where there is no such node,
  \e i.e. the lattice is local in that direction.

  \param rcv_buf A buffer into which the data is copied.
  \param send_buf A buffer from which the data is to be copied.
  \param len The amount of data; the number of floating point numbers.
  \param mu The direction of the transfer, one of {0, 1, 2, 3, 4} corresponding
  to {x, y, z, t, s} respectively.

  \ingroup comms
*/
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


/*!
  Gets a contiguous block of floating point data from the neighbouring
  node in a negative direction.

  This function also handles the case where there is no such node,
  \e i.e. the lattice is local in that direction.

  \param rcv_buf A buffer into which the data is copied.
  \param send_buf A buffer from which the data is to be copied.
  \param len The amount of data; the number of floating point numbers.
  \param mu The direction of the transfer, one of {0, 1, 2, 3, 4} corresponding
  to {x, y, z, t, s} respectively.

  \ingroup comms
*/

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
/*!
  Gets a contiguous block of floating point data from the second-nearest
  neighbouring node in a negative direction, \e i.e.  the node in the
  (-mu, -nu) position relative to this one.

  This function also handles the case where there is no such node
  in either of these directions,
  \e i.e. the lattice is local in that direction.

  \param rcv_buf A buffer into which the data is copied.
  \param send_buf A buffer from which the data is to be copied.
  \param len The amount of data; the number of floating point numbers.
  \param mu The direction of the transfer, one of {0, 1, 2, 3} corresponding
  to {x, y, z, t} respectively.
  \param nu The other direction of the transfer.

  \ingroup comms
*/
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
/*!
  Gets a contiguous block of floating point data from the 
  node in a relative position to this one of (-mu, -nu, -rho)
  where mu, nu and rho are any of the directions {x, y, z, t}.
  
  This function also handles the case where there is no such node
  in any of these directions,
  \e i.e. the lattice is local in that direction.

  \param rcv_buf A buffer into which the data is copied.
  \param send_buf A buffer from which the data is to be copied.
  \param len The amount of data; the number of floating point numbers.
  \param dir The direction perpendicular to all of the directions in transfer
  is to take place, one of {0, 1, 2, 3} corresponding to {x, y, z, t}
  respectively. 

  \ingroup comms
*/
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
