#include<config.h>
#include<stdio.h>
#include<qalloc.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief  Definitions of communications routines

  $Id: get_data.C,v 1.3 2004-03-15 15:09:02 cwj Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: cwj $
//  $Date: 2004-03-15 15:09:02 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdoc/scu/get_data.C,v 1.3 2004-03-15 15:09:02 cwj Exp $
//  $Id: get_data.C,v 1.3 2004-03-15 15:09:02 cwj Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Log: not supported by cvs2svn $
//  Revision 1.2.2.2  2004/03/15 08:15:57  cwj
//  Update for QOS v2.1
//
//  Revision 1.2.2.1  2004/03/15 00:51:58  cwj
//  Updated communication calls for QOS v2.1
//
//  Revision 1.2  2004/01/13 20:39:07  chulwoo
//  Merging with multibuild
//
//  Revision 1.1.2.1.2.2  2003/11/28 21:44:27  cwj
//   Copy data into non-cachable memory to solve memory coherency problem
//
//  Revision 1.1.2.1.2.1  2003/11/06 20:22:20  cwj
//  *** empty log message ***
//
//  Revision 1.1.1.1  2003/11/04 05:05:03  chulwoo
//
//  starting again
//
//
//  Revision 1.1  2003/10/08 18:38:39  chulwoo
//  start from vanilla_comms
//  start from QCDSP comms
//  start from QCDSP comms
//
//  Revision 1.1.1.1  2003/09/18 22:30:45  chulwoo
//  Mike's files for single node QCDOC + Parallel transport
//  I added some hacks for PARALLEL without MPI_SCU
//  PARALLEL=2 set PARALLEL without MPI_SCU
//
//
//  Revision 1.2  2003/07/24 16:53:54  zs
//  Addition of documentation via doxygen:
//  doxygen-parsable comment blocks added to many source files;
//  New target in makefile and consequent alterations to configure.in;
//  New directories and files under the doc directory.
//
//  Revision 1.1.1.1  2003/06/22 13:34:47  mcneile
//  This is the cleaned up version of the Columbia Physics System.
//  The directory structure has been changed.
//  The include paths have been updated.
//
//
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
//  $Revision: 1.3 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdoc/scu/get_data.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
CPS_END_NAMESPACE
#include<util/gjp.h>
#include<sysfunc.h>
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
const int MAX_LENGTH = 1024;
static IFloat *rcv_noncache=NULL;
static IFloat *send_noncache=NULL;
void getPlusData(IFloat *rcv_buf, IFloat *send_buf, int len, int mu)
{
  if(rcv_noncache==NULL) rcv_noncache = (IFloat *)qalloc(QNONCACHE|QFAST,sizeof(IFloat)*MAX_LENGTH);
  if(send_noncache==NULL) send_noncache = (IFloat *)qalloc(QNONCACHE|QFAST,sizeof(IFloat)*MAX_LENGTH);
	
  int i = 0;
  if(gjp_local_axis[mu] == 0) {
    if ( len > MAX_LENGTH){
        ERR.General("","getPlusData","len>MAX_LENGTH(%d)\n",MAX_LENGTH);
    }
    for(i=0;i<len;i++) send_noncache[i] = send_buf[i];
    SCUDirArg send(send_noncache, gjp_scu_dir[2*mu+1], SCU_SEND, len*sizeof(IFloat));
    SCUDirArg rcv(rcv_noncache, gjp_scu_dir[2*mu], SCU_REC, len*sizeof(IFloat));
    send.StartTrans();
    rcv.StartTrans();
    send.TransComplete();
    rcv.TransComplete();
    for(i=0;i<len;i++) rcv_buf[i] = rcv_noncache[i];
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
  int i;
  if(rcv_noncache==NULL) rcv_noncache = (IFloat *)qalloc(QNONCACHE|QFAST,sizeof(IFloat)*MAX_LENGTH);
  if(send_noncache==NULL) send_noncache = (IFloat *)qalloc(QNONCACHE|QFAST,sizeof(IFloat)*MAX_LENGTH);
  if(gjp_local_axis[mu] == 0) {
    if ( len > MAX_LENGTH){
        ERR.General("","getMinusData","len>MAX_LENGTH(%d)\n",MAX_LENGTH);
    }
    for(i=0;i<len;i++) send_noncache[i] = send_buf[i];
    SCUDirArg send(send_noncache, gjp_scu_dir[2*mu], SCU_SEND, len*sizeof(IFloat));
    SCUDirArg rcv(rcv_noncache, gjp_scu_dir[2*mu+1], SCU_REC, len*sizeof(IFloat));
    send.StartTrans();
    rcv.StartTrans();
    send.TransComplete();
    rcv.TransComplete();
    for(i=0;i<len;i++) rcv_buf[i] = rcv_noncache[i];
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
    int i;
    IFloat *tmp_buf = (IFloat *)qalloc(QFAST|QNONCACHE,len*sizeof(IFloat));
  if(rcv_noncache==NULL) rcv_noncache = (IFloat *)qalloc(QNONCACHE|QFAST,sizeof(IFloat)*MAX_LENGTH);
  if(send_noncache==NULL) send_noncache = (IFloat *)qalloc(QNONCACHE|QFAST,sizeof(IFloat)*MAX_LENGTH);

    for(i=0;i<len;i++) send_noncache[i] = send_buf[i];
    SCUDirArg send1(send_buf, pos_dir[mu], SCU_SEND, len*sizeof(IFloat));
    SCUDirArg recv1(tmp_buf, neg_dir[mu], SCU_REC, len*sizeof(IFloat));
    send1.StartTrans();
    recv1.StartTrans();
    send1.TransComplete();
    recv1.TransComplete();

    SCUDirArg send2(tmp_buf, pos_dir[nu], SCU_SEND, len*sizeof(IFloat));
    SCUDirArg recv2(rcv_buf, neg_dir[nu], SCU_REC, len*sizeof(IFloat));
    send2.StartTrans();
    recv2.StartTrans();
    send2.TransComplete();
    recv2.TransComplete();
    for(i=0;i<len;i++) rcv_buf[i] = rcv_noncache[i];

    qfree(tmp_buf);
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
    IFloat *tmp_buf = (IFloat *)qalloc(QFAST|QNONCACHE,len*sizeof(IFloat));
  if(rcv_noncache==NULL) rcv_noncache = (IFloat *)qalloc(QNONCACHE|QFAST,sizeof(IFloat)*MAX_LENGTH);
  if(send_noncache==NULL) send_noncache = (IFloat *)qalloc(QNONCACHE|QFAST,sizeof(IFloat)*MAX_LENGTH);

    int i = (dir+1)%4;
    int j = (dir+2)%4;
    int k = (dir+3)%4;
    int l;

    //--------------------------------------------------------------
    // send_buf --> rcv_buf(as a temporary buffer) 
    //--------------------------------------------------------------
    for(l=0;l<len;l++) send_noncache[l] = send_buf[l];
    SCUDirArg send1(send_noncache, pos_dir[i], SCU_SEND, len*sizeof(IFloat) );
    SCUDirArg recv1(rcv_noncache, neg_dir[i], SCU_REC, len*sizeof(IFloat) );
    send1.StartTrans();
    recv1.StartTrans();
    send1.TransComplete();
    recv1.TransComplete();

    //--------------------------------------------------------------
    // rcv_buf --> tmp_buf 
    //--------------------------------------------------------------
    SCUDirArg send2(rcv_noncache, pos_dir[j], SCU_SEND, len*sizeof(IFloat) );
    SCUDirArg recv2(tmp_buf, neg_dir[j], SCU_REC, len*sizeof(IFloat));
    send2.StartTrans();
    recv2.StartTrans();
    send2.TransComplete();
    recv2.TransComplete();

    //--------------------------------------------------------------
    // tmp_buf --> rcv_buf 
    //--------------------------------------------------------------
    SCUDirArg send3(tmp_buf, pos_dir[k], SCU_SEND, len*sizeof(IFloat));
    SCUDirArg recv3(rcv_noncache, neg_dir[k], SCU_REC, len*sizeof(IFloat));
    send3.StartTrans();
    recv3.StartTrans();
    send3.TransComplete();
    recv3.TransComplete();
    for(l=0;l<len;l++) rcv_buf[l] = rcv_noncache[l];

    qfree(tmp_buf);
}






CPS_END_NAMESPACE
