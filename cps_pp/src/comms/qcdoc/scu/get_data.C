#include<config.h>
#include<util/qcdio.h>
#include<qalloc.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief  Definitions of communications routines

  $Id: get_data.C,v 1.15 2008-11-05 21:22:42 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2008-11-05 21:22:42 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qcdoc/scu/get_data.C,v 1.15 2008-11-05 21:22:42 chulwoo Exp $
//  $Id: get_data.C,v 1.15 2008-11-05 21:22:42 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: get_data.C,v $
//  $Revision: 1.15 $
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
const int MAX_LENGTH = 4096;
const int MAX_COPY = 128;
static IFloat *rcv_noncache=NULL;
static IFloat *send_noncache=NULL;

void get1Data(IFloat *rcv_buf, IFloat *send_buf, int len, int mu, int plus)
{
  char *fname = "get1Data";
//  printf("get1Data(%p %p %d %d %d\n\n",rcv_buf,send_buf,len,mu,plus);
  if(rcv_noncache==NULL) rcv_noncache = (IFloat *)qalloc(QNONCACHE|QFAST,sizeof(IFloat)*MAX_LENGTH);
  if(rcv_noncache==NULL) rcv_noncache = (IFloat *)qalloc(QFAST,sizeof(IFloat)*MAX_LENGTH);
  if(rcv_noncache==NULL)
    ERR.Pointer("",fname,"rcv_noncache");
  if(send_noncache==NULL) send_noncache = (IFloat *)qalloc(QNONCACHE|QFAST,sizeof(IFloat)*MAX_LENGTH);
  if(send_noncache==NULL) send_noncache = (IFloat *)qalloc(QFAST,sizeof(IFloat)*MAX_LENGTH);
  if(send_noncache==NULL)
    ERR.Pointer("",fname,"send_noncache");
//  printf("rcv_buf=%p send_buf=%p len=%d \n",rcv_buf,send_buf,len);
	
  int i = 0;
  while (len >0){
    int copy_len = len;
    if (copy_len >MAX_LENGTH) copy_len = MAX_LENGTH;
    if(gjp_local_axis[mu] == 0) {
      if ( len > MAX_LENGTH){
//          ERR.General("",fname,"len(%d)>MAX_LENGTH(%d)\n",len,MAX_LENGTH);
      }
//      for(i=0;i<copy_len;i++) send_noncache[i] = send_buf[i];
	  memcpy(send_noncache,send_buf,copy_len*sizeof(IFloat));
      SCUDirArgIR send(send_noncache, gjp_scu_dir[2*mu+plus], SCU_SEND, copy_len*sizeof(IFloat));
      SCUDirArgIR rcv(rcv_noncache, gjp_scu_dir[2*mu+(1-plus)], SCU_REC, copy_len*sizeof(IFloat));
	if (qalloc_is_noncache(send_noncache)) send.StartTrans();
	else send.SlowStartTrans();
	if (qalloc_is_noncache(rcv_noncache)) rcv.StartTrans();
	else rcv.SlowStartTrans();

      send.TransComplete();
      rcv.TransComplete();
//      for(i=0;i<copy_len;i++) rcv_buf[i] = rcv_noncache[i];
	  memcpy(rcv_buf,rcv_noncache,copy_len*sizeof(IFloat));
    }
    else {
      for(int i = 0; i < copy_len; ++i) {
//        *rcv_buf++ = *send_buf++;
      }
	  memcpy(rcv_buf,send_buf,copy_len*sizeof(IFloat));
    }
    send_buf += MAX_LENGTH;
    rcv_buf += MAX_LENGTH;
    len = len - MAX_LENGTH;
  }
//  printf("get1Data() done\n");
}

void getPlusData(IFloat *rcv_buf, IFloat *send_buf, int len, int mu){
  get1Data(rcv_buf,send_buf,len,mu,1);
}
void getMinusData(IFloat *rcv_buf, IFloat *send_buf, int len, int mu){
  get1Data(rcv_buf,send_buf,len,mu,0);
}


#if 0
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
//  printf("rcv_buf=%p send_buf=%p len=%d \n",rcv_buf,send_buf,len);
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
#endif


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

static SCUDirArg scu_send;
static SCUDirArg scu_recv;
void getData(IFloat* rcv_buf, int rblklen, int rnumblk,int rstr,
        IFloat* send_buf, int sblklen, int snumblk,int sstr,
        int mu, int sign){
	char *fname = "getData()";
	int total_length;
	if (rblklen !=sblklen)
	ERR.General("",fname,
"send and receive block size different: not implemented yet\n");
	if ((total_length = rblklen*rnumblk) !=sblklen*snumblk)
	ERR.General("",fname, "send and receive total size different\n");
	if(gjp_local_axis[mu] == 0) {
		IFloat *saddr = send_buf;
		IFloat *raddr = rcv_buf;
		int r_i=0, r_j=0;
		for(int i = 0;i< snumblk;i++){
			for(int j = 0;j< sblklen;j++){
				*raddr++= *saddr++;
				r_j++;
				if (r_j>=rblklen){ r_i++;r_j=0; raddr += rstr;}
			}
			saddr +=sstr;
		}
		return;
	}
	int s_dir = 2*mu+(1+sign)/2; printf("%s:signe = %d send = %d\n",fname,sign,s_dir);
	int r_dir = 2*mu+(1-sign)/2;
	int r_copy = 0;

	if (!qalloc_is_communicable(send_buf) ){
		if (total_length > MAX_LENGTH)
		ERR.General("",fname, "too big from non-communicable memory\n");
		if (total_length > MAX_COPY ){
			IFloat *saddr = send_buf;
			IFloat *raddr = send_noncache;
			int r_i=0, r_j=0;
			for(int i = 0;i< snumblk;i++){
				for(int j = 0;j< sblklen;j++){
					*raddr++= *saddr++;
				}
				saddr +=sstr;
			}
    		scu_send.Init(send_noncache, gjp_scu_dir[s_dir], SCU_SEND, 
			total_length*sizeof(IFloat));
		} else 
    	scu_send.Init(send_buf, gjp_scu_dir[s_dir], SCU_SEND,
		sblklen*sizeof(IFloat),sblklen,sstr);
	} 

	if ( !qalloc_is_communicable(rcv_buf) ){
		if (total_length > MAX_LENGTH)
		ERR.General("",fname, "too big from non-communicable memory\n");
		if (total_length > MAX_COPY ){
			r_copy = 1;
    		scu_recv.Init(rcv_noncache, gjp_scu_dir[r_dir], SCU_SEND, 
				total_length*sizeof(IFloat));
		} else 
    		scu_recv.Init(rcv_buf, gjp_scu_dir[r_dir], SCU_SEND, 
			rblklen*sizeof(IFloat),rblklen,rstr);
	} 
	scu_send.StartTrans();
	scu_recv.StartTrans();
	scu_send.TransComplete();
	scu_recv.TransComplete();
	if (r_copy){
			IFloat *saddr = rcv_noncache;
			IFloat *raddr = rcv_buf;
			int r_i=0, r_j=0;
			for(int i = 0;i< rnumblk;i++){
				for(int j = 0;j< rblklen;j++){
					*raddr++= *saddr++;
				}
				raddr +=rstr;
			}
	}
}




CPS_END_NAMESPACE
