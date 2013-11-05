#include<config.h>
#ifdef USE_QMP
#include<qmp.h>
#include<util/qcdio.h>
//#include<qalloc.h>
CPS_START_NAMESPACE
//--------------------------------------------------------------------
/*!\file
  \brief  Definitions of communications routines

  $Id: get_data.C,v 1.10 2012-03-26 13:50:11 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2012-03-26 13:50:11 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qmp/scu/get_data.C,v 1.10 2012-03-26 13:50:11 chulwoo Exp $
//  $Id: get_data.C,v 1.10 2012-03-26 13:50:11 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: get_data.C,v $
//  $Revision: 1.10 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/comms/qmp/scu/get_data.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
CPS_END_NAMESPACE
#include<util/gjp.h>
//#include<sysfunc_cps.h>
CPS_START_NAMESPACE
#ifndef USE_QMP
#define USE_QMP
#endif

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
static IFloat *tmp_buf=NULL;
static QMP_mem_t *rcv_mem, *send_mem, *tmp_mem;

inline static void allocate_buffer(){
  if(rcv_noncache==NULL) {
    rcv_mem = QMP_allocate_memory(sizeof(Float)*MAX_LENGTH);
    rcv_noncache = (IFloat *)QMP_get_memory_pointer(rcv_mem);
  }
  if(send_noncache==NULL) {
    send_mem = QMP_allocate_memory(sizeof(Float)*MAX_LENGTH);
    send_noncache = (IFloat *)QMP_get_memory_pointer(send_mem);
  }

}

inline static void allocate_tmp(){
  if(tmp_buf==NULL) {
    tmp_mem = QMP_allocate_memory(sizeof(Float)*MAX_LENGTH);
    tmp_buf = (IFloat *)QMP_get_memory_pointer(tmp_mem);
  }
}

inline static void check_length(int len){
    if ( len > MAX_LENGTH){
        QMP_printf("getData::len>MAX_LENGTH(%d)\n",MAX_LENGTH);
        QMP_abort(QMP_BAD_MESSAGE);
    }
    if ( len < 0){
        QMP_printf("getData::len<0\n");
        QMP_abort(QMP_BAD_MESSAGE);
    }
}


static const int MAX_DIR=6;
static inline int DIR(int mu,int sign){ return ((sign+1)/2)*MAX_DIR+mu; }

static void PassData(IFloat *rcv_noncache, IFloat *send_noncache, int len_i, int mu, int sign){


   static int initted=0, msglen[2*MAX_DIR];
   static QMP_msgmem_t sndmem[2*MAX_DIR];
   static QMP_msgmem_t rcvmem[2*MAX_DIR];
   static QMP_msghandle_t sndhandle[2*MAX_DIR];
   static QMP_msghandle_t rcvhandle[2*MAX_DIR];
   if (!initted){
       for(int i = 0;i<2*MAX_DIR;i++) msglen[i]=0;
       initted=1;
   }


  if(gjp_local_axis[mu] == 1) {
//    for(int i=0;i<len_i;i++) rcv_noncache[i] = send_noncache[i];
    memcpy(rcv_noncache,send_noncache,len_i*sizeof(IFloat));
    return;
  }

  int len = len_i + (len_i%2);
  int dir = DIR(mu,sign);

  if (len != msglen[dir]){
  if(!UniqueID())printf("PassData(%p %p %d %d %d)\n",rcv_noncache,send_noncache,len_i,mu,sign);
    if (msglen[dir]>0){ 	// previously allocated
      QMP_free_msghandle(sndhandle[dir]);
      QMP_free_msghandle(rcvhandle[dir]);
      QMP_free_msgmem(sndmem[dir]);
      QMP_free_msgmem(rcvmem[dir]);
      if(!UniqueID())printf("msglen[%d](%d) deleted\n",dir,msglen[dir]); 
    }
    sndmem[dir] = QMP_declare_msgmem(send_noncache, len*sizeof(IFloat));
    rcvmem[dir] = QMP_declare_msgmem(rcv_noncache, len*sizeof(IFloat));
    sndhandle[dir] = QMP_declare_send_relative(sndmem[dir], mu,-sign, 0);
    rcvhandle[dir] = QMP_declare_receive_relative(rcvmem[dir], mu, sign, 0);
    msglen[dir]=len;
  }
    QMP_start(rcvhandle[dir]);
    QMP_start(sndhandle[dir]);
    QMP_status_t send_status = QMP_wait(sndhandle[dir]);
    if (send_status != QMP_SUCCESS) 
      QMP_error("Send failed in PassData: %s\n", QMP_error_string(send_status));
    QMP_status_t rcv_status = QMP_wait(rcvhandle[dir]);
    if (rcv_status != QMP_SUCCESS) 
      QMP_error("Receive failed in PassData: %s\n", QMP_error_string(rcv_status));

}

#if 0
static void getData(IFloat *rcv_buf, IFloat *send_buf, int len, int mu, int sign)
{

  allocate_buffer();
  int i = 0;
//  printf("gjp_local_axis[%d]=%d\n",mu,gjp_local_axis[mu]);
  if(gjp_local_axis[mu] == 0) {
    check_length(len);
//    for(i=0;i<len;i++) send_noncache[i] = send_buf[i];
    memcpy(send_noncache,send_buf,len*sizeof(IFloat));
    PassData(rcv_noncache,send_noncache,len,mu,sign);
//    for(i=0;i<len;i++) rcv_buf[i] = rcv_noncache[i];
    memcpy(rcv_buf,rcv_noncache,len*sizeof(IFloat));
  }
  else {
//    for(int i = 0; i < len; ++i) *rcv_buf++ = *send_buf++;
    memcpy(rcv_buf,send_buf,len*sizeof(IFloat));
  }
}

void getPlusData(IFloat *rcv_buf, IFloat *send_buf, int len, int mu){
  IFloat *rcv_p = rcv_buf;
  IFloat *send_p = send_buf;
  int  len_t = len;
  while (len_t > MAX_LENGTH){
    getData(rcv_p,send_p,MAX_LENGTH,mu,+1);
    len_t -= MAX_LENGTH;
    rcv_p += MAX_LENGTH;
    send_p += MAX_LENGTH;
  }
  getData(rcv_p,send_p,len_t,mu,+1);
}

void getMinusData(IFloat *rcv_buf, IFloat *send_buf, int len, int mu){
  IFloat *rcv_p = rcv_buf;
  IFloat *send_p = send_buf;
  int  len_t = len;
  while (len_t > MAX_LENGTH){
    getData(rcv_p,send_p,MAX_LENGTH,mu,-1);
    len_t -= MAX_LENGTH;
    rcv_p += MAX_LENGTH;
    send_p += MAX_LENGTH;
  }
  getData(rcv_p,send_p,len_t,mu,-1);
}
#else
static void getData(void *rcv_buf, void *snd_buf, size_t size, int mu, int sign)
{
    if(gjp_local_axis[mu] == 1) {
        memcpy(rcv_buf, snd_buf, size);
        return;
    }

    QMP_msgmem_t snd_msg = QMP_declare_msgmem(snd_buf, size);
    QMP_msgmem_t rcv_msg = QMP_declare_msgmem(rcv_buf , size);
    QMP_msghandle_t snd_hnd = QMP_declare_send_relative   (snd_msg, mu, -sign, 0);
    QMP_msghandle_t rcv_hnd = QMP_declare_receive_relative(rcv_msg, mu,  sign, 0);

    QMP_start(snd_hnd);
    QMP_start(rcv_hnd);

    QMP_status_t snd = QMP_wait(snd_hnd);
    QMP_status_t rcv = QMP_wait(rcv_hnd);

    if(snd != QMP_SUCCESS || rcv != QMP_SUCCESS) {
      QMP_error("Send/Receive failed in getData: %s %s\n",
                QMP_error_string(snd),
                QMP_error_string(rcv));
    }

    QMP_free_msghandle(snd_hnd);
    QMP_free_msghandle(rcv_hnd);
    QMP_free_msgmem(snd_msg);
    QMP_free_msgmem(rcv_msg);
}

void getPlusData(IFloat *rcv_buf, IFloat *send_buf, int len, int mu)
{
    getData(rcv_buf, send_buf, sizeof(Float) * len, mu, +1);

  // IFloat *rcv_p = rcv_buf;
  // IFloat *send_p = send_buf;
  // int  len_t = len;
  // while (len_t > MAX_LENGTH){
  //   getData(rcv_p,send_p,MAX_LENGTH,mu,+1);
  //   len_t -= MAX_LENGTH;
  //   rcv_p += MAX_LENGTH;
  //   send_p += MAX_LENGTH;
  // }
  // getData(rcv_p,send_p,len_t,mu,+1);
}

void getMinusData(IFloat *rcv_buf, IFloat *send_buf, int len, int mu)
{
    getData(rcv_buf, send_buf, sizeof(Float) * len, mu, -1);

  // IFloat *rcv_p = rcv_buf;
  // IFloat *send_p = send_buf;
  // int  len_t = len;
  // while (len_t > MAX_LENGTH){
  //   getData(rcv_p,send_p,MAX_LENGTH,mu,-1);
  //   len_t -= MAX_LENGTH;
  //   rcv_p += MAX_LENGTH;
  //   send_p += MAX_LENGTH;
  // }
  // getData(rcv_p,send_p,len_t,mu,-1);
}
#endif

//====================================================================
//*  SUI
//*  used in NLocalProp class
//====================================================================

//const SCUDir pos_dir[] = { SCU_XP, SCU_YP, SCU_ZP, SCU_TP };
//const SCUDir neg_dir[] = { SCU_XM, SCU_YM, SCU_ZM, SCU_TM };

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
  allocate_buffer();
  allocate_tmp();
  check_length(len);
  int i;
  for(i=0;i<len;i++) send_noncache[i] = send_buf[i];

  PassData(tmp_buf,send_noncache,len,mu,-1);
  PassData(rcv_noncache,tmp_buf,len,nu,-1);

  for(i=0;i<len;i++) rcv_buf[i] = rcv_noncache[i];
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
  allocate_buffer();
  allocate_tmp();
  check_length(len);

    int i = (dir+1)%4;
    int j = (dir+2)%4;
    int k = (dir+3)%4;
    int l;

    //--------------------------------------------------------------
    // send_buf --> rcv_buf(as a temporary buffer) 
    //--------------------------------------------------------------
    for(l=0;l<len;l++) send_noncache[l] = send_buf[l];
    PassData(rcv_noncache,send_noncache,len,i,-1);
    //--------------------------------------------------------------
    // rcv_buf --> tmp_buf 
    //--------------------------------------------------------------
    PassData(tmp_buf,rcv_noncache,len,j,-1);
    //--------------------------------------------------------------
    // tmp_buf --> rcv_buf 
    //--------------------------------------------------------------
    PassData(rcv_noncache,tmp_buf,len,k,-1);
    for(l=0;l<len;l++) rcv_buf[l] = rcv_noncache[l];
}

CPS_END_NAMESPACE
#endif
