#ifdef USE_QMP
/*! \file
  \brief  Definition of parallel transport definitions for QCDOC.
  
  $Id: pt_vec.C,v 1.8 2012-08-02 21:20:01 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2012-08-02 21:20:01 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_base/qmp/pt_vec.C,v 1.8 2012-08-02 21:20:01 chulwoo Exp $
//  $Id: pt_vec.C,v 1.8 2012-08-02 21:20:01 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: pt_vec.C,v $
//  $Revision: 1.8 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_base/qmp/pt_vec.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
#include "asq_data_types.h"
#include "pt_int.h"
#include "pt_qcdoc.h"
#ifndef SCIDAC 
#include <config.h>
#include <util/time_cps.h>
#endif

#if 1
void inline PRINT(char *format,...) {}
#else
#define PRINT printf
#endif

static int MAX_DIR=10;


//Parallel transport of a vector defined on single parity sites
//
//Parameters
//
//n - number of directions in which to parallel transport
//vout - Transported vector
//vin - Initial vector
//dir - list of directions in which to transport the vectors
//cb - Parity of the sites where the vectors are defined
//gauge - Pointer to block of gauge fields in STAG order

#ifndef SCIDAC 

//Normal parallel transport with normal gauge fields
#undef PROFILE
void PT::vec_cb(int n, IFloat **vout, IFloat **vin, const int *dir, int
 parity)
{
  vec_cb_norm(n,vout,vin,dir,parity,gauge_field_addr);
}

//Normal parallel transport, but with user-specified gauge fields
void PT::vec_cb(int n, IFloat **vout, IFloat **vin, const int *dir, int
parity, IFloat * new_gauge_field)
{
  vec_cb_norm(n,vout,vin,dir,parity,new_gauge_field);
}

//Padded parallel transport with normal gauge fields
void PT::vec_cb(int n, IFloat *vout, IFloat **vin, const int *dir, int
parity, int pad)
{
  vec_cb_pad(n,vout,vin,dir,parity,gauge_field_addr);
}

//Padded parallel transport, but with user-specified gauge fields
void PT::vec_cb(int n, IFloat *vout, IFloat **vin, const int *dir, int
parity, int pad, IFloat * new_gauge_field)
{
  vec_cb_pad(n,vout,vin,dir,parity,new_gauge_field);
}

#undef PROFILE
//inline double dclock(){return CPS_NAMESPACE::dclock();}
inline double dclock(){return 0.;}
void PT::vec_cb_norm(int n, IFloat **vout, IFloat **vin, const int *dir,int parity, IFloat * gauge)
{

  //Name our function
//  char *fname="pt_1vec_cb_norm()";
  //List of the different directions
  int wire[MAX_DIR];
  int i;
  int vlen = VECT_LEN;
  int vlen2 = VECT_LEN;

//  Float dtime;
  static int call_num = 0;
  static double setup_time=0.;
  static double total_time=0.;
  static double comp_time=0.;
  static double comm_init_time=0.;
  static double comm_complete_time=0.;

  total_time -=dclock();
  setup_time -=dclock();
  QMP_msgmem_t msg_mem_p[4*NDIM];
  QMP_msghandle_t msg_handle_p[4*NDIM];
//  QMP_msghandle_t multiple;
  setup_time +=dclock();

  call_num++;
  
  
  //Set the transfer directions
  //If wire[i] is even, then we have communication in the negative direction
  //If wire[i] is odd, then we have communication in the positive direction
  for(i=0;i<n;i++)
    wire[i]=dir[i];



  //If wire[i] is odd, then we have parallel transport in the
  //positive direction.  In this case, the matrix multiplication is
  //done before the field is transferred over to the adjacent node
  //
  //If we have transfer in the negative T direction (wire[i] = 6) then
  //we have to copy the appropriate fields into the send buffer

  comp_time -= dclock();
  if(conjugated)
    for(i=0;i<n;i++)
      {
	if(!local[wire[i]/2])
	  {
	    if(wire[i]%2)
	      {
		//printf("dir = %d, pre-mulitply\n", wire[i]);

  partrans_cmv(non_local_chi_cb[wire[i]]/2,uc_nl_cb_pre[parity][wire[i]/2],gauge,vin[i],snd_buf_cb[wire[i]/2]);

	      }
	    else if((wire[i] == 6))
	      {
            partrans_copy_buffer(non_local_chi_cb[6],vin[i],snd_buf_t_cb,Toffset[parity]);
	      }
	  }
      }
  else
    for(i=0;i<n;i++)
      {
	if(!local[wire[i]/2])
	  {
	    if(wire[i]%2)
	      {
 
  partrans_cmv_dag(non_local_chi_cb[wire[i]]/2,uc_nl_cb_pre[parity][wire[i]/2],gauge,vin[i],snd_buf_cb[wire[i]/2]);	 

	      }
	    else if(wire[i] == 6)
	      {
            partrans_copy_buffer(non_local_chi_cb[6],vin[i],snd_buf_t_cb,Toffset[parity]);
	      }
	  }
      }
  comp_time +=dclock();

  setup_time  -= dclock();
  int comms = 0;
  for(i=0;i<n;i++)
    {
      if(!local[wire[i]/2])
	{
	  //Calculate the starting address for the data to be sent
	  IFloat *addr = vin[i] + VECT_LEN * offset_cb[wire[i]];
      msg_mem_p[2*comms]=msg_mem_p[2*comms+1]=NULL;

      msg_handle_p[2*comms] = msg_handle_cb[wire[i]*2];

	  //Initialize the msg_mem for sends
	  if(wire[i]%2) {
        msg_handle_p[2*comms+1] = msg_handle_cb[wire[i]*2+1];
	  } else if(wire[i] == 6) {
        msg_handle_p[2*comms+1] = msg_handle_cb[wire[i]*2+1];
	  } else{
	    msg_mem_p[2*comms+1] = QMP_declare_strided_msgmem((void *)addr, (size_t)blklen_cb[wire[i]], numblk_cb[wire[i]], (ptrdiff_t)(stride_cb[wire[i]]+blklen_cb[wire[i]]));
	    msg_handle_p[2*comms+1] = QMP_declare_send_relative(msg_mem_p[2*comms+1], wire[i]/2, 2*(wire[i]%2)-1, 0);
      }

	  comms++;
	}
    }
//  if(comms) 
//    multiple = QMP_declare_multiple(msg_handle_p, 2*comms);
  setup_time +=dclock();
  

  comm_init_time -=dclock();
  if(comms) 
  for(i=0;i<2*comms;i++)
    QMP_start(msg_handle_p[i]);
    
//    QMP_start(multiple);
  comm_init_time +=dclock();

  comp_time -=dclock();
  //Do local calculations
  PRINT("//Do local calculations\n");
  if(conjugated)
    {
      for(i=0;i<n;i++)
	{
PRINT("wire=%d\n",i);
	
  if(wire[i]%2)
    partrans_cmv(local_chi_cb[wire[i]]/2,uc_l_cb[parity][wire[i]],gauge,vin[i],vout[i]);
  else
    partrans_cmv_dag(local_chi_cb[wire[i]]/2,uc_l_cb[parity][wire[i]],gauge,vin[i],vout[i]);
  
	}
    }
  else
    {
    for(i=0;i<n;i++)
      {
PRINT("wire=%d\n",i);

  if(!(wire[i]%2))
    partrans_cmv(local_chi_cb[wire[i]]/2,uc_l_cb[parity][wire[i]],gauge,vin[i],vout[i]);
  else
    partrans_cmv_dag(local_chi_cb[wire[i]]/2,uc_l_cb[parity][wire[i]],gauge,vin[i],vout[i]);

      }
    }
  comp_time +=dclock();
  
  //End transmission
  if(comms) {
    comm_complete_time -=dclock();
    for(i=0;i<2*comms;i++){
      QMP_status_t qmp_complete_status = QMP_wait(msg_handle_p[i]);
      if (qmp_complete_status != QMP_SUCCESS)
        QMP_error("Send failed in vec_cb_norm: %s\n", QMP_error_string(qmp_complete_status));
    }
    comm_complete_time +=dclock();
    setup_time -=dclock();
//    QMP_free_msghandle(multiple);
    for(int i = 0; i < 2*comms; i++)
      if(msg_mem_p[i]!=NULL){
        QMP_free_msgmem(msg_mem_p[i]);
        QMP_free_msghandle(msg_handle_p[i]);
      }
    setup_time +=dclock();
  }

  comp_time -=dclock();
  //If wire[i] is even, then we have transport in the negative direction.
  //In this case, the vector field is multiplied by the SU(3) link matrix
  //after all communication is complete
  IFloat *fp0, *fp1;
  for(i=0;i<n;i++)
    {
      if(!local[wire[i]/2])
      	{
	  if(!(wire[i]%2))
	    {
	    if(conjugated)
	      partrans_cmv_dag(non_local_chi_cb[wire[i]]/2,uc_nl_cb[parity][wire[i]],gauge,rcv_buf[wire[i]],vout[i]);
	    else
	      partrans_cmv(non_local_chi_cb[wire[i]]/2,uc_nl_cb[parity][wire[i]],gauge,rcv_buf[wire[i]],vout[i]);

	    }
	  //Otherwise we have parallel transport in the positive direction.
	  //In this case, the received data has already been pre-multiplied
	  //All we need to do is to put the transported field in the correct place
	  else
	    {
              partrans_copy(non_local_chi_cb[wire[i]]/2,uc_nl_cb[parity][wire[i]],rcv_buf[wire[i]],vout[i]);
	    }
	}
    }
  comp_time +=dclock();
  total_time +=dclock();
//  ParTrans::PTflops +=33*n*vol;
#ifdef PROFILE
  if(call_num %100==0){
    double fac = 0.01;
    CPS_NAMESPACE::print_time("vec_cb_norm","total_time",total_time*fac);
    CPS_NAMESPACE::print_time("vec_cb_norm","comp_time",comp_time*fac);
    CPS_NAMESPACE::print_time("vec_cb_norm","setup_time",setup_time*fac);
    CPS_NAMESPACE::print_time("vec_cb_norm","comm_init_time",comm_init_time*fac);
    CPS_NAMESPACE::print_time("vec_cb_norm","comm_complete_time",comm_complete_time*fac);
    call_num=0;
    total_time=comp_time=setup_time=comm_init_time=comm_complete_time=0.;
  }
#endif
}

#define PROFILE
#undef PROFILE

void PT::vec_cb_pad(int n, IFloat *vout, IFloat **vin, const int *dir,int parity, IFloat * gauge)
{
  //List of the different directions
  int wire[MAX_DIR];
  int i;

  QMP_msgmem_t msg_mem_p[4*NDIM];
  QMP_msghandle_t msg_handle_p[4*NDIM];
//  QMP_msghandle_t multiple;
  static int call_num = 0;
#if 0
  printf("gauge=%p parity =%d\n",gauge,parity);
  for(i=0;i<n;i++){
    printf("%d: vin=%p vout=%p\n",i,vin[i],vout);
  }
#endif

  call_num++;
  
  //Name our function
//  char *fname="pt_1vec_cb_pad()";
  //VRB.Func("",fname);
  
  //Set the transfer directions
  //If wire[i] is even, then we have communication in the negative direction
  //If wire[i] is odd, then we have communication in the positive direction
  for(i=0;i<n;i++)
    wire[i]=dir[i];

  Float dtime;

  //If wire[i] is odd, then we have parallel transport in the
  //positive direction.  In this case, the matrix multiplication is
  //done before the field is transferred over to the adjacent node
  //
  //If we have transfer in the negative T direction (wire[i] = 6) then
  //we have to copy the appropriate fields into the send buffer
  for(i=0;i<n;i++)
    {
      if(!local[wire[i]/2])
      	{
	  if(wire[i]%2)
	    {
	      //printf("dir = %d, pre-mulitply \n", wire[i]);
#ifdef PROFILE
  dtime  = - dclock();
#endif

  if(conjugated)
    {
      partrans_cmv(non_local_chi_cb[wire[i]]/2,uc_nl_cb_pre[parity][wire[i]/2],gauge,vin[i],snd_buf_cb[wire[i]/2]);
    }
  else
    partrans_cmv_dag(non_local_chi_cb[wire[i]]/2,uc_nl_cb_pre[parity][wire[i]/2],gauge,vin[i],snd_buf_cb[wire[i]/2]);
  
#ifdef PROFILE
  dtime +=dclock();
//  print_flops("",fname,66*non_local_chi_cb[wire[i]],dtime);
#endif
	    }
	  else if(wire[i] == 6)
	    {
#ifdef PROFILE
  dtime  = - dclock();
#endif
#if 1
            partrans_copy_buffer(non_local_chi_cb[6],vin[i],snd_buf_t_cb,Toffset[parity]);
#else
	      for(int j = 0; j < non_local_chi_cb[6];j++)
		  for(k = 0; k < VECT_LEN;k++)
		    *(snd_buf_t_cb+j*VECT_LEN+k) = *(vin[i] + *(Toffset[parity]+j)+ k);
#endif
#ifdef PROFILE
  dtime +=dclock();
  print_flops(fname,"pt_copy_buffer()",0,dtime);
#endif
	    }
	}
    }

#ifdef PROFILE
  dtime  = - dclock();
#endif
  int comms = 0;
  for(i=0;i<n;i++)
    {
      if(!local[wire[i]/2])
	{
	  //Calculate the starting address for the data to be sent
	  IFloat *addr = vin[i] + VECT_LEN * offset_cb[wire[i]];
      msg_mem_p[2*comms]=msg_mem_p[2*comms+1]=NULL;

      msg_handle_p[2*comms] = msg_handle_cb[wire[i]*2];

	  //Initialize the msg_mem for sends
	  if(wire[i]%2) {
        msg_handle_p[2*comms+1] = msg_handle_cb[wire[i]*2+1];
	  } else if(wire[i] == 6) {
        msg_handle_p[2*comms+1] = msg_handle_cb[wire[i]*2+1];
	  } else{
	    msg_mem_p[2*comms+1] = QMP_declare_strided_msgmem((void *)addr, (size_t)blklen_cb[wire[i]], numblk_cb[wire[i]], (ptrdiff_t)(stride_cb[wire[i]]+blklen_cb[wire[i]]));
	    msg_handle_p[2*comms+1] = QMP_declare_send_relative(msg_mem_p[2*comms+1], wire[i]/2, 2*(wire[i]%2)-1, 0);
      }
	  comms++;
	}
    }
#ifdef PROFILE
  dtime +=dclock();
//  print_flops(fname,"Addr",0,dtime);
#endif

#ifdef PROFILE
  dtime  = - dclock();
#endif

  if(comms) {
//    multiple = QMP_declare_multiple(msg_handle_p, 2*comms);
//    QMP_start(multiple);
  for(i=0;i<2*comms;i++)
    QMP_start(msg_handle_p[i]);
  }

#ifdef PROFILE
  dtime +=dclock();
//  print_flops(fname,"StartTrans()",0,dtime);
#endif


  PRINT("//Do local calculations\n");
  //Do local calculations
  for(i=0;i<n;i++)
    {
#ifdef PROFILE
  dtime  = - dclock();
#endif

PRINT("wire=%d\n",i);
      if((wire[i]%2 && conjugated) || ((wire[i]%2 == 0) && (conjugated == 0)))
	{
	partrans_cmv_pad(local_chi_cb[wire[i]]/2,uc_l_pad_cb[parity][wire[i]],gauge,vin[i],vout);
	}
      else
	partrans_cmv_dag_pad(local_chi_cb[wire[i]]/2,uc_l_pad_cb[parity][wire[i]],gauge,vin[i],vout);
  PRINT("//local done\n");

#ifdef PROFILE
  dtime +=dclock();
//  print_flops("",fname,66*local_chi_cb[wire[i]],dtime);
#endif
    }
#ifdef PROFILE
  dtime  = - dclock();
#endif
  if(comms) {
//    QMP_status_t qmp_complete_status = QMP_wait(multiple);
    for(i=0;i<2*comms;i++){
      QMP_status_t qmp_complete_status = QMP_wait(msg_handle_p[i]);
      if (qmp_complete_status != QMP_SUCCESS)
	  QMP_error("Send failed in vec_cb_pad: status(%d)=%s\n", i,QMP_error_string(qmp_complete_status));
    }
//    QMP_free_msghandle(multiple);
    for(int i = 0; i < 2*comms; i++)
      if(msg_mem_p[i]!=NULL){
        QMP_free_msgmem(msg_mem_p[i]);
        QMP_free_msghandle(msg_handle_p[i]);
      }
  }
#ifdef PROFILE
  dtime +=dclock();
//  print_flops(fname,"TransComplete()",0,dtime);
#endif

  //If wire[i] is even, then we have transport in the negative direction.
  //In this case, the vector field is multiplied by the SU(3) link matrix
  //after all communication is complete
//  IFloat *fp0,*fp1;
  for(i=0;i<n;i++)
    {
      if(!local[wire[i]/2])
	{
	  if(!(wire[i]%2))
	    {
#ifdef PROFILE
  dtime  = - dclock();
#endif
 
  if(conjugated)
    partrans_cmv_dag_pad(non_local_chi_cb[wire[i]]/2,uc_nl_pad_cb[parity][wire[i]],gauge,rcv_buf[wire[i]],vout);
  else
    partrans_cmv_pad(non_local_chi_cb[wire[i]]/2,uc_nl_pad_cb[parity][wire[i]],gauge,rcv_buf[wire[i]],vout);

#ifdef PROFILE
  dtime +=dclock();
//  print_flops("",fname,66*non_local_chi_cb[wire[i]],dtime);
#endif
	    }
	  //Otherwise we have parallel transport in the positive direction.
	  //In this case, the received data has already been pre-multiplied
	  //All we need to do is to put the transported field in the correct place
	  else
	    {
#ifdef PROFILE
  dtime  = - dclock();
#endif
              partrans_copy_pad(non_local_chi_cb[wire[i]]/2,uc_nl_pad_cb[parity][wire[i]],rcv_buf[wire[i]],vout);
#ifdef PROFILE
  dtime +=dclock();
//  print_flops(fname,"pt_copy_pad()",0,dtime);
#endif
	    }
	}
    }
//  ParTrans::PTflops +=33*n*vol;
}

#endif // #ifndef SCIDAC

#undef PROFILE
//Parallel transport of a vector through one hop
void PT::vec(int n, IFloat **vout, IFloat **vin, const int *dir){
  int i;
  static int call_num=0;
  #ifdef USE_QMP
  QMP_msgmem_t *msg_mem_p = (QMP_msgmem_t *)Alloc("","vec_cb_norm", "msg_mem_p", 2*non_local_dirs*sizeof(QMP_msgmem_t));
  QMP_msghandle_t* msg_handle_p = (QMP_msghandle_t *)Alloc("","vec_cb_norm", "msg_handle_p", 2*non_local_dirs*sizeof(QMP_msghandle_t));
  QMP_msghandle_t multiple;
  #else
  SCUDirArgIR *SCUarg_p[2*n];
  SCUDirArgMulti SCUmulti;
  #endif
  call_num++;

#ifdef PROFILE
  Float dtime  = - dclock();
#endif
  int wire[MAX_DIR];

  char *fname="pt_1vec";
//  VRB.Func("",fname);
	
  int non_local_dir=0;
  for(i=0;i<n;i++) wire[i] = dir[i]; // from (x,y,z,t) to (t,x,y,z)
//  for(i=0;i<n;i++) printf("wire[%d]=%d\n",i,dir[i]);
  for(i=0;i<n;i++)
  if (!local[wire[i]/2]){
    IFloat * addr = (vin[i]+VECT_LEN*offset[wire[i]]);
    msg_mem_p[2*non_local_dir] = QMP_declare_msgmem((void *)rcv_buf[wire[i]], non_local_chi[wire[i]]*VECT_LEN*sizeof(IFloat));
    msg_mem_p[2*non_local_dir+1] = QMP_declare_strided_msgmem((void *)addr, (size_t)blklen[wire[i]], numblk[wire[i]], (ptrdiff_t)(stride[wire[i]]+blklen[wire[i]]));
    
    msg_handle_p[2*non_local_dir] = QMP_declare_receive_relative(msg_mem_p[2*non_local_dir], wire[i]/2, 1-2*(wire[i]%2), 0);
    msg_handle_p[2*non_local_dir+1] = QMP_declare_send_relative(msg_mem_p[2*non_local_dir+1], wire[i]/2, 2*(wire[i]%2)-1, 0);
    non_local_dir++;
  }

  if(non_local_dir) {
    multiple = QMP_declare_multiple(msg_handle_p, 2*non_local_dir);
    QMP_start(multiple);
  }
	
  for(i=0;i<n;i++) 
    partrans_cmv_agg(local_chi[wire[i]],(long)uc_l[wire[i]], (long)vin[i],(long)vout[i]);
	
  if(non_local_dir) {
    QMP_status_t qmp_complete_status = QMP_wait(multiple);
    if (qmp_complete_status != QMP_SUCCESS)
	  QMP_error("Send failed in vec_cb_norm: %s\n", QMP_error_string(qmp_complete_status));
    QMP_free_msghandle(multiple);
    for(int i = 0; i < 2*non_local_dir; i++)
      QMP_free_msgmem(msg_mem_p[i]);
    Free(msg_handle_p);
    Free(msg_mem_p);
  }

  for(i=0;i<n;i++) 
    partrans_cmv_agg(non_local_chi[wire[i]],(long)uc_nl[wire[i]], (long)rcv_buf[wire[i]],(long)vout[i]);

#ifdef PROFILE
  dtime +=dclock();
  print_flops("",fname,66*n*vol,dtime);
#endif
  Flops +=66*n*vol;
}
#endif
