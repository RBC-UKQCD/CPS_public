#ifdef USE_QMP
/*! \file
  \brief  Definition of parallel transport definitions for QCDOC.
  
  $Id: pt_shift.C,v 1.4 2012-08-02 21:20:01 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2012-08-02 21:20:01 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_base/qmp/pt_shift.C,v 1.4 2012-08-02 21:20:01 chulwoo Exp $
//  $Id: pt_shift.C,v 1.4 2012-08-02 21:20:01 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: pt_shift.C,v $
//  $Revision: 1.4 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_base/qmp/pt_shift.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
#include "asq_data_types.h"
#include "pt_int.h"
#include "pt_qcdoc.h"


//! u[x] = v[x+dir] for n_dir forward or backward directions dir.
void PT::shift_field(IFloat **v, const int *dir, int n_dir,
		     int hop, IFloat **u){
  
  int i, length;
  int wire[n_dir];
  for (i=0; i<n_dir;i++) wire[i] = dir[i];
#ifdef USE_QMP
  QMP_msgmem_t msg_mem_p[20];
  QMP_msghandle_t msg_handle_p[20];
  QMP_msghandle_t multiple;
#else
  SCUDirArgMulti SCUmulti;
  SCUDirArgIR *SCUarg_p[2*n_dir];
#endif
  
  
  int comms=0;
  for (i=0; i<n_dir; i++) 
  if (!local[wire[i]/2]){
#ifndef USE_QMP
    SCUarg_p[2*comms] = SCUarg_mat[hop-1][2*wire[i]];
    SCUarg_p[2*comms+1] = SCUarg_mat[hop-1][2*wire[i]+1];
    SCUarg_p[2*comms+1]->Addr((void *)(v[i]+GAUGE_LEN*set_offset(wire[i], hop)));
#else
    msg_mem_p[2*comms] = QMP_declare_msgmem((void *)rcv_buf[wire[i]], 3*hop*non_local_chi[wire[i]]*VECT_LEN*sizeof(IFloat));
    msg_mem_p[2*comms+1] = QMP_declare_strided_msgmem((void *)(v[i]+GAUGE_LEN*set_offset(wire[i], hop)), (size_t)(3*hop*blklen[wire[i]]), numblk[wire[i]], (ptrdiff_t)(3*stride[wire[i]]+3*blklen[wire[i]]));
    
    msg_handle_p[2*comms] = QMP_declare_receive_relative(msg_mem_p[2*comms], wire[i]/2, 1-2*(wire[i]%2), 0);
    msg_handle_p[2*comms+1] = QMP_declare_send_relative(msg_mem_p[2*comms+1], wire[i]/2, 2*(wire[i]%2)-1, 0);
#endif
   
    comms++;
  }

#ifndef USE_QMP
  if (comms) SCUmulti.Init(SCUarg_p,2*comms);
  if (comms) SCUmulti.SlowStartTrans();
#else
  if(comms) {
    multiple = QMP_declare_multiple(msg_handle_p, 2*comms);
    QMP_start(multiple);
  }
#endif
  
//  SCUmulti.TransComplete();
  
  for (i=0; i<n_dir; i++) {
    length = vol-hop*non_local_chi[wire[i]];
    copy_matrix(u[i],v[i],&length,dest_l[hop-1][wire[i]],
		src_l[hop-1][wire[i]]);
  }
#ifndef USE_QMP
  if (comms) SCUmulti.TransComplete();
#else
  if(comms) {
    QMP_status_t qmp_complete_status = QMP_wait(multiple);
    if (qmp_complete_status != QMP_SUCCESS)
          QMP_error("Send failed in shift_field: %s\n", QMP_error_string(qmp_complete_status));
    QMP_free_msghandle(multiple);
    for(int i = 0; i < 2*comms; i++)
      QMP_free_msgmem(msg_mem_p[i]);
  }
#endif


  for (i=0; i<n_dir; i++) {
    length = hop*non_local_chi[wire[i]];
    copy_matrix(u[i],(IFloat*)rcv_buf[wire[i]],&length,
		dest_nl[hop-1][wire[i]],src_nl[hop-1][wire[i]]);
  }
}

//! u[x] = v[x+dir] for n_dir forward or backward directions dir for vector fields
void PT::shift_field_vec(IFloat **v, const int *dir, int n_dir,
		    int hop, IFloat **u){
#if 0
  printf("PT::shift_field_vec() not implemented\n"); exit(-1);
#else
  //printf("shift_field_vec() called\n");
  int i, length;
  int wire[n_dir];
  for (i=0; i<n_dir;i++) wire[i] = dir[i];

#ifndef USE_QMP
  SCUDirArgMulti SCUmulti;
  SCUDirArgIR *SCUarg_p[2*n_dir];
#else
  QMP_msgmem_t msg_mem_p[20];
  QMP_msghandle_t msg_handle_p[20];
  QMP_msghandle_t multiple;
#endif
  
  int comms=0;
  for (i=0; i<n_dir; i++) 
  if (!local[wire[i]/2]){
#ifndef USE_QMP
    SCUarg_p[2*comms] = SCUarg[hop-1][2*wire[i]];
    SCUarg_p[2*comms+1] = SCUarg[hop-1][2*wire[i]+1];
    SCUarg_p[2*comms+1]->Addr((void *)(v[i]+VECT_LEN*set_offset(wire[i], hop)));
#else
    msg_mem_p[2*comms] = QMP_declare_msgmem((void *)rcv_buf[wire[i]], hop*non_local_chi[wire[i]]*VECT_LEN*sizeof(IFloat));
    msg_mem_p[2*comms+1] = QMP_declare_strided_msgmem((void *)(v[i]+VECT_LEN*set_offset(wire[i], hop)), (size_t)(hop*blklen[wire[i]]), numblk[wire[i]], (ptrdiff_t)(stride[wire[i]]+blklen[wire[i]]));
    
    msg_handle_p[2*comms] = QMP_declare_receive_relative(msg_mem_p[2*comms], wire[i]/2, 1-2*(wire[i]%2), 0);
    msg_handle_p[2*comms+1] = QMP_declare_send_relative(msg_mem_p[2*comms+1], wire[i]/2, 2*(wire[i]%2)-1, 0);
#endif
    comms++;
  }

#ifndef USE_QMP
  if (comms) SCUmulti.Init(SCUarg_p,2*comms);
  if (comms) SCUmulti.SlowStartTrans();
//  SCUmulti.TransComplete();
#else
  if(comms) {
    multiple = QMP_declare_multiple(msg_handle_p, 2*comms);
    QMP_start(multiple);
  }
#endif
  
  for (i=0; i<n_dir; i++) {
    length = vol-hop*non_local_chi[wire[i]];
    copy_vector(u[i],v[i],&length,dest_l[hop-1][wire[i]],
		src_l[hop-1][wire[i]]);
  }

#ifndef USE_QMP
  if (comms) SCUmulti.TransComplete();
#else
  if(comms) {
    QMP_status_t qmp_complete_status = QMP_wait(multiple);
    if (qmp_complete_status != QMP_SUCCESS)
          QMP_error("Send failed in shift_field_vec: %s\n", QMP_error_string(qmp_complete_status));
    QMP_free_msghandle(multiple);
    for(int i = 0; i < 2*comms; i++)
      QMP_free_msgmem(msg_mem_p[i]);
  }
#endif

  for (i=0; i<n_dir; i++) {
    length = hop*non_local_chi[wire[i]];
    copy_vector(u[i],(IFloat*)rcv_buf[wire[i]],&length,
		dest_nl[hop-1][wire[i]],src_nl[hop-1][wire[i]]);
  }
#endif
}

//! u[-/+nu](x) = U_[-/+nu](x) 
void PT::shift_link(IFloat **u, const int *dir, int n_dir){

//  char *fname = "pt_shift_link()";
  int length;
#if 0
printf("sizeof(int)=%d\n",sizeof(int));
printf("sizeof(Float)=%d\n",sizeof(Float));
printf("sizeof(int*)=%d\n",sizeof(int*));
printf("sizeof(Float*)=%d\n",sizeof(Float*));
printf("sizeof(unsigned long)=%d\n",sizeof(unsigned long));
printf("sizeof(gauge_agg )=%d\n",sizeof(gauge_agg));
#endif
  for (int i=0; i<n_dir; i++) {
    
    length = local_chi[dir[i]];
    copy_gauge(u[i],uc_l[dir[i]],&length,dest_l[0][dir[i]]);
    length = non_local_chi[dir[i]];
    copy_gauge(u[i],uc_nl[dir[i]],&length,dest_nl[0][dir[i]]);
    
  }

}
#endif
