/*! \file
  \brief  Definition of parallel transport definitions for QCDOC.
  
  $Id: pt_shift.C,v 1.1 2006-07-05 18:13:49 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2006-07-05 18:13:49 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_base/qcdoc/pt_shift.C,v 1.1 2006-07-05 18:13:49 chulwoo Exp $
//  $Id: pt_shift.C,v 1.1 2006-07-05 18:13:49 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: pt_shift.C,v $
//  $Revision: 1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_base/qcdoc/pt_shift.C,v $
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
  SCUDirArgMulti SCUmulti;
  SCUDirArgIR *SCUarg_p[2*n_dir];
  
  int comms=0;
  for (i=0; i<n_dir; i++) 
  if (!local[wire[i]/2]){
    SCUarg_p[2*comms] = SCUarg_mat[hop-1][2*wire[i]];
    SCUarg_p[2*comms+1] = SCUarg_mat[hop-1][2*wire[i]+1];
    SCUarg_p[2*comms+1]->Addr((void *)(v[i]+GAUGE_LEN*set_offset(wire[i], hop)));
    comms++;
  }

  if (comms) SCUmulti.Init(SCUarg_p,2*comms);
  if (comms) SCUmulti.SlowStartTrans();
//  SCUmulti.TransComplete();
  
  for (i=0; i<n_dir; i++) {
    length = vol-hop*non_local_chi[wire[i]];
    copy_matrix(u[i],v[i],&length,dest_l[hop-1][wire[i]],
		src_l[hop-1][wire[i]]);
  }
  if (comms) SCUmulti.TransComplete();

  for (i=0; i<n_dir; i++) {
    length = hop*non_local_chi[wire[i]];
    copy_matrix(u[i],(IFloat*)rcv_buf[wire[i]],&length,
		dest_nl[hop-1][wire[i]],src_nl[hop-1][wire[i]]);
  }
}

//! u[x] = v[x+dir] for n_dir forward or backward directions dir for vector fields
void PT::shift_field_vec(IFloat **v, const int *dir, int n_dir,
		    int hop, IFloat **u){
  //printf("shift_field_vec() called\n");
  int i, length;
  int wire[n_dir];
  for (i=0; i<n_dir;i++) wire[i] = dir[i];
  SCUDirArgMulti SCUmulti;
  SCUDirArgIR *SCUarg_p[2*n_dir];
  
  int comms=0;
  for (i=0; i<n_dir; i++) 
  if (!local[wire[i]/2]){
    SCUarg_p[2*comms] = SCUarg[hop-1][2*wire[i]];
    SCUarg_p[2*comms+1] = SCUarg[hop-1][2*wire[i]+1];
    SCUarg_p[2*comms+1]->Addr((void *)(v[i]+VECT_LEN*set_offset(wire[i], hop)));
    comms++;
  }

  if (comms) SCUmulti.Init(SCUarg_p,2*comms);
  if (comms) SCUmulti.SlowStartTrans();
//  SCUmulti.TransComplete();
  
  for (i=0; i<n_dir; i++) {
    length = vol-hop*non_local_chi[wire[i]];
    copy_vector(u[i],v[i],&length,dest_l[hop-1][wire[i]],
		src_l[hop-1][wire[i]]);
  }
  if (comms) SCUmulti.TransComplete();

  for (i=0; i<n_dir; i++) {
    length = hop*non_local_chi[wire[i]];
    copy_vector(u[i],(IFloat*)rcv_buf[wire[i]],&length,
		dest_nl[hop-1][wire[i]],src_nl[hop-1][wire[i]]);
  }

}

//! u[-/+nu](x) = U_[-/+nu](x) 
void PT::shift_link(IFloat **u, const int *dir, int n_dir){

//  char *fname = "pt_shift_link()";
  int length;
  for (int i=0; i<n_dir; i++) {
    
    length = local_chi[dir[i]];
    copy_gauge(u[i],uc_l[dir[i]],&length,dest_l[0][dir[i]]);
    length = non_local_chi[dir[i]];
    copy_gauge(u[i],uc_nl[dir[i]],&length,dest_nl[0][dir[i]]);
    
  }

}
