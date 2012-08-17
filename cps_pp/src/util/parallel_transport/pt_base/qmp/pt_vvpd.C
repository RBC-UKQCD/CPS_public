#ifdef USE_QMP
/*! \file
  \brief  Definition of parallel transport definitions for QCDOC.
  
  $Id: pt_vvpd.C,v 1.5 2012-08-17 05:53:54 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2012-08-17 05:53:54 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_base/qmp/pt_vvpd.C,v 1.5 2012-08-17 05:53:54 chulwoo Exp $
//  $Id: pt_vvpd.C,v 1.5 2012-08-17 05:53:54 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: pt_vvpd.C,v $
//  $Revision: 1.5 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_base/qmp/pt_vvpd.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
#include "asq_data_types.h"
#include "pt_int.h"
#include "pt_qcdoc.h"
#ifndef USE_QMP
#define USE_QMP
#endif

static int MAX_DIR=10;

/*! 
  Computes sum[x] = vect[x] vect[x + hop dir]^dagger
  where the sum is over n_vect vectors and the hop is in a forward direction.
*/
void PT::vvpd(IFloat **vect, int n_vect, const int *dir, 
	     int n_dir, int hop, IFloat **sum){

  //---------------------------------------------------------------------
  //Adapt old vvpd for usage with new interface

  //Convert direction array
  int newdir[MAX_DIR];
  for(int i = 0; i < n_dir; i++)
    newdir[i] = 2*dir[i];
  IFloat ***vect_shift = (IFloat ***)FastAlloc("","","vect_shift",n_vect*sizeof(IFloat**));
  for(int i = 0; i < n_vect; i++)
    {
      vect_shift[i] = (IFloat**)FastAlloc("","","vect_shift[i]",n_dir*sizeof(IFloat *));
      for(int j = 0; j < n_dir; j++)
	vect_shift[i][j] = vect[i];
    }
  vvpd(vect, vect_shift, n_vect, newdir, n_dir, hop, sum, 1);
  for(int i = 0; i < n_vect; i++)
    Free(vect_shift[i]);
  Free(vect_shift);
  //--------------------------------------------------------------------
}

/*! 
  Computes sum[x] = vect2[x] vect[x + hop dir]^dagger
  where the sum is over n_vect vectors and the hop is in a forward direction.
*/
void PT::vvpd(IFloat **vect2, IFloat ***vect, int n_vect, const int *dir, int n_dir, int hop, IFloat **sum, int overwrite){

  char *fname = "pt_vvpd()";
#if 1
//  ERR.NotImplemented(cname,fname);
   QMP_error("%s""%s Not implemented\n");
#else
//  VRB.Func("",fname);
  int i, s, v;
  Float f = 2.0;
  int wire[MAX_DIR];
  for(i=0;i<n_dir;i++) wire[i] = dir[i]; // from (x,y,z,t) to (t,x,y,z)

  QMP_msgmem_t *msg_mem_p = (QMP_msgmem_t *)Alloc("","vvpd", "msg_mem_p", 2*non_local_dirs*sizeof(QMP_msgmem_t));
  QMP_msgmem_t *msg_mem_p2 = (QMP_msgmem_t *)Alloc("","vvpd", "msg_mem_p", 2*non_local_dirs*sizeof(QMP_msgmem_t));
  QMP_msghandle_t* msg_handle_p = (QMP_msghandle_t *)Alloc("","vvpd", "msg_handle_p", 2*non_local_dirs*sizeof(QMP_msghandle_t));
  QMP_msghandle_t* msg_handle_p2 = (QMP_msghandle_t *)Alloc("","vvpd", "msg_handle_p", 2*non_local_dirs*sizeof(QMP_msghandle_t));
  QMP_msghandle_t multiple;

  //Setup communciation
  int comms=0;
  for(i=0;i<n_dir;i++)
  if( !local[wire[i]/2]) {
    if ( size[wire[i]/2] <hop)
      fprintf(stderr, 
		"%s:size(%d) in direction %d is smaller than the hop(%d)\n",
		fname,size[wire[i]],wire[i],hop);

    comms++;
  }

  for(v=0; v<n_vect; v++){


    if (v%2==0) {
      comms=0;
      for(i=0;i<n_dir;i++)
        if( !local[wire[i]/2]){ 
	  msg_mem_p[2*comms] = QMP_declare_msgmem((void *)rcv_buf[wire[i]], hop*non_local_chi[wire[i]]*VECT_LEN*sizeof(IFloat));    
	  msg_handle_p[2*comms] = QMP_declare_receive_relative(msg_mem_p[2*comms], wire[i]/2, 1-2*(wire[i]%2), 0);
	  msg_mem_p[2*comms+1] = QMP_declare_strided_msgmem((void *)(vect[v][i]+VECT_LEN*set_offset(wire[i], hop)), (size_t)(hop*blklen[wire[i]]), numblk[wire[i]], (ptrdiff_t)(stride[wire[i]] + blklen[wire[i]]));
	  msg_handle_p[2*comms+1] =  QMP_declare_send_relative(msg_mem_p[2*comms+1], wire[i]/2, 2*(wire[i]%2)-1, 0);
          comms++;
        }

      // Start communication
      if(comms) {
	multiple = QMP_declare_multiple(msg_handle_p, 2*comms);
      }
      if (comms) {
	QMP_start(multiple);
	QMP_status_t qmp_complete_status = QMP_wait(multiple);
	if (qmp_complete_status != QMP_SUCCESS)
	  QMP_error("Send failed in vvpd: %s\n", QMP_error_string(qmp_complete_status));
	QMP_free_msghandle(multiple);
	for(int i = 0; i < 2*comms; i++) 
	  QMP_free_msgmem(msg_mem_p[i]);
      }

    } else {
      comms=0;
      for(i=0;i<n_dir;i++)
        if( !local[wire[i]/2]){ 
	  msg_mem_p2[2*comms] = QMP_declare_msgmem((void *)rcv_buf2[wire[i]], hop*non_local_chi[wire[i]]*VECT_LEN*sizeof(IFloat));  
	  msg_handle_p2[2*comms] = QMP_declare_receive_relative(msg_mem_p2[2*comms], wire[i]/2, 1-2*(wire[i]%2), 0);  
	  msg_mem_p2[2*comms+1] = QMP_declare_strided_msgmem((void *)(vect[v][i]+VECT_LEN*set_offset(wire[i], hop)), (size_t)(hop*blklen[wire[i]]), numblk[wire[i]], (ptrdiff_t)(stride[wire[i]] + blklen[wire[i]]));
	  msg_handle_p2[2*comms+1] =  QMP_declare_send_relative(msg_mem_p2[2*comms+1], wire[i]/2, 2*(wire[i]%2)-1, 0);
          comms++;
	}

      // Start communication
      if(comms) {
	multiple = QMP_declare_multiple(msg_handle_p2, 2*comms);
      }
      if (comms) {
	QMP_start(multiple);
	QMP_status_t qmp_complete_status = QMP_wait(multiple);
	if (qmp_complete_status != QMP_SUCCESS)
	  QMP_error("Send failed in vvpd: %s\n", QMP_error_string(qmp_complete_status));
      QMP_free_msghandle(multiple);
      for(int i = 0; i < 2*comms; i++) 
	QMP_free_msgmem(msg_mem_p2[i]);
      }
    }

     
    

    // Perform non-local calculation for previous v
    if (v>0)
      if (v==1 && overwrite==1) {
	for(i=0; i<n_dir; i++)
	  if(non_local_chi[wire[i]]>0)
	    cross_over_lin(sum[i], &f, vect2[v-1],rcv_buf[wire[i]], hop*non_local_chi[wire[i]],
			   src_nl[hop-1][wire[i]], dest_nl[hop-1][wire[i]]);
      } else if (v%2==1) {
	for(i=0; i<n_dir; i++) 
	  if(non_local_chi[wire[i]]>0)
	    cross_lin(sum[i], &f, vect2[v-1],rcv_buf[wire[i]], hop*non_local_chi[wire[i]],
		      src_nl[hop-1][wire[i]], dest_nl[hop-1][wire[i]]);
      } else {
	for(i=0; i<n_dir; i++) 
	  if(non_local_chi[wire[i]]>0)
	    cross_lin(sum[i], &f,vect2[v-1],rcv_buf2[wire[i]], hop*non_local_chi[wire[i]],
		      src_nl[hop-1][wire[i]], dest_nl[hop-1][wire[i]]);
      }
    
    // Perform local calculation for current v
    if (v==0 && overwrite==1)
      {
	for(i=0; i<n_dir; i++)
	  if((vol-hop*non_local_chi[wire[i]])>0)
	    cross_over_look(sum[i], &f, vect2[v], vect[v][i], vol-hop*non_local_chi[wire[i]], src_l[hop-1][wire[i]], dest_l[hop-1][wire[i]]);
      }
    else
      {
	for(i=0; i<n_dir; i++)
	  if((vol-hop*non_local_chi[wire[i]])>0)
	    cross_look(sum[i], &f, vect2[v], vect[v][i], vol-hop*non_local_chi[wire[i]], src_l[hop-1][wire[i]], dest_l[hop-1][wire[i]]);
      }

  }


  if (v==1 && overwrite==1) {
    for(i=0; i<n_dir; i++) 
      if(non_local_chi[wire[i]]>0)
	cross_over_lin(sum[i], &f, vect2[v-1],rcv_buf[wire[i]], hop*non_local_chi[wire[i]],
		       src_nl[hop-1][wire[i]], dest_nl[hop-1][wire[i]]);
  } else if (v%2==1) {
    for(i=0; i<n_dir; i++)
      if(non_local_chi[wire[i]]>0)
	cross_lin(sum[i], &f, vect2[v-1],rcv_buf[wire[i]], hop*non_local_chi[wire[i]],
		  src_nl[hop-1][wire[i]], dest_nl[hop-1][wire[i]]);
  } else {
    for(i=0; i<n_dir; i++) 
      if(non_local_chi[wire[i]]>0)
	cross_lin(sum[i], &f,vect2[v-1],rcv_buf2[wire[i]], hop*non_local_chi[wire[i]],
		  src_nl[hop-1][wire[i]], dest_nl[hop-1][wire[i]]);
  }
#endif
  
  //  ParTrans::PTflops += 90*n_vect*n_dir*vol;
}
#endif
