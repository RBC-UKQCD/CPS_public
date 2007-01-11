/*! \file
  \brief  Definition of parallel transport definitions for QCDOC.
  
  $Id: pt_mat.C,v 1.2 2007-01-11 22:45:57 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2007-01-11 22:45:57 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_base/qmp/pt_mat.C,v 1.2 2007-01-11 22:45:57 chulwoo Exp $
//  $Id: pt_mat.C,v 1.2 2007-01-11 22:45:57 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: pt_mat.C,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_base/qmp/pt_mat.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
#include <string.h>
#include "asq_data_types.h"
#include "pt_int.h"
#include "pt_qcdoc.h"

#ifndef USE_QMP
#define USE_QMP
#endif

//Parallel transport of a matrix defined on one half of the
//checkerboaded lattice
//
//Parameters
//
//n - The number of direction in which to perform the parallel transport
//mout - Result of the parallel transport, on sites with opposite parity of min
//min - Initial field, defined on sites with only one parity
//dir - a list of the n directions in which the field will be transported
//cb - Checkerboard parity of the vector min

#undef PROFILE
void PT::mat_cb(int n, IFloat **mout, IFloat **min, const int *dir, int
parity, IFloat * new_gauge_field)
{
  mat_cb_norm(n,mout,min,dir,parity,new_gauge_field);
}

#undef PROFILE
void PT::mat_cb(int n, IFloat **mout, IFloat **min, const int *dir, int
parity)
{
  mat_cb_norm(n,mout,min,dir,parity,gauge_field_addr);
}

#define PROFILE
#undef PROFILE
void PT::mat_cb_norm(int n, IFloat **mout, IFloat **min, const int *dir, int
parity, IFloat * gauge)
{
  //List of the different directions
  int wire[n];
  int i;

  #ifdef USE_QMP
  QMP_msgmem_t *msg_mem_p = (QMP_msgmem_t *)Alloc("","vec_cb_norm", "msg_mem_p", 2*non_local_dirs*sizeof(QMP_msgmem_t));
  QMP_msghandle_t* msg_handle_p = (QMP_msghandle_t *)Alloc("","vec_cb_norm", "msg_handle_p", 2*non_local_dirs*sizeof(QMP_msghandle_t));
  QMP_msghandle_t multiple;
  #else
  //SCUDirArgs for sending and receiving in the n directions
  SCUDirArgIR *SCUarg_p[2*n];
  SCUDirArgMulti SCUmulti;
  #endif
  static int call_num = 0;
  int vlen = VECT_LEN;
  int vlen2 = VECT_LEN;

  call_num++;
  
  //Name our function
  char *fname="pt_mat_cb()";
  //  VRB.Func("",fname);
  
  //Set the transfer directions
  //If wire[i] is even, then we have communication in the negative direction
  //If wire[i] is odd, then we have communication in the positive direction
  for(i=0;i<n;i++)
    wire[i]=dir[i];

#ifdef PROFILE
  Float dtime  = - dclock();
#endif

  //If wire[i] is odd, then we have parallel transport in the
  //positive direction.  In this case, multiplication by the link matrix is
  //done before the field is transferred over to the adjacent node
  //
  //If we have transfer in the negative T direction (wire[i] = 6), then
  //we have to copy the appropriate fields to a send buffer
  for(i=0;i<n;i++)
    {
      if(!local[wire[i]/2])
      {
	if(wire[i]%2)
	  {
	    if(conjugated)
	      pt_cmm_cpp(non_local_chi_cb[wire[i]],(long)uc_nl_cb_pre[parity][wire[i]/2],(long)min[i],(long)snd_buf_cb[wire[i]/2],(long)gauge);
	    else
	      pt_cmm_dag_cpp(non_local_chi_cb[wire[i]],(long)uc_nl_cb_pre[parity][wire[i]/2],(long)min[i],(long)snd_buf_cb[wire[i]/2],(long)gauge);
	  }
	else if((wire[i] == 6))
	  {
	    for(int j = 0; j < non_local_chi_cb[6];j++)
//	      moveMem(snd_buf_t_cb + j*GAUGE_LEN,min[i] + 3 * *(Toffset[parity]+j)*3,GAUGE_LEN*sizeof(IFloat));
	      memcpy(snd_buf_t_cb + j*GAUGE_LEN,min[i] + 3 * *(Toffset[parity]+j)*3,GAUGE_LEN*sizeof(IFloat));
	  }
      }
    }

  int non_local_dir=0;
  for(i=0;i<n;i++)
    if(!local[wire[i]/2])
    {
      //Calculate the starting address for the data to be sent
      IFloat *addr = min[i] + GAUGE_LEN * offset_cb[wire[i]];
      #ifndef USE_QMP
      //This points to the appropriate SCUDirArg for receiving
      SCUarg_p[2*non_local_dir] = SCUarg_mat_cb[2*wire[i]];
      //This points to the appropriate SCUDirArg for sending
      SCUarg_p[2*i+non_local_dir] = SCUarg_mat_cb[2*wire[i]+1];
      
      //Set the send address
      if(wire[i]%2)
	SCUarg_p[2*non_local_dir+1]->Addr((void *)snd_buf_cb[wire[i]/2]);
      else if(wire[i] == 6)
	SCUarg_p[2*non_local_dir+1]->Addr((void *)snd_buf_t_cb);
      else
	SCUarg_p[2*non_local_dir+1]->Addr((void *)addr);
      #else

      msg_mem_p[2*non_local_dir] = QMP_declare_msgmem((void *)rcv_buf[wire[i]], 3*non_local_chi_cb[wire[i]]*VECT_LEN*sizeof(IFloat));
      
      //Initialize the msg_mem for sends
      if(wire[i]%2) 
	msg_mem_p[2*non_local_dir+1] = QMP_declare_msgmem((void *)snd_buf_cb[wire[i]/2], 3*non_local_chi_cb[wire[i]]*VECT_LEN*sizeof(IFloat));
      else if(wire[i] == 6)
	msg_mem_p[2*non_local_dir+1] = QMP_declare_msgmem((void *)snd_buf_t_cb, 3*non_local_chi_cb[wire[i]]*VECT_LEN*sizeof(IFloat));
      else
	msg_mem_p[2*non_local_dir+1] = QMP_declare_strided_msgmem((void *)addr, (size_t)(3*blklen_cb[wire[i]]), numblk_cb[wire[i]], (ptrdiff_t)(3*stride_cb[wire[i]]+3*blklen_cb[wire[i]]));
      
      msg_handle_p[2*non_local_dir] = QMP_declare_receive_relative(msg_mem_p[2*non_local_dir], wire[i]/2, 1-2*(wire[i]%2), 0);
      msg_handle_p[2*non_local_dir+1] = QMP_declare_send_relative(msg_mem_p[2*non_local_dir+1], wire[i]/2, 2*(wire[i]%2)-1, 0);
#endif

      non_local_dir++;

    }

  #ifndef USE_QMP
  if(non_local_dir){
    SCUmulti.Init(SCUarg_p,2*non_local_dir);

  //Begin transmission
    SCUmulti.SlowStartTrans();
  }
  #else
  if(non_local_dir) {
    multiple = QMP_declare_multiple(msg_handle_p, 2*non_local_dir);
    QMP_start(multiple);
  }
  #endif

  //Do local calculations
  for(i=0;i<n;i++)
    {
      if((wire[i]%2 && conjugated) || ((wire[i]%2 == 0) && (conjugated == 0)))
	pt_cmm_cpp(local_chi_cb[wire[i]],(long)uc_l_cb[parity][wire[i]],(long)min[i],(long)mout[i],(long)gauge);
      else
	pt_cmm_dag_cpp(local_chi_cb[wire[i]],(long)uc_l_cb[parity][wire[i]],(long)min[i],(long)mout[i],(long)gauge);
    }

  #ifndef USE_QMP
  if(non_local_dir) {
    //End transmission
    SCUmulti.TransComplete();
  }
 #else
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
  #endif  

  //If wire[i] is even, then we have transport in the negative direction
  //In this case, the vector field is multiplied by the SU(3) link matrix
  //after all communication is complete
  IFloat *fp0,*fp1;
  for(i=0;i<n;i++)
    {
      if(!local[wire[i]/2])
      	{
	  if(!(wire[i]%2))
	    {
	      if(conjugated)
		pt_cmm_dag_cpp(non_local_chi_cb[wire[i]],(long)uc_nl_cb[parity][wire[i]],(long)rcv_buf[wire[i]],(long)mout[i],(long)gauge);
	      else
		pt_cmm_cpp(non_local_chi_cb[wire[i]],(long)uc_nl_cb[parity][wire[i]],(long)rcv_buf[wire[i]],(long)mout[i],(long)gauge);
	    }
	  //Otherwise we have parallel transport in the positive direction.
	  //In this case, the received data has already been pre-multiplied
	  //All we need to do is to put the transported field in the correct place
	  else
	    {
	      //int destination, source;
	      //Place the data in the receive buffer into the result vector
	      for(int s=0;s<non_local_chi_cb[wire[i]];s++)
		{
		  //source = uc_nl_cb[parity][wire[i]][s].src;
		  fp0 = (IFloat *)((long)rcv_buf[wire[i]]+3*uc_nl_cb[parity][wire[i]][s].src);
		  //destination = uc_nl_cb[parity][wire[i]][s].dest;
		  fp1 = (IFloat *)(mout[i]+3*uc_nl_cb[parity][wire[i]][s].dest);
		  //for(int d = 0;d<GAUGE_LEN;d++)
		  //*(fp1+d) = *(fp0+d);
//		  moveMem(fp1,fp0,GAUGE_LEN*sizeof(IFloat));
		  memcpy(fp1,fp0,GAUGE_LEN*sizeof(IFloat));
		}
	    }
	}
    }
#ifdef PROFILE
  dtime +=dclock();
  print_flops("",fname,99*vol*n,dtime);
#endif
//  ParTrans::PTflops +=99*n*vol;
}

//-----------------------------------------------------------------------------

//Parallel transport of a matrix. through one hop.
//The matrix min is parallel transported and the result is placed in mout
#undef PROFILE
void PT::mat(int n, matrix **mout, matrix **min, const int *dir){
    
  int wire[n];
  int i;
  #ifdef USE_QMP
  QMP_msgmem_t *msg_mem_p = (QMP_msgmem_t *)Alloc("","vec_cb_norm", "msg_mem_p", 2*non_local_dirs*sizeof(QMP_msgmem_t));
  QMP_msghandle_t* msg_handle_p = (QMP_msghandle_t *)Alloc("","vec_cb_norm", "msg_handle_p", 2*non_local_dirs*sizeof(QMP_msghandle_t));
  QMP_msghandle_t multiple;
  #else
  SCUDirArgIR *SCUarg_p[2*n];
  SCUDirArgMulti SCUmulti;
  #endif
  static int call_num = 0;

  call_num++;
//  char *fname="pt_mat()";
//  VRB.Func("",fname);
  
  for(i=0;i<n;i++) wire[i] = dir[i]; 
#ifdef PROFILE
  Float dtime  = - dclock();
#endif
  int non_local_dir=0;
  for(i=0;i<n;i++)
  if (!local[wire[i]/2]) {
    //Calculate the address for transfer in a particular direction
    Float * addr = ((Float *)min[i]+GAUGE_LEN*offset[wire[i]]);
    #ifndef USE_QMP
    //This should point to the appropriate SCUDirArg for receiving
    SCUarg_p[2*non_local_dir] = SCUarg_mat[0][2*wire[i]];
    //This points to the appropriate SCUDirArg for sending
    SCUarg_p[2*non_local_dir+1] = SCUarg_mat[0][2*wire[i]+1];
    //Reset the send address
    SCUarg_p[2*non_local_dir+1]->Addr((void *)addr);
   #else
    msg_mem_p[2*non_local_dir] = QMP_declare_msgmem((void *)rcv_buf[wire[i]], 3*non_local_chi[wire[i]]*VECT_LEN*sizeof(IFloat));
    msg_mem_p[2*non_local_dir+1] = QMP_declare_strided_msgmem((void *)addr, (size_t)(3*blklen[wire[i]]), numblk[wire[i]], (ptrdiff_t)(3*stride[wire[i]]+3*blklen[wire[i]]));
    
    msg_handle_p[2*non_local_dir] = QMP_declare_receive_relative(msg_mem_p[2*non_local_dir], wire[i]/2, 1-2*(wire[i]%2), 0);
    msg_handle_p[2*non_local_dir+1] = QMP_declare_send_relative(msg_mem_p[2*non_local_dir+1], wire[i]/2, 2*(wire[i]%2)-1, 0);
    #endif

    non_local_dir++;
  }

  #ifndef USE_QMP
if (non_local_dir){
  SCUmulti.Init(SCUarg_p,non_local_dir*2);
  //Start transmission
  SCUmulti.SlowStartTrans();
}
  #else
  if(non_local_dir) {
    multiple = QMP_declare_multiple(msg_handle_p, 2*non_local_dir);
    QMP_start(multiple);
  }
  #endif
  //Interleaving of local computation of matrix multiplication
  for(i=0;i<n;i++){
    partrans_cmm_agg(uc_l[wire[i]],min[i],mout[i],local_chi[wire[i]]/2);
  }

  #ifndef USE_QMP
  if (non_local_dir)
  SCUmulti.TransComplete();
  #else
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
  #endif

  //Do non-local computations
  for(i=0;i<n;i++) 
  if (!local[wire[i]/2]) {
    partrans_cmm_agg(uc_nl[wire[i]],(matrix *)rcv_buf[wire[i]],mout[i],non_local_chi[wire[i]]/2);
  }
#ifdef PROFILE
  dtime +=dclock();
  print_flops("",fname,198*vol*n,dtime);
#endif
//  ParTrans::PTflops +=198*n*vol;
}

