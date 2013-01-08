#ifdef USE_QMP
#include <util/omp_wrapper.h>
/*! \file
  \brief  Definition of parallel transport definitions for QCDOC.
  
  $Id: pt_mat.C,v 1.9 2013-01-08 21:09:25 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2013-01-08 21:09:25 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_base/qmp/pt_mat.C,v 1.9 2013-01-08 21:09:25 chulwoo Exp $
//  $Id: pt_mat.C,v 1.9 2013-01-08 21:09:25 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: pt_mat.C,v $
//  $Revision: 1.9 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_base/qmp/pt_mat.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
#include <string.h>
#include "asq_data_types.h"
#include "pt_int.h"
#include "pt_qcdoc.h"


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

static const int MAX_DIR=10;

#define PROFILE
#undef PROFILE
void PT::mat_cb_norm(int n, IFloat **mout, IFloat **min, const int *dir, int
parity, IFloat * gauge)
{
  //List of the different directions
  int wire[MAX_DIR];
  int i;
//  printf("PT::mat_cb_norm\n");

  QMP_msgmem_t *msg_mem_p = (QMP_msgmem_t *)Alloc("","vec_cb_norm", "msg_mem_p", 2*non_local_dirs*sizeof(QMP_msgmem_t));
  QMP_msghandle_t* msg_handle_p = (QMP_msghandle_t *)Alloc("","vec_cb_norm", "msg_handle_p", 2*non_local_dirs*sizeof(QMP_msghandle_t));
  QMP_msghandle_t multiple;
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
  int non_local_dir=0;

//#pragma omp parallel default(shared)
{

  //If wire[i] is odd, then we have parallel transport in the
  //positive direction.  In this case, multiplication by the link matrix is
  //done before the field is transferred over to the adjacent node
  //
  //If we have transfer in the negative T direction (wire[i] = 6), then
  //we have to copy the appropriate fields to a send buffer
//#pragma omp for
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
	      memcpy(snd_buf_t_cb + j*GAUGE_LEN,min[i] + 3 * *(Toffset[parity]+j)*3,GAUGE_LEN*sizeof(IFloat));
	  }
      }
    }

//#pragma omp barrier
//#pragma omp master 
{
  for(i=0;i<n;i++)
    if(!local[wire[i]/2])
    {
      //Calculate the starting address for the data to be sent
      IFloat *addr = min[i] + GAUGE_LEN * offset_cb[wire[i]];

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

      non_local_dir++;

    }

  if(non_local_dir) {
    multiple = QMP_declare_multiple(msg_handle_p, 2*non_local_dir);
    QMP_start(multiple);
  }
} //#pragma omp master {

  //Do local calculations
//#pragma omp for
  for(i=0;i<n;i++)
    {
      if((wire[i]%2 && conjugated) || ((wire[i]%2 == 0) && (conjugated == 0)))
	pt_cmm_cpp(local_chi_cb[wire[i]],(long)uc_l_cb[parity][wire[i]],(long)min[i],(long)mout[i],(long)gauge);
      else
	pt_cmm_dag_cpp(local_chi_cb[wire[i]],(long)uc_l_cb[parity][wire[i]],(long)min[i],(long)mout[i],(long)gauge);
    }

//#pragma omp barrier
//#pragma omp master 
{
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
} //#pragma omp master {

  //If wire[i] is even, then we have transport in the negative direction
  //In this case, the vector field is multiplied by the SU(3) link matrix
  //after all communication is complete
  IFloat *fp0,*fp1;
//#pragma omp for
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
		  memcpy(fp1,fp0,GAUGE_LEN*sizeof(IFloat));
		}
	    }
	}
    }

} //#pragma omp parallel
#ifdef PROFILE
  dtime +=dclock();
  print_flops("",fname,99*vol*n,dtime);
#endif
//  ParTrans::PTflops +=99*n*vol;
}

//-----------------------------------------------------------------------------

//Parallel transport of a matrix. through one hop.
//The matrix min is parallel transported and the result is placed in mout
#if 1
#undef PROFILE
void PT::mat(int n, matrix **mout, matrix **min, const int *dir){
    
  int wire[MAX_DIR];
  int i;
  QMP_msgmem_t msg_mem_p[2*MAX_DIR];
  QMP_msghandle_t msg_handle_p[2*MAX_DIR];
  QMP_msghandle_t multiple;
  static double setup=0.,qmp=0.,localt=0.,nonlocal=0.;
  static int call_num = 0;

  call_num++;
//  char *fname="pt_mat()";
//  VRB.Func("",fname);
//  if (call_num%100==1) printf("PT:mat()\n");

  
  for(i=0;i<n;i++) wire[i] = dir[i]; 
#ifdef PROFILE
  Float dtime2  = - dclock();
#endif

  double dtime = -dclock();

  int non_local_dir=0;
  
  for(i=0;i<n;i++)
  if (!local[wire[i]/2]) {
    //Calculate the address for transfer in a particular direction
    Float * addr = ((Float *)min[i]+GAUGE_LEN*offset[wire[i]]);
    msg_mem_p[2*non_local_dir] = QMP_declare_msgmem((void *)rcv_buf[wire[i]], 3*non_local_chi[wire[i]]*VECT_LEN*sizeof(IFloat));
    msg_mem_p[2*non_local_dir+1] = QMP_declare_strided_msgmem((void *)addr, (size_t)(3*blklen[wire[i]]), numblk[wire[i]], (ptrdiff_t)(3*stride[wire[i]]+3*blklen[wire[i]]));
    
    msg_handle_p[2*non_local_dir] = QMP_declare_receive_relative(msg_mem_p[2*non_local_dir], wire[i]/2, 1-2*(wire[i]%2), 0);
    msg_handle_p[2*non_local_dir+1] = QMP_declare_send_relative(msg_mem_p[2*non_local_dir+1], wire[i]/2, 2*(wire[i]%2)-1, 0);

    non_local_dir++;
  }
  if (call_num==1 && !QMP_get_node_number())
	printf("non_local_dir=%d\n",non_local_dir);

  if(non_local_dir) {
    multiple = QMP_declare_multiple(msg_handle_p, 2*non_local_dir);
    QMP_start(multiple);
  }

  dtime += dclock();
  setup +=dtime;
  dtime = -dclock();
  int if_print = 0;
//  if ( (call_num%10000==1) && (!QMP_get_node_number()) ) if_print=1;

#define USE_TEST2
#ifdef USE_TEST2
//assume nt > n!
    static char *cname="mat()";
#pragma omp parallel default(shared)
{
  int iam,nt,ipoints,istart,offset;
  iam = omp_get_thread_num();
  nt = omp_get_num_threads();
  int nt_dir = nt/n;
  int n_t = iam/nt_dir;
  int i_t = iam%nt_dir;
  if (n_t >= n ){  n_t = n-1;
    i_t = iam - (n-1)*nt_dir;
    nt_dir = nt -(n-1)*nt_dir;
  }
  int w_t = wire[n_t];
  ipoints = (local_chi[w_t]/2)/nt_dir;
  offset = ipoints*i_t;
  if (i_t == (nt_dir-1)) ipoints = (local_chi[w_t]/2)-offset;
    if ( if_print )
      printf("thread %d of %d nt_dir n_t i_t ipoints offset= %d %d %d %d %d\n",iam,nt,nt_dir,n_t,i_t,ipoints,offset);
  //Interleaving of local computation of matrix multiplication
  partrans_cmm_agg((uc_l[w_t]+offset*2),min[n_t],mout[n_t],ipoints);
    if ( if_print )
      printf("thread %d of %d done\n",iam,nt);
}
#else
{
  //Interleaving of local computation of matrix multiplication
#pragma omp parallel for default(shared)
  for(i=0;i<n;i++){
    partrans_cmm_agg(uc_l[wire[i]],min[i],mout[i],local_chi[wire[i]]/2);
  }
}
#endif

  dtime += dclock();
  localt +=dtime;
  dtime = -dclock();
//#pragma omp barrier
//#pragma omp master 
{
  if(non_local_dir) {
    QMP_status_t qmp_complete_status = QMP_wait(multiple);
    if (qmp_complete_status != QMP_SUCCESS)
          QMP_error("Send failed in vec_cb_norm: %s\n", QMP_error_string(qmp_complete_status));
    QMP_free_msghandle(multiple);
    for(int i = 0; i < 2*non_local_dir; i++)
      QMP_free_msgmem(msg_mem_p[i]);
//    Free(msg_handle_p);
//    Free(msg_mem_p);
  }
} //#pragma omp master {
  dtime += dclock();
  qmp +=dtime;
  dtime = -dclock();

  //Do non-local computations
#ifdef USE_TEST2
//assume nt > n!
#pragma omp parallel default(shared)
{
  int iam,nt,ipoints,istart,offset;
  iam = omp_get_thread_num();
  nt = omp_get_num_threads();
  int nt_dir = nt/n;
  int n_t = iam/nt_dir;
  int i_t = iam%nt_dir;
  if (n_t >= n ){  n_t = n-1;
    i_t = iam - (n-1)*nt_dir;
    nt_dir = nt -(n-1)*nt_dir;
  }
  int w_t = wire[n_t];
  ipoints = (non_local_chi[w_t]/2)/nt_dir;
  offset = ipoints*i_t;
  if (i_t == (nt_dir-1)) ipoints = (non_local_chi[w_t]/2)-offset;
    if ( if_print )
      printf("thread %d of %d nt_dir n_t i_t ipoints offset= %d %d %d %d %d\n",iam,nt,nt_dir,n_t,i_t,ipoints,offset);
  //Non-local computation
  if (ipoints>0)
  partrans_cmm_agg((uc_nl[w_t]+offset*2),(matrix *)rcv_buf[w_t],mout[n_t],ipoints);
    if ( if_print )
      printf("thread %d of %d done\n",iam,nt);
}
#else
{
#pragma omp parallel for
  for(i=0;i<n;i++) 
  if (!local[wire[i]/2]) {
#ifdef USE_OMP
    if (call_num%10000==1 && !QMP_get_node_number() ) 
      printf("thread %d of %d i=%d\n",omp_get_thread_num(),omp_get_num_threads(),i);
#endif
    partrans_cmm_agg(uc_nl[wire[i]],(matrix *)rcv_buf[wire[i]],mout[i],non_local_chi[wire[i]]/2);
  }

}//#pragma omp parallel
#endif

  dtime += dclock();
  nonlocal +=dtime;

  if (call_num%100==0){
    static char *cname="mat()";
    if (!QMP_get_node_number() ) {
    print_flops("mat():local*100",0,localt);
    print_flops("mat():nonlocal*100",0,nonlocal);
    print_flops("mat():qmp*100",0,qmp);
    print_flops("mat():setup*100",0,setup);
    }
    localt=nonlocal=qmp=setup=0.;
  }


#ifdef PROFILE
  dtime2 +=dclock();
  print_flops("",fname,198*vol*n,dtime2);
#endif
//  ParTrans::PTflops +=198*n*vol;
}
#endif

#endif
