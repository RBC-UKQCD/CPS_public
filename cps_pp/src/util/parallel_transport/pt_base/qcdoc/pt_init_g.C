/*! \file
  \brief  Definition of parallel transport definitions for QCDOC.
  
  $Id: pt_init_g.C,v 1.1 2006-07-05 18:13:49 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2006-07-05 18:13:49 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_base/qcdoc/pt_init_g.C,v 1.1 2006-07-05 18:13:49 chulwoo Exp $
//  $Id: pt_init_g.C,v 1.1 2006-07-05 18:13:49 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: pt_init_g.C,v $
//  $Revision: 1.1 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_base/qcdoc/pt_init_g.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
#include "asq_data_types.h"
#include "pt_int.h"

static unsigned long PEC = 0xb0000000;
static unsigned long PLB = 0xb0000000;

#define TESTING
#undef TESTING
#undef CPP
//CPS_START_NAMESPACE
//void dirac_cmv_jcw_agg_cpp( int sites, long chi, long u,long in, long out);

//External function definitions
extern "C"{

  //matrix multiply for checkerboarded fields
  void pt_cmm_cpp(int sites, long u, long in, long out, long gauge_field);
  void pt_cmm_dag_cpp(int sites, long u, long in, long out, long gauge_field);

  //------------------------------------------------------------------------
  //C++ routines
#ifdef CPP
  //Matrix multiplication for full matrix fields
  void cmm_agg_cpp(gauge_agg *chi, matrix *phi, matrix *result, int counter);
  void cmv_agg_cpp( int sites, long u,long in, long out);
  #define partrans_cmm_agg(A,B,C,D) cmm_agg_cpp(A,B,C,D)
  #define partrans_cmv_agg(A,B,C,D) cmv_agg_cpp(A,B,C,D)

  //matrix vector multiply for checkerboarded fields
  void pt_cmv_cpp(int sites, ind_agg *agg, double *gauge_field, double *src, double *dest);
  void pt_cmv_dag_cpp(int sites, ind_agg *agg, double *gauge_field, double *src, double *dest);
  void pt_cmv_pad_cpp(int sites, ind_agg *agg, double *gauge_field, double *src, double *dest);
  void pt_cmv_dag_pad_cpp(int sites, ind_agg *agg, double *gauge_field, double *src, double *dest);
  #define partrans_cmv(A,B,C,D,E) pt_cmv_cpp(A,B,C,D,E)
  #define partrans_cmv_dag(A,B,C,D,E) pt_cmv_dag_cpp(A,B,C,D,E)
  #define partrans_cmv_pad(A,B,C,D,E) pt_cmv_pad_cpp(A,B,C,D,E)
  #define partrans_cmv_dag_pad(A,B,C,D,E) pt_cmv_dag_pad_cpp(A,B,C,D,E)
  //--------------------------------------------------------------------------
  //Assembly Routines
#else
  //Matrix multiplication for full matrix fields
  void pt_cmm_agg(gauge_agg *chi, matrix *phi,matrix *result, int counter);
  //void cmm_agg(gauge_agg *chi, matrix *phi,matrix *result, int counter);
  void pt_asqtad_agg( int sites, long chi, long u,long in, long out);
  void pt_asqtad_agg_s( int sites, long chi, long u,long in, long out);
  #define partrans_cmm_agg(A,B,C,D) pt_cmm_agg(A,B,C,D)
  #define partrans_cmv_agg(A,B,C,D) pt_asqtad_agg(A,0,B,C,D)

  void pt_cmv(int count, ind_agg *ind, double *gauge, double *src, double *dest);
  void pt_cmv_pad(int count, ind_agg *ind, double *gauge, double *src, double *dest);
  void pt_cmv_dag(int count, ind_agg *ind, double *gauge, double *src, double *dest);
  void pt_cmv_dag_pad(int count, ind_agg *ind, double *gauge, double *src, double *dest);
  void pt_cmv_s(int count, ind_agg *ind, float *gauge, float *src, float *dest);
  void pt_cmv_pad_s(int count, ind_agg *ind, float *gauge, float *src, float *dest);
  void pt_cmv_dag_s(int count, ind_agg *ind, float *gauge, float *src, float *dest);
  void pt_cmv_dag_pad_s(int count, ind_agg *ind, float *gauge, float *src, float *dest);
  #define partrans_cmv(A,B,C,D,E) pt_cmv(A,B,C,D,E)
  #define partrans_cmv_dag(A,B,C,D,E) pt_cmv_dag(A,B,C,D,E)
  #define partrans_cmv_pad(A,B,C,D,E) pt_cmv_pad(A,B,C,D,E)
  #define partrans_cmv_dag_pad(A,B,C,D,E) pt_cmv_dag_pad(A,B,C,D,E)
#endif
  //--------------------------------------------------------------------------

  void pt_copy(int count, ind_agg *ind, double *src, double *dest);
  void pt_copy_pad(int count, ind_agg *ind, double *src, double *dest);
  void pt_copy_s(int count, ind_agg *ind, float *src, float *dest);
  void pt_copy_pad_s(int count, ind_agg *ind, float *src, float *dest);

  void pt_copy_buffer(int n, long src, long dest, long ptable);
  // Assembler copying routines
  void copy_matrix(IFloat *res, IFloat *src, int *length, 
		   unsigned long *res_ptr, unsigned long *src_ptr);
  void copy_gauge(IFloat *res, struct gauge_agg *src, int *length,
		  unsigned long *res_ptr);
  // This is perhaps overkill but gives a couple of extra flops
  // cross_look - all input fields are lookup and sum to result
  void cross_look(IFloat *result, Float *fac, const IFloat *chi, const IFloat *phi,  
	     int counter, unsigned long *src, unsigned long *dest);
  // cross_lin - one input field is linear and sum to result
  void cross_lin(IFloat *result, Float *fac, const IFloat *chi, const IFloat *phi,  
	     int counter, unsigned long *dest, unsigned long *dest);
  // cross_over_look - all input fields are lookup and overwrite result
  void cross_over_look(IFloat *result, Float *fac, const IFloat *chi, const IFloat *phi,  
	     int counter, unsigned long *src, unsigned long *dest);
  // cross_over_lin - one input field is linear and overwrite result
  void cross_over_lin(IFloat *result, Float *fac, const IFloat *chi, const IFloat *phi,  
	     int counter, unsigned long *dest, unsigned long *dest);
  //Copies a vectors from v to u
  void copy_vector(IFloat *u, IFloat *v, int *length, unsigned long *dest, unsigned long *src);

  //---------------------------------------------------------------------------
  
  //---------------------------------------------------------------------------

  void m1m2_lookup(matrix *result, matrix *m1, matrix *m2, int length,
		   unsigned long *dest, unsigned long *dest, unsigned long *src);
  void m1m2_lookup_copy(matrix *result2, matrix *result, matrix *m1, matrix *m2, 
			int length, unsigned long *dest2,  
			unsigned long *dest, unsigned long *dest, 
			unsigned long *src);
  void m1m2_lin_copy(matrix *result2, matrix *result, matrix *m1, matrix *m2, 
		     int length, unsigned long *dest2,
		     unsigned long *dest, unsigned long *dest);
  
}
inline  void pt_cmm_agg_print(gauge_agg *chi, matrix *phi,matrix *result, int counter){
   printf("pt_cmm_agg(%p %p %p %d)\n",chi,phi,result,counter);
//    for(int i =0;i<2*counter;i++){
//      printf("%d: %d %d\n",i,chi[i].src,chi[i].dest);
//    }
   printf("pt_cmm_agg(%p %p %p %d) done \n",chi,phi,result,counter);
}

inline  void cross_over_lin_cpp(IFloat *result, Float *fac, const IFloat *chi,
 const IFloat *phi,  int counter, unsigned long *src, unsigned long *dest){
    printf("cross_over_lin(%p %0.4f %p %p %d %p %p)\n",
    result,*fac,chi,phi,counter,src,dest);
    for(int i =0;i<counter;i++){
      printf("%d: %d %d\n",i,src[i],dest[i]);
    }
    cross_over_lin(result,fac,chi,phi,counter,src,dest);
    printf("cross_over_lin(%p %0.4f %p %p %d %p %p) done\n");
}

inline  void cross_over_look_cpp(IFloat *result, Float *fac, const IFloat *chi,
 const IFloat *phi,  int counter, unsigned long *src, unsigned long *dest){
    printf("cross_over_look(%p %0.4f %p %p %d %p %p)\n",
    result,*fac,chi,phi,counter,src,dest);
    for(int i =0;i<counter;i++){
      printf("%d: %d %d\n",i,src[i],dest[i]);
    }
    cross_over_look(result,fac,chi,phi,counter,src,dest);
    printf("cross_over_look(%p %0.4f %p %p %d %p %p) done\n");
}

#ifdef ASQD_SINGLE
#define pt_asqtad_agg(A,B,C,D,E) pt_asqtad_agg_s(A,B,C,D,E)
#define pt_cmv(A,B,C,D,E) pt_cmv_s(A,B,C,D,E)
#define pt_cmv_dag(A,B,C,D,E) pt_cmv_dag_s(A,B,C,D,E)
#define pt_cmv_pad(A,B,C,D,E) pt_cmv_pad_s(A,B,C,D,E)
#define pt_cmv_dag_pad(A,B,C,D,E) pt_cmv_dag_pad_s(A,B,C,D,E)
#define pt_copy_pad(A,B,C,D) pt_copy_pad_s(A,B,C,D)
#define pt_copy(A,B,C,D) pt_copy_s(A,B,C,D)
#endif


//Free memory associated with gauge parallel transport
void PT::delete_g_buf(){
  char *fname = "pt_delete_g()";
//  VRB.Func("",fname);
//  printf("gauge_txyz=%p\n",gauge_txyz);
//  Free(gauge_txyz);
  for(int hop = 0; hop < MAX_HOP; hop++)
    for(int i = 0; i < 4*NDIM; i++)
      if(!local[i/4])
      {
	delete SCUarg[hop][i];
	delete SCUarg2[hop][i];
	delete SCUarg_mat[hop][i];
      }

  //----------------------------------------------------------------------
  //Checkerboarding
  for(int i = 0; i < 4*NDIM; i++)
    if(!local[i/4])
    {
      delete SCUarg_cb[i];
      delete SCUarg_mat_cb[i];
      }
  //---------------------------------------------------------------------
}

void PT::init_g(Float * g_addr){
  int x[NDIM], nei[NDIM];
  int local_count[2*NDIM];
  int non_local_count[2*NDIM];
  int i;

  char *fname = "init_g()";
//printf("%s\n",fname);

  for(i=0; i<2*NDIM;i++){
    local_count[i]=non_local_count[i]=0;
  }
  //Location of gauge field
  IFloat *u = gauge_field_addr;
  if (g_addr) u = g_addr;

  //For staggered parallel transport, we need to re-order the gauge fields
  //to match the ordering of the vector field



  //Send and receive directions
  SCUDir rcv_dir[]={SCU_XP, SCU_XM, SCU_YP, SCU_YM, SCU_ZP, SCU_ZM,SCU_TP,SCU_TM};
  SCUDir snd_dir[]={SCU_XM, SCU_XP, SCU_YM, SCU_YP, SCU_ZM, SCU_ZP,SCU_TM,SCU_TP};

  //Temporary buffer (allocated on cache) that receives an SU(3) matrix
  IFloat *rcv_mat = (IFloat *)Alloc(cname,fname,"rcv_mat",18*sizeof(IFloat),QFAST|QNONCACHE);

  sys_cacheflush(0);

  for(i=0;i<NDIM;i++){
    //Initialize SCUDirArg for sending and receiving the one-hop term
    SCUDirArgIR snd;
    SCUDirArgIR rcv;
    if(!local[i]){
      snd.Init(u,snd_dir[i*2+1],SCU_SEND,sizeof(matrix));
      rcv.Init(rcv_mat,rcv_dir[i*2+1],SCU_REC,sizeof(matrix));
    }

    for(x[3]=0,nei[3]=0;x[3]<size[3];x[3]++,nei[3]++)
      for(x[2]=0,nei[2]=0;x[2]<size[2];x[2]++,nei[2]++)
	for(x[1]=0,nei[1]=0;x[1]<size[1];x[1]++,nei[1]++)
	  for(x[0]=0,nei[0]=0;x[0]<size[0];x[0]++,nei[0]++){
//printf("%d %d %d %d %d\n",i,x[0],x[1],x[2],x[3]);
	    // positive direction
	    //this is a hop in the positive direction, meaning data must be sent
	    //in the negative direction.
	    if((x[i] == 0) && (!local[i]) ){
	      //Calculate the appropriate coordinate on the adjacent node
	      nei[i] = size[i]-1;  
	      //Copy the appropriate matrix from u to uc_nl
	      Copy((uc_nl[2*i]+non_local_count[2*i])->mat, u+LexGauge(nei,i)*GAUGE_LEN);
	      non_local_count[i*2]++;
	      if (non_local_count[i*2]>non_local_chi[i*2])
		fprintf(stderr,"%s:non_local_count[%d](%d)>non_local_chi[%d](%d)\n",
			fname,2*i,non_local_count[2*i],2*i,non_local_chi[2*i]);
	    } else {
	      //Calculate the appropriate neighbor coordinate on the local node
	      nei[i] = (x[i]-1+size[i])%size[i];
	      //Copy from u to uc_l
	      Copy((uc_l[2*i]+local_count[2*i])->mat, u+LexGauge(nei,i)*GAUGE_LEN);
	      local_count[i*2]++;
	      if (local_count[i*2+1]>local_chi[i*2+1])
		fprintf(stderr,"%s:local_count[%d](%d)>local_chi[%d](%d)\n",
			fname,2*i+1,local_count[2*i+1],2*i+1,local_chi[2*i+1]);
	    }
	    // negative direction
	    if((x[i] == (size[i]-1) ) && (!local[i]) ){
	      nei[i] = 0;
	      //Send the appropriate matrix
	      snd.Addr(u+LexGauge(x,i)*GAUGE_LEN);
	      //Send the transmission, prepare to receive
	      snd.StartTrans();rcv.StartTrans();
	      //Complete the send and receive
	      snd.TransComplete();rcv.TransComplete();
	      //Copy to uc_nl from the received matrix
	      DagCopy((uc_nl[2*i+1]+non_local_count[2*i+1])->mat, rcv_mat);
	      non_local_count[i*2+1]++;
	      if (non_local_count[i*2]>non_local_chi[i*2])
		fprintf(stderr,"%s:non_local_count[%d](%d)>non_local_chi[%d](%d)\n",
			fname,2*i,non_local_count[2*i],2*i,non_local_chi[2*i]);
	    } else {
	      //Calculate the appropriate neighbor coordinate on the local volume
	      nei[i] = (x[i]+1)%size[i];
	      //Copy from u to uc_l
	      DagCopy((uc_l[2*i+1]+local_count[2*i+1])->mat, u+LexGauge(x,i)*GAUGE_LEN);
	      local_count[i*2+1]++;
	      if (local_count[i*2+1]>local_chi[i*2+1])
		fprintf(stderr,"%s:local_count[%d](%d)>local_chi[%d](%d)\n",
			fname,2*i+1,local_count[2*i+1],2*i+1,local_chi[2*i+1]);
	    }
	    nei[i] = x[i];
	  } // x[]
  } // i
  //Loop over all possible communication directions
  for(i=0;i<2*NDIM;i++) 
    if (!local[i/2]) {
      for (int hop=1; hop<=MAX_HOP; hop++) {
//printf("%d %d\n",i,hop);
      //Initialize SCUArg to receive fermion fields
      SCUarg[hop-1][i*2] = new SCUDirArgIR;
      SCUarg[hop-1][i*2]->Init((void *)rcv_buf[i],rcv_dir[i],SCU_REC,
			       hop*non_local_chi[i]*VECT_LEN*sizeof(IFloat),1,0,IR_9);
      //Initialize SCUArg to send fermion field
      SCUarg[hop-1][i*2+1] = new SCUDirArgIR;
      SCUarg[hop-1][i*2+1]->Init((void *)rcv_buf[i],snd_dir[i],SCU_SEND,hop*blklen[i],
				 numblk[i],stride[i]+(1-hop)*blklen[i],IR_9);
    
      //Receive for Matrices
      SCUarg_mat[hop-1][i*2] = new SCUDirArgIR;
      SCUarg_mat[hop-1][i*2]->Init((void *)rcv_buf[i],rcv_dir[i],SCU_REC,
				   3*hop*non_local_chi[i]*VECT_LEN*sizeof(IFloat),1,0,IR_9);
      //send for matrices
      SCUarg_mat[hop-1][i*2+1] = new SCUDirArgIR;
      SCUarg_mat[hop-1][i*2+1]->Init((void *)rcv_buf[i],snd_dir[i],SCU_SEND,
				     3*hop*blklen[i],numblk[i],
				     3*(stride[i]+(1-hop)*blklen[i]),IR_9);
      //Receive
      SCUarg2[hop-1][i*2] = new SCUDirArgIR;
      SCUarg2[hop-1][i*2]->Init((void *)rcv_buf2[i],rcv_dir[i],SCU_REC,
				hop*non_local_chi[i]*VECT_LEN*sizeof(IFloat),1,0,IR_9);
      //Send
      SCUarg2[hop-1][i*2+1] = new SCUDirArgIR;
      SCUarg2[hop-1][i*2+1]->Init((void *)rcv_buf2[i],snd_dir[i],SCU_SEND,hop*blklen[i],
				  numblk[i],stride[i]+(1-hop)*blklen[i],IR_9);
      }

      //-----------------------------------------------------------------------
      //Initialize SCUArg to receive fermion fields
      SCUarg_cb[i*2] = new SCUDirArgIR;
//printf("SCUarg_cb\n");
      SCUarg_cb[i*2]->Init((void *)rcv_buf[i],rcv_dir[i],SCU_REC,
			     non_local_chi_cb[i]*VECT_LEN*sizeof(IFloat),1,0,IR_9);

      //Initialize SCUArg to send fermion field
      SCUarg_cb[i*2+1] = new SCUDirArgIR;
//printf("SCUarg_cb\n");
      
      if(i%2)
	SCUarg_cb[i*2+1]->Init((void *)rcv_buf[i],snd_dir[i],SCU_SEND,
			       non_local_chi_cb[i]*VECT_LEN*sizeof(IFloat),1,0,IR_9);
      else 
	SCUarg_cb[i*2+1]->Init((void *)rcv_buf[i],snd_dir[i],SCU_SEND,blklen_cb[i],
				   numblk_cb[i],stride_cb[i],IR_9);
    
      //Receive for Matrices
      SCUarg_mat_cb[i*2] = new SCUDirArgIR;
//printf("SCUarg_mat_cb\n");
	SCUarg_mat_cb[i*2]->Init((void *)rcv_buf[i],rcv_dir[i],SCU_REC,
				   3*non_local_chi_cb[i]*VECT_LEN*sizeof(IFloat),1,0,IR_9);

      //send for matrices
      SCUarg_mat_cb[i*2+1] = new SCUDirArgIR;
//printf("SCUarg_mat_cb\n");
      if(i%2)
	SCUarg_mat_cb[i*2+1]->Init((void *)rcv_buf[i],snd_dir[i],SCU_SEND,
				   3*non_local_chi_cb[i]*VECT_LEN*sizeof(IFloat),1,0,IR_9);
      else
	SCUarg_mat_cb[i*2+1]->Init((void *)rcv_buf[i],snd_dir[i],SCU_SEND,
				   3*blklen_cb[i],numblk_cb[i],3*stride_cb[i],IR_9);
      //-----------------------------------------------------------------------
    }
  

  Free(rcv_mat);
//printf("%s done\n",fname);
}

//CPS_END_NAMESPACE
