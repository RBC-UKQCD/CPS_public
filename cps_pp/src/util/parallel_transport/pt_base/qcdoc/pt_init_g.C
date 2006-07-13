/*! \file
  \brief  Definition of parallel transport definitions for QCDOC.
  
  $Id: pt_init_g.C,v 1.2 2006-07-13 20:14:44 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2006-07-13 20:14:44 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_base/qcdoc/pt_init_g.C,v 1.2 2006-07-13 20:14:44 chulwoo Exp $
//  $Id: pt_init_g.C,v 1.2 2006-07-13 20:14:44 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $RCSfile: pt_init_g.C,v $
//  $Revision: 1.2 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/parallel_transport/pt_base/qcdoc/pt_init_g.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------
#include "asq_data_types.h"
#include "pt_int.h"

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
