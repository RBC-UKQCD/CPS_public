#ifdef USE_QMP
/*! \file
  \brief  Definition of parallel transport definitions for QCDOC.
  
*/
#include "asq_data_types.h"
#include "pt_int.h"
//#include <qmp.h>

#include <util/gjp.h>

//Flag to determine whether to use QMP calls
#ifndef USE_QMP
#define USE_QMP
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

	//Free msg_mem created during initialization (receive only)
	if (i%2 == 0) {
	  QMP_free_msgmem(*(msg_mem[hop][i]));
	  QMP_free_msgmem(*(msg_mem2[hop][i]));
	  QMP_free_msgmem(*(msg_mem_mat[hop][i]));
	}
	delete msg_mem[hop][i];
	delete msg_mem2[hop][i];
	delete msg_mem_mat[hop][i];
      }

  //----------------------------------------------------------------------
  //Checkerboarding
  for(int i = 0; i < 2*NDIM; i++)
    if(!local[i/2])
    {
      
      //Free msg_mem created during initialization (receive only)
      QMP_free_msgmem(msg_mem_cb[i*2]);
      QMP_free_msghandle(msg_handle_cb[i*2]);
      QMP_free_msgmem(msg_mem_mat_cb[i*2]);
      if(i%2){
        QMP_free_msgmem(msg_mem_cb[i*2+1]);
        QMP_free_msghandle(msg_handle_cb[i*2+1]);
      } else if(i==6){
        QMP_free_msgmem(msg_mem_cb[i*2+1]);
        QMP_free_msghandle(msg_handle_cb[i*2+1]);
      }

    }
  //---------------------------------------------------------------------
}

void PT::init_g(Float * g_addr){
  int x[NDIM], nei[NDIM];
  int local_count[2*NDIM];
  int non_local_count[2*NDIM];
  int local_gpbound_count[2*NDIM];
  int i;

  char *fname = "init_g()";
//printf("%s\n",fname);

  for(i=0; i<2*NDIM;i++){
    local_count[i]=non_local_count[i]=local_gpbound_count[i]=0;
  }
  //Location of gauge field
  IFloat *u = gauge_field_addr;
  if (g_addr) u = g_addr;

  //For staggered parallel transport, we need to re-order the gauge fields
  //to match the ordering of the vector field



  //Temporary buffer (allocated on cache) that receives an SU(3) matrix
#ifdef USE_QALLOC
  IFloat *rcv_mat = (IFloat *)Alloc(cname,fname,"rcv_mat",18*sizeof(IFloat),QFAST|QNONCACHE);
#else
  IFloat *rcv_mat = (IFloat *)Alloc(cname,fname,"rcv_mat",18*sizeof(IFloat),0);
#endif


#if TARGET == QCDOC
  sys_cacheflush(0);
#endif

  for(i=0;i<NDIM;i++){


// or 'backwards' psi(x) = U_mu(x)psi(x+mu), requires backwards comms of psi only
//    'forwards'  psi(x) = U_mu^dag(x-mu)psi(x-mu), requires forwards comms of psi and link

    for(x[3]=0,nei[3]=0;x[3]<size[3];x[3]++,nei[3]++)
      for(x[2]=0,nei[2]=0;x[2]<size[2];x[2]++,nei[2]++)
	for(x[1]=0,nei[1]=0;x[1]<size[1];x[1]++,nei[1]++)
	  for(x[0]=0,nei[0]=0;x[0]<size[0];x[0]++,nei[0]++){

//printf("%d %d %d %d %d\n",i,x[0],x[1],x[2],x[3]);
	    // positive direction
	    //this is a hop in the positive direction, meaning data must be sent
	    //in the negative direction.

	    /*
	      CK:
	      Here is the code (from pt_init.C) that sets up the src and dest positions
	      if( (x[i] == 0) && (!local[i])){
	      nei[i] = size[i]-1;
	      (uc_nl[2*i]+non_local_count[2*i])->src = non_local_count[2*i]*vlen;
	      (uc_nl[2*i]+non_local_count[2*i])->dest = LexVector(nei)*vlen2;
	      non_local_count[i*2]++;
	      if (non_local_count[i*2]>non_local_chi[i*2])
	      fprintf(stderr,"%s:non_local_count[%d](%d)>non_local_chi[%d](%d)\n",
	      fname,2*i,non_local_count[2*i],2*i,non_local_chi[2*i]);
	      } 
	      so 'nei' is the site we are parallel transporting *to*
	         'x' is the site we are transporting *from*?
	      This is a backwards transport:
	      psi(nei) = U_mu(nei)psi(nei+mu = x)
	      
	      note link site should be 'nei'
	      requires negative comms of psi only
	    */


	    if((x[i] == 0) && (!local[i]) ){
	      //'wire = 2*i' is even, meaning comms (psi only here) is in the backwards direction, receive from node in front

	      //Calculate the appropriate coordinate on the adjacent node
	      nei[i] = size[i]-1; //is the site one before x. 
	      //Copy the appropriate matrix from u to uc_nl
	      Copy((uc_nl[2*i]+non_local_count[2*i])->mat, u+LexGauge(nei,i)*GAUGE_LEN); //copy link from nei to storage
	      non_local_count[i*2]++;
	      if (non_local_count[i*2]>non_local_chi[i*2])
		fprintf(stderr,"%s:non_local_count[%d](%d)>non_local_chi[%d](%d)\n",
			fname,2*i,non_local_count[2*i],2*i,non_local_chi[2*i]);
	    } 
	    else if(x[i] == 0 && cps::GJP.Bc(i)==cps::BND_CND_GPARITY && local[i]){
	      nei[i] = (x[i]-1+size[i])%size[i];
	      printf("For PT from x=(%d,%d,%d,%d) to y=(%d,%d,%d,%d), U(y) stored\n",x[0],x[1],x[2],x[3],nei[0],nei[1],nei[2],nei[3]);

	      Copy((uc_l_gpbound[2*i]+local_gpbound_count[2*i])->mat, u+LexGauge(nei,i)*GAUGE_LEN); //doesnt need complex conjugating as link is on dest site (L-1)
	      //for(int j=0;j<18;j++) (uc_l_gpbound[2*i]+local_gpbound_count[2*i])->mat[j]*=-1;

	      local_gpbound_count[i*2]++;
	    }
	    else {
	      //Calculate the appropriate neighbor coordinate on the local node
	      nei[i] = (x[i]-1+size[i])%size[i]; //gauge link from site 
	      //Copy from u to uc_l
	      Copy((uc_l[2*i]+local_count[2*i])->mat, u+LexGauge(nei,i)*GAUGE_LEN);
	      
	      local_count[i*2]++;
	      if (local_count[i*2+1]>local_chi[i*2+1])
		fprintf(stderr,"%s:local_count[%d](%d)>local_chi[%d](%d)\n",
			fname,2*i+1,local_count[2*i+1],2*i+1,local_chi[2*i+1]);
	    }
	    // negative direction (CK: what?? cf. below)

	    /*
	      CK:
	      Here is the code that sets up  the src and dest positions
	      if((x[i] == (size[i] -1))  && (!local[i])){
	      nei[i] = 0;
	      (uc_nl[2*i+1]+non_local_count[2*i+1])->src = non_local_count[2*i+1]*vlen;
	      (uc_nl[2*i+1]+non_local_count[2*i+1])->dest = LexVector(nei)*vlen2;
	      non_local_count[i*2+1]++;
	      if (non_local_count[i*2+1]>non_local_chi[i*2+1])
	      fprintf(stderr,"%s:non_local_count[%d](%d)>non_local_chi[%d](%d)\n",
	      fname,2*i+1,non_local_count[2*i+1],2*i+1,non_local_chi[2*i+1]);
	      }
	      forwards transport from x to nei
	      psi(nei) = U^dag_mu(nei-mu = x)psi(nei-mu = x)
	      
	      note link site should be 'x'
	      requires forwards comms of psi and U^dag
	    */

	    if((x[i] == (size[i]-1) ) && (!local[i]) ){
	      nei[i] = 0;

	      #ifdef USE_QMP
	      {
		QMP_msgmem_t snd_msgmem = QMP_declare_msgmem((void *)(u+LexGauge(x,i)*GAUGE_LEN), sizeof(PTmatrix));
		QMP_msgmem_t rcv_msgmem = QMP_declare_msgmem((void *)rcv_mat, sizeof(PTmatrix));
		QMP_msghandle_t snd_msghandle = QMP_declare_send_relative(snd_msgmem, i, +1, 0); //CK: data is sent in *positive* direction
		QMP_msghandle_t rcv_msghandle = QMP_declare_receive_relative(rcv_msgmem, i, -1, 0); //CK: data is received from *negative* direction
		QMP_start(snd_msghandle);
		QMP_start(rcv_msghandle);
		QMP_status_t snd_status = QMP_wait(snd_msghandle);
		QMP_status_t rcv_status = QMP_wait(rcv_msghandle);
		if (snd_status != QMP_SUCCESS) {
		  QMP_error("Send failed: %s\n", QMP_error_string(snd_status));
		}
		if (rcv_status != QMP_SUCCESS) {
		  QMP_error("Send failed: %s\n", QMP_error_string(rcv_status));
		}
		QMP_free_msghandle(snd_msghandle);
		QMP_free_msghandle(rcv_msghandle);
		QMP_free_msgmem(snd_msgmem);
		QMP_free_msgmem(rcv_msgmem);				
	      }
	      #else
	      //Send the appropriate matrix
	      snd.Addr(u+LexGauge(x,i)*GAUGE_LEN);
	      //Send the transmission, prepare to receive
	      snd.StartTrans();rcv.StartTrans();
	      //Complete the send and receive
	      snd.TransComplete();rcv.TransComplete();
	      #endif
	      //Copy to uc_nl from the received matrix
	      if(cps::GJP.Bc(i)==cps::BND_CND_GPARITY && cps::GJP.NodeCoor(i) == 0) //received data from node behind
		TransCopy((uc_nl[2*i+1]+non_local_count[2*i+1])->mat, rcv_mat); //do transpose instead of dagger
	      else
		DagCopy((uc_nl[2*i+1]+non_local_count[2*i+1])->mat, rcv_mat);
	      non_local_count[i*2+1]++;
	      if (non_local_count[i*2]>non_local_chi[i*2])
		fprintf(stderr,"%s:non_local_count[%d](%d)>non_local_chi[%d](%d)\n",
			fname,2*i,non_local_count[2*i],2*i,non_local_chi[2*i]);
	    } 
	    else if(x[i] == size[i]-1 && cps::GJP.Bc(i) == cps::BND_CND_GPARITY && local[i]){
	      nei[i] = (x[i]+1)%size[i];
	      //local comms in G-parity direction. Link pulled across G-parity boundary from upper-edge needs conjugating as well as daggering
	      //Also need the minus sign for the fermion G-parity flavour twist
	      //printf("For PT from x=(%d,%d,%d,%d) to y=(%d,%d,%d,%d), -(U^dag)^*(x) stored\n",x[0],x[1],x[2],x[3],nei[0],nei[1],nei[2],nei[3]);

	      TransCopy((uc_l_gpbound[2*i+1]+local_gpbound_count[2*i+1])->mat, u+LexGauge(x,i)*GAUGE_LEN);
	      //for(int j=0;j<18;j++) (uc_l_gpbound[2*i+1]+local_gpbound_count[2*i+1])->mat[j]*=-1;
	      
	      local_gpbound_count[i*2+1]++;
	    }
	    else {
	      //Calculate the appropriate neighbor coordinate on the local volume
	      nei[i] = (x[i]+1)%size[i];
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
	
	//Allocate memory for all msg_mem
	//Initialize msg_mem but for receive only

      msg_mem[hop-1][i*2] = new QMP_msgmem_t;
      *(msg_mem[hop-1][i*2]) = QMP_declare_msgmem((void*)rcv_buf[i], hop*non_local_chi[i]*VECT_LEN*sizeof(IFloat));

      msg_mem[hop-1][i*2+1] = new QMP_msgmem_t;

      msg_mem_mat[hop-1][i*2] = new QMP_msgmem_t;
      *(msg_mem_mat[hop-1][i*2]) = QMP_declare_msgmem((void*)rcv_buf[i], 3*hop*non_local_chi[i]*VECT_LEN*sizeof(IFloat));

      msg_mem_mat[hop-1][i*2+1] = new QMP_msgmem_t;
      //Receive

      msg_mem2[hop-1][i*2] = new QMP_msgmem_t;
      *(msg_mem2[hop-1][i*2]) = QMP_declare_msgmem((void*)rcv_buf[i], hop*non_local_chi[i]*VECT_LEN*sizeof(IFloat));

      msg_mem2[hop-1][i*2+1] = new QMP_msgmem_t;
      //#endif
      }

      //-----------------------------------------------------------------------
      //Initialize SCUArg to receive fermion fields
      //#ifdef USE_QMP
      //Declare QMP routines where they will be used

//      msg_mem_cb[i*2] = new QMP_msgmem_t;
      msg_mem_cb[i*2] = QMP_declare_msgmem((void*)rcv_buf[i], non_local_chi_cb[i]*VECT_LEN*sizeof(IFloat));
      msg_handle_cb[i*2] = QMP_declare_receive_relative(msg_mem_cb[i*2], i/2, 1-2*(i%2), 0); 


//      msg_mem_cb[i*2+1] = new QMP_msgmem_t;
    if(i%2){
      msg_mem_cb[i*2+1] = QMP_declare_msgmem((void*)snd_buf_cb[i/2], non_local_chi_cb[i]*VECT_LEN*sizeof(IFloat));
      msg_handle_cb[i*2+1] = QMP_declare_send_relative(msg_mem_cb[i*2+1], i/2, 2*(i%2)-1, 0); 
    } else if(i==6){
      msg_mem_cb[i*2+1] = QMP_declare_msgmem((void*)snd_buf_t_cb, non_local_chi_cb[i]*VECT_LEN*sizeof(IFloat));
      msg_handle_cb[i*2+1] = QMP_declare_send_relative(msg_mem_cb[i*2+1], i/2, 2*(i%2)-1, 0); 
    }
    
//      msg_mem_mat_cb[i*2] = new QMP_msgmem_t;
      msg_mem_mat_cb[i*2] = QMP_declare_msgmem((void*)rcv_buf[i], 3*non_local_chi_cb[i]*VECT_LEN*sizeof(IFloat));
      

//      msg_mem_mat_cb[i*2+1] = new QMP_msgmem_t;
      //-----------------------------------------------------------------------
    }
  

  Free(rcv_mat);
//printf("%s done\n",fname);
}

//CPS_END_NAMESPACE
#endif
