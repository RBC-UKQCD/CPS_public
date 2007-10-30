#include "wfm.h"
#include "wfm_internal.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#if TARGET == BGP
void spi_init(){}
#endif

#ifdef USE_COMMS_NONE
void wfm::comm_init(void)
{
  for(int mu=0;mu<4;mu++) {
    if ( !local_comm[mu] ) {
     if ( isBoss() ) printf("wfm::comm_init() nonlocal comms requested"
	    " and COMMS_NONE compiled");
     exit(-1);
    }
  }
}

void wfm::comm_end(void)
{
}

/*Need to blast wilson_p->send_buf[pm][mu] into wilson_p->recv_buf[pm][mu]*/
void wfm::comm_start(int cb)
{
}

/*Need to face_scatter wilson_p->recv_buf[pm][mu]*/
void wfm::comm_complete(int cb)
{
}
int wfm::isBoss(void)
{
  return 1;
}
#endif

#ifdef USE_COMMS_FAKE
/*
 * A fake local copying implementation.
 * If ( local_comm[mu]== 1 ) these routines do not get called
 * If ( local_comm[mu]== 0 ) these routines get called but loop back
 */
#include <string.h>

void wfm::comm_init(void)
{
  return;
}

void wfm::comm_end(void)
{
  return;
}

/*Need to blast wilson_p->send_buf[pm][mu] into wilson_p->recv_buf[pm][mu]*/
void wfm::comm_start(int cb)
{
  int idx,atom;
  int base,offset;
  char *Base;
  int stride;

  char *ptr = (char *) two_spinor;

  for ( int mu=0;mu<ND;mu++ ) { 

    if ( !local_comm[mu] ) { 

      // Special case. We receive the Tface directly into the two spinor
      if( mu == 3) {
	
	// This is a little involved. Since we are simulating sends,
	// being sent with a particular stride (stridenormal) but are
	// the receives are done  with some potentially different 
	// stride (stride).
	//
	unsigned int block_size=HALF_SPINOR_SIZE*TwoSpinSize();
	unsigned int stridenormal = PAD_HALF_SPINOR_SIZE*TwoSpinSize();
	
	
	// -----------------------------------------------------------------
	// Receive from -ve direction into directly 
	// -----------------------------------------------------------------
	// Location where we receive the data
	Base = &ptr[face_table[cb][Plus][mu][0]*PAD_HALF_SPINOR_SIZE*TwoSpinSize()];
	
	// Receive stride. (QMP_conventions, stride from the start of 1 block
	//                  to the start of the next block)
	stride=face_table[cb][Plus][mu][1]-face_table[cb][Plus][mu][0];
	stride*=PAD_HALF_SPINOR_SIZE*TwoSpinSize();
	
	// Send nbound[mu] blocks of block size.
	for(int block=0; block < nbound[mu]; block++) {
	  
	  // Send pointer set with "stridenormal"
	  unsigned char* send_ptr =  ((unsigned char *)send_bufs[Plus][mu])
	    + block*stridenormal;
	  
	  // Receive pointer set with "stride"
	  unsigned char *recv_ptr = (unsigned char *)Base + block*stride;
	  
	  // Do the comms
	  memcpy(recv_ptr, send_ptr, block_size);
	}
	
	// ------------------------------------------------------------------
	// Receive from the +ve direction directly into the two spinor
	// ------------------------------------------------------------------
	Base = &ptr[face_table[cb][Minus][mu][0]*PAD_HALF_SPINOR_SIZE*TwoSpinSize()];
	
	// Communicate n-block 
	for(int block=0; block < nbound[mu]; block++) { 
	  
	  // Send with a stride of "stridenormal"
	  unsigned char *send_ptr = ((unsigned char *)send_bufs[Minus][mu])
	    + block*stridenormal;
	  
	  // Receive with a stride of "stride"
	  unsigned char *recv_ptr = (unsigned char *)Base + block*stride;
	  
	  // Do the "comms"
	  memcpy(recv_ptr, send_ptr, block_size);
	}
      } else { 
	
	// Normal receives
	// Just copy one set of buffers to the other
	for (int  mp=0;mp<2;mp++) {	  
	  
	  // Do the "comms"
	  char *rbuf = (char *)recv_bufs[mp][mu];
	  char *sbuf = (char *)send_bufs[mp][mu];
	  for(idx=0;idx<nbound[mu];idx++){
	    int offset = idx *PAD_HALF_SPINOR_SIZE * TwoSpinSize();
	    memcpy(&rbuf[offset],&sbuf[offset],HALF_SPINOR_SIZE*TwoSpinSize()); 
	  }
	  
	}
      }
    }
  }
}

/* Do the comms, and scatter the necessary faces */
void wfm::comm_complete(int cb)
{
  // Tface is now directly received into the two spinor so we 
  // don't scatter it.
  for ( int mu=0;mu<ND-1;mu++ ) { 
    if ( !local_comm[mu] ) { 
      for (int mp=0;mp<2;mp++) {

		
	  face_scatter(two_spinor,
		       recv_bufs[mp][mu],
		       face_table[cb][mp][mu],
		       nbound[mu]);

	
      }
    }
  }

  return;
}
int wfm::isBoss(void)
{
  return 1;
}
#endif

#ifdef USE_COMMS_SCU
/*
 * QCDOC SCU implementation
 *
 */


void wfm::comm_init(void)
{
  static const SCUDir plus_dirs[] = { SCU_XP,SCU_YP,SCU_ZP,SCU_TP};
  static const SCUDir minus_dirs[]= { SCU_XM,SCU_YM,SCU_ZM,SCU_TM};
  int num_op;
  int cb;
  char *Base;
  unsigned block;
  unsigned stridenormal,stride;
  
  /*----------------------------------------------------------------------*/
  /* T direction - we could receive directly into the two spinor array    */
  /* with a single block-stride                                           */
  /* Optimise for this later                                              */
  /*----------------------------------------------------------------------*/
  
  // CPS like allocation of the send/recv ops (BJ: 18/02/05)
  // Except I check that new succeeds
  for(cb=0; cb < Ncb; cb++) {  
    SendOps[cb] = new SCUDirArgIR[NMinusPlus*ND];
    if( SendOps[cb] == 0x0 ) { 
      if ( isBoss() ) printf("wfm::comm_init(): Failed to allocate SCUDirArgIR, SendOps[%d]\n", 	     cb);
      exit(-1);
    } 
    
    RecvOps[cb] = new SCUDirArgIR[NMinusPlus*ND];
    if( RecvOps[cb] == 0x0 ) { 
      if ( isBoss() ) printf("wfm::comm_init(): Failed to allocate SCUDirArgIR, RecvOps[%d]\n",
	     cb);
      exit(-1);
    }
  }
  
  // CPS like alloction of DA_multi (BJ: 18/02/05)
  // except I check that new succeeds
  DA_multi = new SCUDirArgMulti[2];
  if( DA_multi == 0x0 ) { 
    if ( isBoss() ) printf("wfm::comm_init(): Failed to allocate DA_multi\n");
    exit(-1);
  }
  
  LoadDirArgIRs[0] = 1; 
  LoadDirArgIRs[1] = 1;


  // if ( isBoss() ) printf("Initialising SCU for Dslash IR %d %d\n",IR,IR+1);
 
  block=(HALF_SPINOR_SIZE*TwoSpinSize());
  stridenormal=((PAD_HALF_SPINOR_SIZE-HALF_SPINOR_SIZE )*TwoSpinSize());
  
  char *ptr = (char *)two_spinor;
  for(cb=0;cb<Ncb;cb++){

    num_op = 0;

    // These are used in the cases where various faces can be 
    // received directly into the two spinor.
    XBaseAddr[cb][Plus] =&ptr[face_table[cb][Plus][0][0]*PAD_HALF_SPINOR_SIZE*TwoSpinSize()];
    XBaseAddr[cb][Minus]=&ptr[face_table[cb][Plus][0][0]*PAD_HALF_SPINOR_SIZE*TwoSpinSize()];

    
    //  Define the operations for all the non local comms
    //  In every direction we send the boundaries (the conents
    //  of send_bufs). However the receives are a little more cunning.
    //  with the block strided moves, we can receive them directly
    //  into the two spinor array. We have 3 different receive cases,
    //    i) Tface receive (mu==3) block strided receive directly into
    //       the two spinor array
    //   ii) Xface receive -- block strided receive directly into the 
    //       two spinor array. This is experimental and is currently 
    //       deliberately disabled. The disable is done by marking this 
    //       case as mu=8 (which never gets reached since 0 <= mu <= 3).
    //  iii) Normal receives. Block Strided receives into the recv_bufs
    //       array
     
    for(int mu=0;mu<4;mu++){
      if( !local_comm[mu] ) { 
	
	// Send operations: Everybody sends what is in their 
        //  send buffer. Forwards and backwards.
	
	// Send backwards
	// if ( isBoss() ) printf("Snd mu = %d %d-\n",mu,nbound[mu]);
	SendOps[cb][mu].Init(send_bufs[Minus][mu],
			     minus_dirs[mu],
			     SCU_SEND,
			     block,
			     nbound[mu],
			     stridenormal,
			     IR+cb);
	
        // Send forwards
	// if ( isBoss() ) printf("Snd mu = %d %d+\n",mu,nbound[mu]);
	SendOps[cb][mu+4].Init(send_bufs[Plus][mu],
			       plus_dirs[mu],
			       SCU_SEND,
			       block,
			       nbound[mu],
			       stridenormal,
			       IR+cb);
	
        if ( mu == 3 )  {
          // Special case receive for Tface
	  
	  /*
           * Tface can be pulled directly into
           * the array on receive.
           *
           * Block stride:
           */
          // The offset is the same regardless for T
          Base = &ptr[face_table[cb][Plus][mu][0]*PAD_HALF_SPINOR_SIZE*TwoSpinSize()];
          stride=face_table[cb][Plus][mu][1]-face_table[cb][Plus][mu][0];
          stride*=PAD_HALF_SPINOR_SIZE*TwoSpinSize();
          stride-=block;
	  
          RecvOps[cb][mu].Init((Float *)Base,
                               minus_dirs[mu],
                               SCU_REC,
                               block,
                               nbound[mu],
                               stride,
                               IR+cb);
	  
          Base = &ptr[face_table[cb][Minus][mu][0]*PAD_HALF_SPINOR_SIZE*TwoSpinSize()]; 
          RecvOps[cb][mu+4].Init((Float *)Base,
                                 plus_dirs[mu],
                                 SCU_REC,
                                 (HALF_SPINOR_SIZE*TwoSpinSize()),
                                 nbound[mu],
                                 stride,
                                 IR+cb);
	  
        } else { 
	  // These are the ordinary receives (ie non special cases) 
          // Into the recv bufs
	  // if ( isBoss() ) printf("Rcv mu = %d +\n",mu);
	  RecvOps[cb][mu].Init(recv_bufs[Minus][mu],
			       plus_dirs[mu],
			       SCU_REC,
			       block,
			       nbound[mu],
			       stridenormal,
			       IR+cb);
	  
	  // if ( isBoss() ) printf("mu = %d -\n",mu);
	  RecvOps[cb][mu+4].Init(recv_bufs[Plus][mu],
				 minus_dirs[mu],
				 SCU_REC,
				 block,
				 nbound[mu],
				 stridenormal,
				 IR+cb);
        }
	
        // Accumulate the comms ops into a linear array 
        // which will be used to make a single multi--comms operation
        DA_p[cb][num_op++] = &SendOps[cb][mu];
        DA_p[cb][num_op++] = &SendOps[cb][mu+4];
        DA_p[cb][num_op++] = &RecvOps[cb][mu];
        DA_p[cb][num_op++] = &RecvOps[cb][mu+4];
      } // if ( !local_comm )
      
    } // for(mu=... )
    // if ( isBoss() ) printf("Initialised SCU regs\n");
    
    // Create the multiple comms operation
    DA_multi[cb].Init(DA_p[cb],num_op);
    
  } // for cb
  
  return;
}

void wfm::comm_end(void)
{
  // Free up resources which we newed in comms_init
  LoadDirArgIRs[0] = 1;
  LoadDirArgIRs[1] = 1;
  delete [] SendOps[0];
  delete [] SendOps[1];
  delete [] RecvOps[0];
  delete [] RecvOps[1];
  delete [] DA_multi;
  return;
}


extern "C" { 
  int sys_cacheflush(int);
}

void wfm::comm_start(int cb)
{
  //  sys_cacheflush(0);
  DA_multi[cb].StartTrans();
  LoadDirArgIRs[cb]=0;
  return;
}

void wfm::comm_complete(int cb)
{

  // Complete the communications in progress 
  DA_multi[cb].TransComplete(0);

  // sys_cacheflush commented out by Peter both here and in CPS (BJ: 18/02/05)
  //  sys_cacheflush(0);

  /*
   * Optimisation possible....
   * Could allocate contiguous comms buffers and then
   * issue _one_ call to face_scatter.
   *
   * These copies could be eliminated with chained BS moves.
   * But the chain length grows with Nx and Nt for Y and Z faces,
   * and so I don't want to hardwire a limit.
   */

  // Now that the comms is complete, we need to scatter the faces back
  // into the two spinor.
  // We no longer scatter the Tface because it is received directly 
  // into the two spinor
  for ( int mu=0;mu<ND-1;mu++ ) { 
    if ( !local_comm[mu] ) { 
      for (int mp=0;mp<2;mp++) {
	
	face_scatter(two_spinor,
		     recv_bufs[mp][mu],
		     face_table[cb][mp][mu],
		     nbound[mu]);
      }
    }
  }
  return;
}
int wfm::isBoss(void)
{
  if ( UniqueID() == 0 ) return 1; 
  return 0;
}
#endif

#ifdef USE_COMMS_QMP
/*
 * QMP comms implementation
 * If ( local_comm[mu]== 1 ) these routines do not get called
 * If ( local_comm[mu]== 0 ) these routines get called but loop back
 */
                                                                                
                                                                                
void wfm::comm_init(void)
{
  //static const SCUDir plus_dirs[] = { SCU_XP,SCU_YP,SCU_ZP,SCU_TP};
  // static const SCUDir minus_dirs[]= { SCU_XM,SCU_YM,SCU_ZM,SCU_TM};
  
  static const int dims[] = { 0, 1, 2, 3 };
  static const int plus_signs[] = { +1, +1, +1, +1 };
  static const int minus_signs[] = { -1, -1, -1, -1};
  int num_op;
  int cb;
  Float *Base;
  unsigned block;
  unsigned stridenormal,stride;
  
  // if ( isBoss() ) printf("Initialising QMP for Dslash\n");
  
  // Space for 16 handles. In this I will collect the handles for
  // into a linear array for just 1 checkerboard, 
  // from which I will create a multiple handle
  // The multiple handle creation invalidates all its inputs
  // and frees all its handles  so I can simply reuse this 
  // for each checkerboard.
  QMP_msghandle_t tmp_handles[16];
  
                                                                              
  block=(HALF_SPINOR_SIZE*TwoSpinSize());

  // QMP Stride style (from start of block, rather than from end)
  stridenormal=(PAD_HALF_SPINOR_SIZE*TwoSpinSize());
  
  for(cb=0;cb<Ncb;cb++){
    
    num_op = 0;
    
    // These are used in the cases where various faces can be
    // received directly into the two spinor.
    
    XBaseAddr[cb][Plus] =&two_spinor[face_table[cb][Plus][0][0]
				     *PAD_HALF_SPINOR_SIZE];
    
    XBaseAddr[cb][Minus]=&two_spinor[face_table[cb][Plus][0][0]
				     *PAD_HALF_SPINOR_SIZE];
    
    //  Define the operations for all the non local comms
    //  In every direction we send the boundaries (the conents
    //  of send_bufs). However the receives are a little more cunning.
    //  with the block strided moves, we can receive them directly
    //  into the two spinor array. We have 3 different receive cases,
    //    i) Tface receive (mu==3) block strided receive directly into
    //       the two spinor array
    //   ii) Xface receive -- block strided receive directly into the
    //       two spinor array. This is experimental and is currently
    //       deliberately disabled. The disable is done by marking this
    //       case as mu=8 (which never gets reached since 0 <= mu <= 3).
    //  iii) Normal receives. Block Strided receives into the recv_bufs
    //       array
    
    for(int mu=0;mu<4;mu++){
      
      if( !local_comm[mu] ) {

          /* ----------------------------------------------------------*/
          /* Normal Receives to recv_bufs                              */
          /* ----------------------------------------------------------*/
	  
          /* ----------------------------------------------------------*/
          /* RECV_FROM: POSITIVE direction, into recv_bufs[Minus][mu]  */
          /* ----------------------------------------------------------*/
	  /* Equivalent SCU Code:
                 RecvOps[cb][mu].Init(recv_bufs[Minus][mu],
                                 plus_dirs[mu],
                                 SCU_REC,
                                 block,
                                 nbound[mu],
                                 stridenormal,
                                 IR+cb);
           */

          // if ( isBoss() ) printf("Recv mu = %d +\n",mu);
	  // Declare the msgmemory for it:
          if (stridenormal==block)
          qmp_recv_ops_msgmem_t[cb][mu] = QMP_declare_msgmem(
                                           (void *)recv_bufs[Minus][mu],
                                           block*nbound[mu]);
          else
          qmp_recv_ops_msgmem_t[cb][mu] = QMP_declare_strided_msgmem(
                                           (void *)recv_bufs[Minus][mu],
                                           block,
                                           nbound[mu],
                                           stridenormal);
      	  // Check it 
          if ( qmp_recv_ops_msgmem_t[cb][mu] == 0 ) {
            if ( isBoss() ) printf("QMP_declare_strided_msgmem returned NULL QMP_msgmem_t: qmp_recv_ops_msgmem_t[%d][%d]\n", cb, mu);
        
            exit(-1);
          } 
	  
	  // Declare the handle for it: 
          tmp_handles[num_op] = QMP_declare_receive_relative(
                                             qmp_recv_ops_msgmem_t[cb][mu],
                                             dims[mu],
                                             plus_signs[mu],
                                             0);

	  // Check it
          if( tmp_handles[num_op] == 0x0 ) {
            if ( isBoss() ) printf("QMP_declare_recv_relative returned NULL QMP_msghandle_t: tmp_handle[%d], cb=%d, mu = %d direction=%d\n", num_op, cb, mu, plus_signs[mu]);
            exit(-1);
          }
	  
	  // Increment the counter
	  num_op++;
	  
          /* ----------------------------------------------------------*/
          /* RECV_FROM: POSITIVE DIRECTION DECLARED                    */
          /* ----------------------------------------------------------*/
	  
          /* ----------------------------------------------------------*/
          /* RECV_FROM: NEGATIVE direction, into recv_bufs[Plus][mu]   */
          /* ----------------------------------------------------------*/
          /* Equivalent to SCU code:

              RecvOps[cb][mu+4].Init(recv_bufs[Plus][mu],
                             minus_dirs[mu],
                             SCU_REC,
                             block,
                             nbound[mu],
                             stridenormal,
                             IR+cb);
          */	
                                                                             
          // if ( isBoss() ) printf("Recv mu = %d -\n",mu);
	  // Declare the message memory from recv_bufs[Plus]
          if (stridenormal==block)
          qmp_recv_ops_msgmem_t[cb][mu+4] = QMP_declare_msgmem(
                                              (void *)recv_bufs[Plus][mu],
                                              block*nbound[mu]);
          else
          qmp_recv_ops_msgmem_t[cb][mu+4] = QMP_declare_strided_msgmem(
                                              (void *)recv_bufs[Plus][mu],
                                              block,
                                              nbound[mu],
                                              stridenormal);

	  // Check it
          if ( qmp_recv_ops_msgmem_t[cb][mu+4] == 0x0 ) {
            if ( isBoss() ) printf("QMP_declare_strided_msgmem returned NULL QMP_msgmem_t: qmp_recv_ops_msgmem_t[%d][%d]\n", cb, mu+4);
                                                                                
            exit(-1);
          }
      
	  // Declare the message handle for it
          tmp_handles[num_op]  = QMP_declare_receive_relative(
					    qmp_recv_ops_msgmem_t[cb][mu+4],
                                            dims[mu],
                                            minus_signs[mu],
                                            0);
          // Check it:
          if( tmp_handles[num_op] == 0x0 ) {
            if ( isBoss() ) printf("QMP_declare_recv_relative returned NULL QMP_msghandle_t: tmp_handle[%d], cb=%d, mu = %d, direction = %d\n", num_op, cb, mu, minus_signs[mu]);
            exit(-1);
          }
	  
          num_op++; 


	/* ========================================================= */
        /* Send Operations. As normal in all directions              */
        /* ========================================================= */
	
        /* ----------------------------------------------------------*/
        /* SEND TO: NEGATIVE direction, from send_bufs[Minus][mu]    */
        /* ----------------------------------------------------------*/
	/* Equivalent to SCU Code:
	   SendOps[cb][mu].Init(send_bufs[Minus][mu],
                             minus_dirs[mu],
                             SCU_SEND,
                             block,
                             nbound[mu],
                             stridenormal,
                             IR+cb);
        */
	// if ( isBoss() ) printf("Snd mu = %d %d-\n",mu,nbound[mu]);

        // Declare a message memory from the send buf
        if (stridenormal==block)
        qmp_send_ops_msgmem_t[cb][mu] = QMP_declare_msgmem(
					(void *)send_bufs[Minus][mu],
				        block* nbound[mu]);
	else
	qmp_send_ops_msgmem_t[cb][mu] = QMP_declare_strided_msgmem(
					(void *)send_bufs[Minus][mu],
				        block,  
					nbound[mu],
					stridenormal);

        // Check it -- Assuming that QMP_msgmem_t is really a pointer 
        if ( qmp_send_ops_msgmem_t[cb][mu] == 0x0 ) { 
	   if ( isBoss() ) printf("QMP_declare_strided_msgmem returned NULL QMP_msgmem_t: qmp_send_ops_msgmem_t[%d][%d]\n", cb, mu);
           exit(-1);
        }

	// Declare an individual message handle for the communication
	tmp_handles[num_op] = QMP_declare_send_relative(
							qmp_send_ops_msgmem_t[cb][mu],
							dims[mu],
							minus_signs[mu],
							0);
	
        // Check the message handle has been created 
	if( tmp_handles[num_op] == 0x0 ) {
	  if ( isBoss() ) printf("QMP_declare_send_relative returned NULL QMP_msghandle_t: tmp_handles[%d], cb=%d, mu = %d, direction=%d\n", num_op, cb, mu, minus_signs[mu]);
	  exit(-1);
        }
	// Increase tmp_handles index
        num_op++;
	
	/* ----------------------------------------------------------*/
	/* SEND TO NEGATIVE DIRECTION DECLARED                       */
	/* ----------------------------------------------------------*/
	
        /* ----------------------------------------------------------*/
        /* SEND TO: POSITIVE direction, from send_bufs[Plus][mu]     */
        /* ----------------------------------------------------------*/
        /* Equivalent to SCU Code:
	       SendOps[cb][mu+4].Init(send_bufs[Plus][mu],
	                              plus_dirs[mu],
	                              SCU_SEND,
                                      block,
                                      nbound[mu],
                                      stridenormal,
                                      IR+cb);
         */
        // if ( isBoss() ) printf("Snd mu = %d %d+\n",mu,nbound[mu]);
	// Declare msgmem from sendbuf
        if (stridenormal==block)
        qmp_send_ops_msgmem_t[cb][mu+4] = QMP_declare_msgmem(
					    (void *)send_bufs[Plus][mu],
					     block* nbound[mu]);
	else
        qmp_send_ops_msgmem_t[cb][mu+4] = QMP_declare_strided_msgmem(
					     (void *)send_bufs[Plus][mu],
					     block,
					     nbound[mu],
					     stridenormal);
	
	// Check it
	if ( qmp_send_ops_msgmem_t[cb][mu+4] == 0x0 ) {
	  if ( isBoss() ) printf("QMP_declare_strided_msgmem returned NULL QMP_msgmem_t: qmp_send_ops_msgmem_t[%d][%d]\n", cb, mu+4);
           exit(-1);
        }
	
 	// Declare message handle
        tmp_handles[num_op]  = QMP_declare_send_relative(
							 qmp_send_ops_msgmem_t[cb][mu+4],
							 dims[mu],
							 plus_signs[mu],
							 0);
	// Check it
        if( tmp_handles[num_op] == 0x0 ) {
	  if ( isBoss() ) printf("QMP_declare_send_relative returned NULL QMP_msghandle_t: tmp_handle[%d], cb=%d, mu = %d, direction=%d\n", num_op, cb, mu, plus_signs[mu]);
	  exit(-1);
        }
	
	// Increment counter/tmp_handles index
	num_op++;

        /* ----------------------------------------------------------*/
        /* SEND TO POSITIVE DIRECTION DECLARED                       */
        /* ----------------------------------------------------------*/

   
      } // if ( !local_comm )
      
    } // for(mu=... )
    
    /* -------------------------------------------------------------------*
     * Declare a multiple handle for all the comms on this checkerboard   *
     * Equivalent to SCU Code:                                            *
     *        DA_multi[cb].Init(DA_p[cb],num_op);                         *
     * -------------------------------------------------------------------*
     * Note:    
     *    Our handles are already in the linear array tmp_handles      
     *    we now call QMP_declare_multiple() 
     *
     *  According to the QMP spec:
     * "Use of declare_multiple invalidates the individual I/O handles 
     *  and subsequent operations on the individual I/O handles are undefined.
     *  QMP_declare_multiple will free the individual I/O handles so that the 
     *  user will not have to do this." 
     *  I take this to mean that I don't have to free the individual 
     *  tmp_handles. However I will need to free the msgmem_t-s in comm_end() 
     *  as well as these multi_handles                                       
     * ------------------------------------------------------------------*/
    // if ( isBoss() ) printf("num op = %d\n", num_op);
    if( num_op > 0 ) { 
      qmp_multi_handles[cb] = QMP_declare_multiple(tmp_handles, num_op);
    }

    
  } // for cb
  
  nonlocal_comms = false;
  for(int mu=0; mu < ND; mu++) {
    if ( !local_comm[mu] ) {
      nonlocal_comms = true;
    }
  } 
  return;
}

void wfm::comm_end(void)
{
  /*-----------------------------------------------------------------------* 
   * I need to free the multiple msghandles                                *
   * The individual handles (tmp_handles) used should have been freed      *
   * when the multiple was created and according to the spec I don't need  *
   * to free them                                                          *
   *---------------------------------------------------------------------- */
  if( nonlocal_comms ) {
    for(int cb=0; cb < Ncb; cb++) { 
      QMP_free_msghandle(qmp_multi_handles[cb]);
    }
  }
  /*---------------------------------------------------------------------*
   * I need to free the associated msgmem_t-s here.                      *
   *---------------------------------------------------------------------*/
  for(int cb=0; cb < Ncb; cb++) { 
    for(int mu=0; mu < ND; mu++) { 
      if (!local_comm[mu]) { 
	QMP_free_msgmem(qmp_send_ops_msgmem_t[cb][mu]);
	QMP_free_msgmem(qmp_send_ops_msgmem_t[cb][mu+4]);
	QMP_free_msgmem(qmp_recv_ops_msgmem_t[cb][mu]);
	QMP_free_msgmem(qmp_recv_ops_msgmem_t[cb][mu+4]);
      }
    }
  }
  
  return;
}
                                                                                
void wfm::comm_start(int cb)
{

  /* -------------------------------------------------------------------*
   * I need to start the multiples communications here                  * 
   * -------------------------------------------------------------------*/
  int lcb = (cb + base_parity) & 1;
  if (nonlocal_comms) { 
    QMP_status_t ret_val = QMP_start(qmp_multi_handles[lcb]);

    // Check
    if( ret_val != QMP_SUCCESS ) { 
      if ( isBoss() ) printf("QMP_start failed for cb %d. Ret_val was %d\n", cb, ret_val);
      exit(-1);
    }
  }
  return;
}

void wfm::comm_complete(int cb)
{
  /* ------------------------------------------------------------------*
   * I need to wait for the multiple comms to finish here and then     *
   * scatter some faces                                                *
   * ------------------------------------------------------------------*/
  int lcb = (cb + base_parity) & 1;
  if ( nonlocal_comms ) { 
    QMP_status_t ret_val = QMP_wait(qmp_multi_handles[lcb]);
    
  
    // Check
    if( ret_val != QMP_SUCCESS ) { 
      if ( isBoss() ) printf("QMP_wait returned unsuccessful value: %d\n", ret_val);
      exit(-1);
    }
  }

#if 0
  for ( int pm = 0;pm<2;pm++ ) {
    for ( int mu = 0 ; mu < 4 ; mu++) {
      unsigned long half_spinor_words = PAD_HALF_SPINOR_SIZE * nbound[mu];
      Float *tmp = recv_bufs[pm][mu];
      for(unsigned long i=0; i<half_spinor_words;i++)
        if( fabs(tmp[i])>0.1)
          printf("recv_buf[%d][%d][%d]=%5f\n",pm,mu,i,tmp[i]);
    }
  }
#endif
  
  // Scatter faces if necessary
  for ( int mu=0;mu<ND;mu++ ) {
    if ( !local_comm[mu] ) {
      for (int mp=0;mp<2;mp++) {

	  face_scatter(two_spinor,
		       recv_bufs[mp][mu],
		       face_table[lcb][mp][mu],
		       nbound[mu]);
	/*
        Float *rbp = recv_bufs[mp][mu];
        for (int i = 0; i < nbound[mu]; i++) {
          Float *face_p = (Float *)(two_spinor + PAD_HALF_SPINOR_SIZE * face_table[lcb][mp][mu][i] );
          for(int atom = 0 ; atom < HALF_SPINOR_SIZE; atom++){
            face_p[atom] = rbp[atom];
          }
          rbp += PAD_HALF_SPINOR_SIZE;
        }
	*/
      }
    }
  }
  return;
}
int wfm::isBoss(void)
{
  return QMP_is_primary_node();
}

#endif
void wfm::face_scatter(Float *TwoSpinor,
		       Float *RcvBuf, 
		       unsigned long*FaceTable,
		       unsigned long n)
{
  WFM_BGL=true;
  if ( WFM_BGL && SloppyPrecision ) {
    bgl_s_face_scatter(TwoSpinor,RcvBuf,FaceTable,n);
  } else if ( WFM_BGL ) { 
    bgl_face_scatter(TwoSpinor,RcvBuf,FaceTable,n);
  } else if ( SloppyPrecision ) { 
    qcdoc_s_face_scatter(TwoSpinor,RcvBuf,FaceTable,n);
  } else {
    qcdoc_face_scatter(TwoSpinor,RcvBuf,FaceTable,n);

  }
}

