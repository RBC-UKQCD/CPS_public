#include <util/wfm.h>
#include "wfm_internal.h"
#include <stdio.h>
#include <stdlib.h>

#ifdef USE_COMMS_NONE
void wfm::comm_init(void)
{
  for(int mu=0;mu<4;mu++) {
    if ( !local_comm[mu] ) {
     printf("wfm::comm_init() nonlocal comms requested"
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
#endif

#ifdef USE_COMMS_FAKE
/*
 * A fake local copying implementation.
 * If ( local_comm[mu]== 1 ) these routines do not get called
 * If ( local_comm[mu]== 0 ) these routines get called but loop back
 */
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
  /*
   * This is equivalent of a block strided move, 
   * block HALF_SPINOR_SIZE, stride (PAD_HALF_SPINOR_SIZE - BLOCK)
   */
  for ( int mu=0;mu<ND;mu++ ) { 
    if ( !local_comm[mu] ) { 
      for (int  mp=0;mp<2;mp++) {
	base = 0;
	for(idx=0;idx<nbound[mu];idx++){
	  for(atom=0;atom<HALF_SPINOR_SIZE;atom++){
	    recv_bufs[mp][mu][base] = send_bufs[mp][mu][base]; 
	    base++;
	  }
	  base +=PAD_HALF_SPINOR_SIZE - HALF_SPINOR_SIZE;
	}
      }
    }
  }
}
/*Need to face_scatter wilson_p->recv_buf[pm][mu]*/
void wfm::comm_complete(int cb)
{
  for ( int mu=0;mu<ND;mu++ ) { 
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
#endif

#ifdef USE_COMMS_SCU
/*
 * A fake local copying implementation.
 * If ( local_comm[mu]== 1 ) these routines do not get called
 * If ( local_comm[mu]== 0 ) these routines get called but loop back
 */


void wfm::comm_init(void)
{
  static const SCUDir plus_dirs[] = { SCU_XP,SCU_YP,SCU_ZP,SCU_TP};
  static const SCUDir minus_dirs[]= { SCU_XM,SCU_YM,SCU_ZM,SCU_TM};
  int num_op;
  int cb;
  Float *Base;
  unsigned block;
  unsigned stridenormal,stride;
  
/*----------------------------------------------------------------------*/
/* T direction - we could receive directly into the two spinor array    */
/* with a single block-stride                                           */
/* Optimise for this later                                              */
/*----------------------------------------------------------------------*/
  printf("Initialising SCU for Dslash IR %d %d\n",IR,IR+1);
 
  block=(HALF_SPINOR_SIZE*sizeof(Float));
  stridenormal=((PAD_HALF_SPINOR_SIZE-HALF_SPINOR_SIZE )*sizeof(Float));

  for(cb=0;cb<Ncb;cb++){

    num_op = 0;
    XBaseAddr[cb][Plus] =&two_spinor[face_table[cb][Plus][0][0]
				   *PAD_HALF_SPINOR_SIZE];
    XBaseAddr[cb][Minus]=&two_spinor[face_table[cb][Plus][0][0]
				   *PAD_HALF_SPINOR_SIZE];

    for(int mu=0;mu<4;mu++){
      if( !local_comm[mu] ) { 


	printf("Snd mu = %d %d+\n",mu,nbound[mu]);
	SendOps[cb][mu].Init(send_bufs[Minus][mu],
			     minus_dirs[mu],
			     SCU_SEND,
			     block,
			     nbound[mu],
			     stridenormal,
			     IR+cb);

	printf("Snd mu = %d %d-\n",mu,nbound[mu]);
	SendOps[cb][mu+4].Init(send_bufs[Plus][mu],
			       plus_dirs[mu],
			       SCU_SEND,
			       block,
			       nbound[mu],
			       stridenormal,
			       IR+cb);

	if ( mu == 3 ) { 

	/*
	 * Tface can be pulled directly into
	 * the array on receive.
	 *
	 * Block stride:
	 */
	// The offset is the same regardless for T
	  Base = &two_spinor[face_table[cb][Plus][mu][0]*PAD_HALF_SPINOR_SIZE];
	  stride=face_table[cb][Plus][mu][1]-face_table[cb][Plus][mu][0];
	  stride*=PAD_HALF_SPINOR_SIZE*sizeof(Float);
	  stride-=block;

	  printf("Tplus %x %d %d\n",Base,block,stride);

	  RecvOps[cb][mu].Init(Base,
			       minus_dirs[mu],
			       SCU_REC,
			       block,
			       nbound[mu],
			       stride,
			       IR+cb);

	  Base = &two_spinor[face_table[cb][Minus][mu][0]*PAD_HALF_SPINOR_SIZE];
	  printf("Tminus %x %d %d\n",Base,block,stride);
	  
	  RecvOps[cb][mu+4].Init(Base,
			     plus_dirs[mu],
			     SCU_REC,
			     (HALF_SPINOR_SIZE*sizeof(Float)),
			     nbound[mu],
			     stride,
			     IR+cb);

	} else if ( mu==8 ) { 
	/*
	 * Xface can be pulled directly into
	 * the array on receive.... This guy depends
	 * on the checker board, and the base address for
	 * these two instructions only must be reloaded
	 * on every application of dslash (unless we take
	 * two instructionregisters, one for odd one for even)
	 */
	  Base = &two_spinor[face_table[cb][Plus][mu][0]*PAD_HALF_SPINOR_SIZE];
	  stride = (face_table[cb][Plus][mu][1] - face_table[cb][Plus][mu][0]) 
	    * PAD_HALF_SPINOR_SIZE * sizeof(Float);
	  stride -=block;
	  printf("Face table %d %d\n",
		 face_table[cb][Plus][mu][0],
		 face_table[cb][Plus][mu][1]
		 );
	  printf("Xplus %x %d %d\n",Base,block,stride);
	  RecvOps[cb][mu].Init(Base,
			   plus_dirs[mu],
			   SCU_REC,
			   block,
			   nbound[mu],
			   stride,
			   IR+cb);

	  Base = &two_spinor[face_table[cb][Minus][mu][0]*PAD_HALF_SPINOR_SIZE];
	  printf("Xminus %x %d %d\n",Base,block,stride);
	  RecvOps[cb][mu+4].Init(Base,
				 minus_dirs[mu],
				 SCU_REC,
				 block,
				 nbound[mu],
				 stride,
				 IR+cb);

	} else { 

	  printf("mu = %d +\n",mu);
	  RecvOps[cb][mu].Init(recv_bufs[Minus][mu],
			       plus_dirs[mu],
			       SCU_REC,
			       block,
			       nbound[mu],
			       stridenormal,
			       IR+cb);

	  printf("mu = %d -\n",mu);
	  RecvOps[cb][mu+4].Init(recv_bufs[Plus][mu],
			     minus_dirs[mu],
			     SCU_REC,
			     block,
			     nbound[mu],
			     stridenormal,
			     IR+cb);
	}


	DA_p[cb][num_op++] = &SendOps[cb][mu];
	DA_p[cb][num_op++] = &SendOps[cb][mu+4];
	DA_p[cb][num_op++] = &RecvOps[cb][mu];
	DA_p[cb][num_op++] = &RecvOps[cb][mu+4];
      }
      
    }
    printf("Initialised SCU regs\n");

    DA_multi[cb].Init(DA_p[cb],num_op);

  }

  return;
}

void wfm::comm_end(void)
{
  return;
}

void wfm::comm_start(int cb)
{
  //  sys_cacheflush(0);
  DA_multi[cb].StartTrans();
  return;
}

extern "C" { 
int sys_cacheflush(int);
}
void wfm::comm_complete(int cb)
{
  DA_multi[cb].TransComplete();
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
  for ( int mu=0;mu<3;mu++ ) { 
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
#endif

#ifdef USE_COMMS_QMP
#error QMP NOT SUPPORTED YET
#endif
