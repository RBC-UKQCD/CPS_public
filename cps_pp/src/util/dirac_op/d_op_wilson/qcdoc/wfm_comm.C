/*--------------------------------------------------------------------*/
/* 1/2002                                                             */
/* PAB... wilson fermion dslash communications source code.           */
/*--------------------------------------------------------------------*/
#include <qcdocos/scu_enum.h>
#include <qcdocos/scu_dir_arg.h>
#include <qcdocos.h>
#include <util/wilson.h>
#include <stdlib.h>
/*Needed for cache_touch*/
#include <qcdoc_align.h>
/*
 * PAB
 * Annoying to have all these globals...
 * Need to think carefully about what to put in the wilson class.
 * Would be nice to compartmentalise all this so that external
 * C compiling/linkage is possible, which means keeping C++ 
 * objects out of the wilson class.
 */

SCUDirArgIR DirArg_send_f[4];
SCUDirArgIR DirArg_send_b[4];
SCUDirArgIR DirArg_recv_f[4];
SCUDirArgIR DirArg_recv_b[4];

SCUDirArgIR *DA_f_p[8];
SCUDirArgIR *DA_b_p[8];

SCUDirArgMulti wfm_multi_f;
SCUDirArgMulti wfm_multi_b;

#define WFM_MAX_INSNS 12
#define WFM_SCU_IR_CHAIN

SCUDMAInst *DMA_recv_f[4][WFM_MAX_INSNS];
SCUDMAInst *DMA_recv_b[4][WFM_MAX_INSNS];

/*
 * Should really go into the wilson struct.
 * Can we have more than one d_op in scope at a time?
 */
static long wfm_do_face_copy[ND];
static long wfm_num_insns[ND];

static const SCUDir  plus_dirs[] = { SCU_XP, SCU_YP, SCU_ZP, SCU_TP } ;
static const SCUDir minus_dirs[] = { SCU_XM, SCU_YM, SCU_ZM, SCU_TM } ;

#include <stdio.h>

extern "C" {

/*--------------------------------------------------------------------*/
/* Initialise the scu stuff called by wilson_init()                   */
/*--------------------------------------------------------------------*/

void check_send(Wilson *wp,int pm,int cb);

void wfm_comm_init(Wilson *wp)
{ 
  int mu,cb,pm, insn;
  int lx,ly,lz,lt,sx;
  int block, stride, nblk;
  int offset;
  Float *Chip;
  unsigned long *face_table;

  /*
   * Compile time switch for now - relatively easy to make a runtime.
   * Just compare wfm_num_insns[mu] to WFM_MAX_INSNS.
   * I will be ok to take the slow/safe code paths whenever *any* of the
   * chains exceed WFM_MAX_INSNS.
   */
  lx = wp->local_latt[0];
  ly = wp->local_latt[1];
  lz = wp->local_latt[2];
  lt = wp->local_latt[3];
  sx = lx/2;

#ifdef WFM_SCU_IR_CHAIN

/*--------------------------------------------------------------------*/
/* Send objects:                                                      */
/*                                                                    */
/* Currently incorrect - need to special case Y-face and chain sx     */
/* instructions, of stride sx*PAD_CHI, block PAD_CHI                  */
/* and offset insn*PAD_CHI                                            */
/* Can get away on 2^4 since sx is 1                                  */
/*--------------------------------------------------------------------*/

  for (mu=0 ; mu < ND ; mu ++ ) {

    /*Forwards xfer send buffer*/
      DirArg_send_f[mu].Init((void *)wp->send_f[mu],
			     plus_dirs[mu],
			     SCU_SEND,
			     (HALF_SPINOR_SIZE*sizeof(Float)),
			     wp->nbound[mu],
			     ((PAD_HALF_SPINOR_SIZE-HALF_SPINOR_SIZE )*sizeof(Float)),
			     IR_1
			     );

    /*Backwards xfer send buffer*/
      DirArg_send_b[mu].Init((void *)wp->send_b[mu],
			     minus_dirs[mu],
			     SCU_SEND,
			     (HALF_SPINOR_SIZE*sizeof(Float)),
			     wp->nbound[mu],
			     ((PAD_HALF_SPINOR_SIZE-HALF_SPINOR_SIZE )*sizeof(Float)),
			     IR_1
			     );
  }

/*--------------------------------------------------------------------*/
/* Receive objects:                                                   */
/* For x,y,z we would need to chain                                   */
/*                                                                    */
/*   ly.lz/2                                                          */
/*   lx/2                                                             */
/*   lt                                                               */
/*                                                                    */
/* Block strided moves respectively.                                  */
/* For now just receive into a receive buffer and scatter the face    */
/* manually for these three dimensions.                               */
/*--------------------------------------------------------------------*/


  /*Just face copy for the x face*/
  wfm_do_face_copy[0] = 1 && (! wp->local_comm[0]); 
  wfm_num_insns[0] = 1;

  /*lx/2 chained for y-face recv */
  wfm_do_face_copy[1] = 0 && (! wp->local_comm[1]); 
  wfm_num_insns[1] = wp->local_latt[0]/2;

  /*lt chained for y-face recv */
  wfm_do_face_copy[2] = 0 && (! wp->local_comm[2]); 
  wfm_num_insns[2] = wp->local_latt[3];

  /*1 t-face never face copy */
  wfm_do_face_copy[3] = 0 ;
  wfm_num_insns[3] = 1;

/*--------------------------------------------------------------------*/
/* X-face receive arguments: just one instruction                     */
/*--------------------------------------------------------------------*/
   /*Forwards xfer recv buffer*/
   mu = 0;
   DirArg_recv_f[mu].Init((void *) wp->recv_f[mu],
			     minus_dirs[mu],
			     SCU_REC,
			     (HALF_SPINOR_SIZE*sizeof(Float)),
			     wp->nbound[mu],
			     ((PAD_HALF_SPINOR_SIZE-HALF_SPINOR_SIZE )*sizeof(Float)),
			     IR_1
			     );
   /*Forwards xfer recv buffer*/
   DirArg_recv_b[mu].Init((void *)wp->recv_b[mu],
			     plus_dirs[mu],
			     SCU_REC,
			     (HALF_SPINOR_SIZE*sizeof(Float)),
			     wp->nbound[mu],
			     ((PAD_HALF_SPINOR_SIZE-HALF_SPINOR_SIZE )*sizeof(Float)),
			     IR_1
			     );


/*--------------------------------------------------------------------*/
/* Y-face receive arguments: wp->local_latt[0]/2 instructions         */
/*                                                                    */
/* These are Block  = HALF_SPINOR_SIZE                                */
/*           Stride = PAD_HALF_SPINOR_SIZE*lx*ly/2                    */
/*           Nblk   = lz*lt                                           */ 
/*                                                                    */
/*        Addr   = Chif + ( 4*insn + 1)                  */
/*                                              *PAD_HALF_SPINOR_SIZE */
/*                                                                    */
/*        Addr   = Chib + ( 4*(insn + (ly-1)*lx/2 ) + 1)*/
/*                                              *PAD_HALF_SPINOR_SIZE */
/*                                                                    */
/*        These should come out the same as offsets in the face_table */
/*        array if I have it right.... check this?                    */
/*--------------------------------------------------------------------*/
    mu=1;
    for ( insn = 0 ; insn < wfm_num_insns[mu] ; insn ++){
      block  = HALF_SPINOR_SIZE * sizeof(Float);
      nblk   = lz *lt;
      stride = ((lx*ly*PAD_HALF_SPINOR_SIZE)*2 * sizeof(Float)) -block;

      offset = ( 4* insn  + mu );
      Chip = wp->af  + PAD_HALF_SPINOR_SIZE * offset ;

      /*Cross check my addressing calculation */
      /*
       *      pm = 0 ;
       *      cb = 0 ;
       *
       *      if ( offset != 	wp->face_table[1-cb][pm][mu][insn] ) {
       *	printf("Y+ Face addressing error\n");
       *	printf("Offset %d != %d\n", offset,wp->face_table[1-cb][pm][mu][insn] );
       *	printf("mu = %d, insn = %d\n",mu,insn);
       *	exit(-1);
       *      }
       */
      DMA_recv_f[mu][insn] = (SCUDMAInst *)malloc(sizeof(SCUDMAInst));
      DMA_recv_f[mu][insn]->Init((void *)Chip,
				 block,
				 nblk,
				 stride
				 );
      
      offset = 4*insn + 2*(ly-1)*lx  + mu ;
      Chip = wp->ab  + PAD_HALF_SPINOR_SIZE * offset ;

      /*Cross check my addressing calculation */
      /*
       *      pm = 1 ;
       *      cb = 0 ;
       *      if ( offset != 	wp->face_table[1-cb][pm][mu][insn] ) {
       *	printf("Y- Face addressing error\n");
       *	printf("Offset %d != %d\n", offset,wp->face_table[1-cb][pm][mu][insn*nblk] );
       *	printf("mu = %d, insn = %d\n",mu,insn);
       *	exit(-1);
       *      }
       */

      DMA_recv_b[mu][insn] = (SCUDMAInst *)malloc(sizeof(SCUDMAInst));
      DMA_recv_b[mu][insn]->Init( (void *)Chip,
				  block,
				  nblk,
				  stride
				 );


    }
    DirArg_recv_f[mu].Init(minus_dirs[mu],SCU_REC,DMA_recv_f[mu],wfm_num_insns[mu],IR_1);
    DirArg_recv_b[mu].Init(plus_dirs[mu] ,SCU_REC,DMA_recv_b[mu],wfm_num_insns[mu],IR_1);

/*--------------------------------------------------------------------*/
/* Z-face receive arguments: wp->local_latt[3]   instructions         */
/*                                                                    */
/* These are Block  = HALF_SPINOR_SIZE                                */
/*           Stride = PAD_HALF_SPINOR_SIZE                            */
/*           Nblk   = lx*ly/2                                         */ 
/*                                                                    */
/*        Addr   = Chif + ((insn*lx*ly*lz/2)*4 + 2)                   */
/*                                              *PAD_HALF_SPINOR_SIZE */
/*                                                                    */
/*        Addr   = Chib + ((insn*lx*ly*lz/2 + (lz-1)*ly*lx/2)*4  + 2) */
/*                                              *PAD_HALF_SPINOR_SIZE */
/*                                                                    */
/*        These should come out the same as offsets in the face_table */
/*        array if I have it right.... check this?                    */
/*                                                                    */
/*                                                                    */
/*--------------------------------------------------------------------*/
    mu=2;
    for ( insn = 0 ; insn < wfm_num_insns[mu] ; insn ++){

      block  = HALF_SPINOR_SIZE * sizeof(Float);
      nblk   = (lx *ly)/2;
      stride = (4*PAD_HALF_SPINOR_SIZE * sizeof(Float)) -block;

      offset = ( 2 * ( insn*lx*ly*lz ) + mu );
      Chip = wp->af  + PAD_HALF_SPINOR_SIZE * offset ;

      /*Cross check my addressing calculation */
      /*
       *      pm = 0 ;
       *      cb = 0 ;
       *      if ( offset != 	wp->face_table[1-cb][pm][mu][insn*nblk] ) {
       *	printf("Z+ Face addressing error\n");
       *	printf("Offset %d != %d\n", offset,wp->face_table[1-cb][pm][mu][insn*nblk] );
       *	printf("mu = %d, insn = %d\n",mu,insn);
       *	exit(-1);
       *      }
       */

      DMA_recv_f[mu][insn] = (SCUDMAInst *)malloc(sizeof(SCUDMAInst));
      DMA_recv_f[mu][insn]->Init((void *)Chip,
				 block,
				 nblk,
				 stride
				 );
      
      offset =  2 * ( insn*lx*ly*lz ) +  2*(lz-1)*ly*lx + mu ;
      Chip = wp->ab  + PAD_HALF_SPINOR_SIZE * offset ;

      /*Cross check my addressing calculation */
      /*
       *      pm = 1 ;
       *      cb = 0 ;
       *      if ( offset != 	wp->face_table[1-cb][pm][mu][insn*nblk] ) {
       *	printf("Z- Face addressing error\n");
       *	printf("Offset %d != %d\n", offset,wp->face_table[1-cb][pm][mu][insn*nblk] );
       *	printf("mu = %d, insn = %d\n",mu,insn);
       *	exit(-1);
       *      }
       */

      DMA_recv_b[mu][insn] = (SCUDMAInst *)malloc(sizeof(SCUDMAInst));

      DMA_recv_b[mu][insn]->Init((void *)Chip,
				 block,
				 nblk,
				 stride
				 );

    }
    DirArg_recv_f[mu].Init(minus_dirs[mu],SCU_REC,DMA_recv_f[mu],wfm_num_insns[mu],IR_1);
    DirArg_recv_b[mu].Init(plus_dirs[mu] ,SCU_REC,DMA_recv_b[mu],wfm_num_insns[mu],IR_1);


/*----------------------------------------------------------------------*/
/* T direction - we can receive directly into the two spinor array      */
/* with a single block-stride                                           */
/*----------------------------------------------------------------------*/

    mu = 3;
    pm = 0;
    cb = 0;
    face_table = wp->face_table[1-cb][pm][mu];
    Chip = wp->af  + PAD_HALF_SPINOR_SIZE * face_table[0] ;

    /*Forwards xfer recv buffer*/
    DirArg_recv_f[mu].Init((void *) Chip,
			   minus_dirs[mu],
			   SCU_REC,
			   (HALF_SPINOR_SIZE*sizeof(Float)),
			   wp->nbound[mu],
			   ((4*PAD_HALF_SPINOR_SIZE-HALF_SPINOR_SIZE )*sizeof(Float)),
			   IR_1
			   );
    

    mu = 3;
    pm = 1;
    cb = 0;
    face_table = wp->face_table[1-cb][pm][mu];

    Chip = wp->ab  + PAD_HALF_SPINOR_SIZE * face_table[0] ;
    /*Forwards xfer recv buffer*/
    DirArg_recv_b[mu].Init((void *)Chip,
			   plus_dirs[mu],
			   SCU_REC,
			   (HALF_SPINOR_SIZE*sizeof(Float)),
			   wp->nbound[mu],
			   ((4*PAD_HALF_SPINOR_SIZE-HALF_SPINOR_SIZE )*sizeof(Float)),
			   IR_1
			   );
#else 

  for (mu=0 ; mu < ND ; mu ++ ) {

    /*Forwards xfer send buffer*/
      DirArg_send_f[mu].Init((void *)wp->send_f[mu],
			     plus_dirs[mu],
			     SCU_SEND,
			     (HALF_SPINOR_SIZE*sizeof(Float)),
			     wp->nbound[mu],
			     ((PAD_HALF_SPINOR_SIZE-HALF_SPINOR_SIZE )*sizeof(Float)),
			     IR_1
			     );

    /*Backwards xfer send buffer*/
      DirArg_send_b[mu].Init((void *)wp->send_b[mu],
			     minus_dirs[mu],
			     SCU_SEND,
			     (HALF_SPINOR_SIZE*sizeof(Float)),
			     wp->nbound[mu],
			     ((PAD_HALF_SPINOR_SIZE-HALF_SPINOR_SIZE )*sizeof(Float)),
			     IR_1
			     );

      if ( mu < ND -1 ) {
	/*Forwards xfer recv buffer*/
	DirArg_recv_f[mu].Init((void *) wp->recv_f[mu],
			     minus_dirs[mu],
			     SCU_REC,
			     (HALF_SPINOR_SIZE*sizeof(Float)),
			     wp->nbound[mu],
			     ((PAD_HALF_SPINOR_SIZE-HALF_SPINOR_SIZE )*sizeof(Float)),
			     IR_1
			     );
	/*Forwards xfer recv buffer*/
	DirArg_recv_b[mu].Init((void *)wp->recv_b[mu],
			     plus_dirs[mu],
			     SCU_REC,
			     (HALF_SPINOR_SIZE*sizeof(Float)),
			     wp->nbound[mu],
			     ((PAD_HALF_SPINOR_SIZE-HALF_SPINOR_SIZE )*sizeof(Float)),
			     IR_1
			     );
      }

  }

/*----------------------------------------------------------------------*/
/* T direction - we can receive directly into the two spinor array      */
/* with a single block-stride                                           */
/*----------------------------------------------------------------------*/

    mu = 3;
    pm = 0;
    cb = 0;
    face_table = wp->face_table[1-cb][pm][mu];
    Chip = wp->af  + PAD_HALF_SPINOR_SIZE * face_table[0] ;

    /*Forwards xfer recv buffer*/
    DirArg_recv_f[mu].Init((void *) Chip,
			   minus_dirs[mu],
			   SCU_REC,
			   (HALF_SPINOR_SIZE*sizeof(Float)),
			   wp->nbound[mu],
			   ((4*PAD_HALF_SPINOR_SIZE-HALF_SPINOR_SIZE )*sizeof(Float)),
			   IR_1
			   );
    

    mu = 3;
    pm = 1;
    cb = 0;
    face_table = wp->face_table[1-cb][pm][mu];

    Chip = wp->ab  + PAD_HALF_SPINOR_SIZE * face_table[0] ;
    /*Forwards xfer recv buffer*/
    DirArg_recv_b[mu].Init((void *)Chip,
			   plus_dirs[mu],
			   SCU_REC,
			   (HALF_SPINOR_SIZE*sizeof(Float)),
			   wp->nbound[mu],
			   ((4*PAD_HALF_SPINOR_SIZE-HALF_SPINOR_SIZE )*sizeof(Float)),
			   IR_1
			   );

#endif

  for (mu=0 ; mu < ND ; mu ++ ) {
    DA_f_p[mu] = & DirArg_send_f[mu];
    DA_b_p[mu] = & DirArg_send_b[mu];
    DA_f_p[mu+ND] = & DirArg_recv_f[mu];
    DA_b_p[mu+ND] = & DirArg_recv_b[mu];
  }

  wfm_multi_f.Init(DA_f_p,8);  
  wfm_multi_b.Init(DA_b_p,8);  

  wp->comm_f = & wfm_multi_f;
  wp->comm_b = & wfm_multi_b;
  return;
}

void wfm_comm_forward_start(Wilson *wp)
{
  wp->comm_f->StartTrans();
}

void wfm_comm_forward_complete(Wilson *wp)
{
  wp->comm_f->TransComplete();

}

void wfm_comm_backward_start(Wilson *wp)
{
  wp->comm_b->StartTrans();
}

void wfm_comm_backward_complete(Wilson *wp)
{
  wp->comm_b->TransComplete();
}

void face_scatter(Float *,Float *, Float **,unsigned long);

void wfm_scatter_face(Wilson *wilson_p, int pm, int cb)
{

  Float *rbp;
  Float *face_p;
  Float *Chip;
  Float **face_table;
  int allbound;
  int i;
  int mu;
  int atom;


  allbound = wilson_p->nbound[0] 
           + wilson_p->nbound[1] 
           + wilson_p->nbound[2];

#ifdef WFM_SCU_IR_CHAIN
  if ( pm == 0 ) {

    /*
     * Hack for chained b/s moves on small volumes
     * 
     */
   for ( mu = 0 ; mu < 3; mu ++ ) {

     if ( (!wilson_p->local_comm[mu]) && wfm_do_face_copy[mu] ){

       face_table = (Float **)wilson_p->face_table[1-cb][pm][mu];
       cache_touch(face_table);

       rbp = wilson_p->recv_f[mu];  // recv_f with pm = 0
       cache_touch(rbp);

       Chip = wilson_p->af;

       /*
	* Call the face scatter assembler routine
	*/
       face_scatter(Chip,rbp,face_table,wilson_p->nbound[mu]);
     }
   }

  }  else {

    for ( mu = 0 ; mu < 3; mu ++ ) {

      if ( (!wilson_p->local_comm[mu])&& wfm_do_face_copy[mu] ){
	face_table = (Float **)wilson_p->face_table[1-cb][pm][mu];
        cache_touch(face_table);

	rbp = wilson_p->recv_b[mu];
        cache_touch(rbp);

	Chip = wilson_p->ab;

	/*
	 * Call the face scatter assembler routine
	 */
	face_scatter(Chip,rbp,face_table,wilson_p->nbound[mu]);
      }
    }
  }
#else
  if ( pm == 0 ) {

    face_table = (Float **)wilson_p->face_table[1-cb][pm][0];
    cache_touch(face_table);
    
    rbp = wilson_p->recv_f[0];  // recv_f with pm = 0
    cache_touch(rbp);

    Chip = wilson_p->af;

    /*
     * Call the face scatter assembler routine
     */
    face_scatter(Chip,rbp,face_table,allbound);

  }  else {

    face_table = (Float **)wilson_p->face_table[1-cb][pm][0];
    cache_touch(face_table);

    rbp = wilson_p->recv_b[0];
    cache_touch(rbp);

    Chip = wilson_p->ab;

	/*
	 * Call the face scatter assembler routine
	 */
    face_scatter(Chip,rbp,face_table,allbound);

  }

#endif

  /* DEBUG receive buffers
   *  check_send(wilson_p,pm,cb); 
   */
}

void check_send(Wilson *wilson_p,int pm, int cb)
{
  int mu,atom;
  int bsite, site;
  unsigned long *face_table;
  Float *Chip;
  Float *sbp;
  Float *face_p;
  if ( pm == 0 ) {

   for ( mu = 0 ; mu < 3; mu ++ ) {

     face_table = wilson_p->face_table[1-cb][pm][mu];

     sbp = wilson_p->send_f[mu];  // recv_f with pm = 0

     Chip = wilson_p->af;
     for ( bsite = 0; bsite < wilson_p->nbound[mu] ; bsite ++ ) {
       face_p = (Float *)(Chip + PAD_HALF_SPINOR_SIZE * face_table[bsite] );
       //       for( atom = 0 ; atom < HALF_SPINOR_SIZE; atom++){
       for( atom = 0 ; atom < 1; atom++){
	 if ( face_p[atom] != sbp[atom] ) {
	   printf("mu=%d bsite=%d ivsite=%d %le %le\n",
		  mu,bsite,face_table[bsite],face_p[atom],sbp[atom]);
	 }
       }

       sbp += PAD_HALF_SPINOR_SIZE;

     }
   }
  }  else {

    for ( mu = 0 ; mu < 3; mu ++ ) {

      face_table = wilson_p->face_table[1-cb][pm][mu];
      sbp = wilson_p->send_b[mu];

      Chip = wilson_p->ab;
      for ( bsite = 0; bsite < wilson_p->nbound[mu] ; bsite ++ ) {
	face_p = (Float *)(Chip + PAD_HALF_SPINOR_SIZE * face_table[bsite] );
       //       for( atom = 0 ; atom < HALF_SPINOR_SIZE; atom++){
	for( atom = 0 ; atom < 1; atom++){
	  if ( face_p[atom] != sbp[atom] ) {
	    printf("mu=%d bsite=%d ivsite=%d %le %le\n",
		   mu,bsite,face_table[bsite],face_p[atom],sbp[atom]);
	  }
	}
	sbp += PAD_HALF_SPINOR_SIZE;
      }

    }
  }
}

}




