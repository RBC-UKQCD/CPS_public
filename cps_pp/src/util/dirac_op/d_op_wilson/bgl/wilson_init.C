#include<config.h>
CPS_START_NAMESPACE
/****************************************************************************/
/* 10/6/96                                                                  */
/*                                                                          */
/* wilson_int:                                                              */
/*                                                                          */
/* This routine performs all initializations needed before wilson func      */
/* are called. It sets the addressing related arrays and reserves memory    */
/* for the needed temporary buffers. It only needs to be called             */
/* once at the begining of the program (or after a wilson_end call)         */
/* before any number of calls to wilson funcs are made.                     */
/*                                                                          */
/* The half spinor temporaries are allocated in write_throug                */
/* store-with-out-allocate memory via direct pointer and not malloc.        */
/*                                                                          */
/****************************************************************************/

/*--------------------------------------------------------------------------*/
/* Include header files                                                     */
/*--------------------------------------------------------------------------*/
CPS_END_NAMESPACE
#include<util/wilson.h>
#include<util/gjp.h>
#include<util/smalloc.h>
#include<util/verbose.h>
#include<util/error.h>
#include <sys/bgl/bgl_sys.h>
#include <sys/bgl/bgl_sys_all.h>
CPS_START_NAMESPACE

void wilson_set_sloppy(bool sloppy) {}



/*--------------------------------------------------------------------------*/
/* # defs
/*--------------------------------------------------------------------------*/
//#define SCRATCH_MEM_BASE BGL_MEM_SCRATCHPAD_BASE

/*--------------------------------------------------------------------------*/
/* Externals                                                                */
/*--------------------------------------------------------------------------*/
int wfm_max_numchunk;
int wfm_numchunk[8];
IFloat **wfm_send_ad;
IFloat **wfm_recv_ad;
IFloat *wfm_s_start[8];
IFloat *wfm_r_start[8];
unsigned long wfm_blklen[8];
unsigned long wfm_numblk[8];
unsigned long wfm_stride[8];

/*--------------------------------------------------------------------------*/
/* Function declarations                                                    */
/*--------------------------------------------------------------------------*/
void wfm_sublatt_pointers(int slx, 
			  int sly, 
			  int slz, 
			  int slt, 
			  int slatt_eo, 
			  Wilson *wilson_p);

/*=========================================================================*/
/* wilson_init:                                                            */
/*=========================================================================*/
int wilson_initted=0;
void wilson_init(Wilson *wilson_p)  /* pointer to Wilson type structure    */
{
  char *cname = " ";
  char *fname = "wilson_init(Wilson*)";
  VRB.Func(cname,fname);

  int spinor_words;             /* size of the spinor field on the         */
				/* sublattice checkerboard                 */

  int half_spinor_words;        /* size of the spin-projected "half_spinors*/
                                /* on the sublattice checkerboard including*/
                                /* the communications padding              */

  int slx;                          /* x-direction size of node sublattice */
  int sly;                          /* y-direction size of node sublattice */
  int slz;                          /* z-direction size of node sublattice */
  int slt;                          /* t-direction size of node sublattice */
  int slatt_eo;                     /* =0/1 if the sublattice is even/odd. */
  int size;

  int   mu, i, j, k, dir;
  IFloat *send_ad, *receive_ad;

  IFloat *offset_swoa;
  unsigned scratch_base;
  unsigned scratch_size;
  int pir;


/*--------------------------------------------------------------------------*/
/* Set sublattice direction sizes                                           */
/*--------------------------------------------------------------------------*/
  slx = GJP.XnodeSites();
  sly = GJP.YnodeSites();
  slz = GJP.ZnodeSites();
  slt = GJP.TnodeSites();

/*--------------------------------------------------------------------------*/
/* Determine if the sublattice is even or odd from the node coordinates     */
/* If (px,py,pz,pt) are the coordinates of the node and the node            */
/* mesh is of size (nx,ny,nz,nt) then a node is even/odd if its             */
/* lexicographical number =  px + nx * ( py + ny * ( pz + nz * ( pt )))     */
/* is even/odd.                                                             */
/*--------------------------------------------------------------------------*/
/* A runtime system function is needed here to determine (px,py,pz,pt) and  */
/* (nx,ny,nz,nt). For now we set slat_eo = 0 which is a safe choice if      */
/* slx,sly,slz,slt are all even.                                            */
/* ??? */
  slatt_eo = 0;

/*--------------------------------------------------------------------------*/
/* Reserve memory for the node sublattice pointers                          */
/*--------------------------------------------------------------------------*/
  size = 40*sly*slz*slt*sizeof(int);
  wilson_p->ptr = (int *) smalloc(size);
  if( wilson_p->ptr == 0)
    ERR.Pointer(cname,fname, "ptr");
  VRB.Smalloc(cname,fname,
	      "ptr", wilson_p->ptr, size);

/*--------------------------------------------------------------------------*/
/* Set the node sublattice pointers                                         */
/*--------------------------------------------------------------------------*/
  wfm_sublatt_pointers(slx, sly, slz, slt, slatt_eo, wilson_p);

/*--------------------------------------------------------------------------*/
/* Reserve memory for 2  temporary spinors (nedded by m, mdag and mdagm)    */
/* Memory is allocated into two consecutive chunks. The above routines      */
/* will have to set the pointer of the first temporary to the base          */
/* and of the second in the middle.                                         */
/* WARNING: valid for "even" sublattices only                               */
/*--------------------------------------------------------------------------*/
  spinor_words = 2 * SPINOR_SIZE * wilson_p->vol[0];

  wilson_p->spinor_tmp = (IFloat *) smalloc(spinor_words*sizeof(IFloat));
  if(wilson_p->spinor_tmp == 0)
    ERR.Pointer(cname,fname, "spinor_tmp");
  VRB.Smalloc(cname,fname,
	      "spinor_tmp", wilson_p->spinor_tmp, spinor_words*sizeof(IFloat));




#if BLRTS_SUPERV == 1 
/*--------------------------------------------------------------------------*/
/* Get the pointer to an area of memory with TLB L1 setting of              */
/* write-through, store without allocate (swoa).                            */
/* Each core is assigned its own area.                                      */
/*--------------------------------------------------------------------------*/
  pir = rts_get_processor_id();
  if(rts_get_scratchpad_window((void **) &scratch_base, &scratch_size) != 0){
    ERR.General("wilson_init", " ", "Failed to get pointer to swoa memory\n");
  }
  offset_swoa = (IFloat *) (scratch_base + pir*(1024)*(512));

/*--------------------------------------------------------------------------*/
/* Set the pointers for the 4 forward and 4 backward spin projected half    */ 
/* spinors in the swoa memory area.                                         */
/*--------------------------------------------------------------------------*/  
  wilson_p->af[0] = offset_swoa;
  half_spinor_words = HALF_SPINOR_SIZE * wilson_p->padded_subgrid_vol[0];
  offset_swoa = offset_swoa + half_spinor_words;

  wilson_p->af[1] = offset_swoa;
  half_spinor_words = HALF_SPINOR_SIZE * wilson_p->padded_subgrid_vol[1];
  offset_swoa = offset_swoa + half_spinor_words;

  wilson_p->af[2] = offset_swoa;
  half_spinor_words = HALF_SPINOR_SIZE * wilson_p->padded_subgrid_vol[2];
  offset_swoa = offset_swoa + half_spinor_words;

  wilson_p->af[3] = offset_swoa;
  half_spinor_words = HALF_SPINOR_SIZE * wilson_p->padded_subgrid_vol[3];
  offset_swoa = offset_swoa + half_spinor_words;


  wilson_p->ab[0] = offset_swoa;
  half_spinor_words = HALF_SPINOR_SIZE * wilson_p->padded_subgrid_vol[0];
  offset_swoa = offset_swoa + half_spinor_words;

  wilson_p->ab[1] = offset_swoa;
  half_spinor_words = HALF_SPINOR_SIZE * wilson_p->padded_subgrid_vol[1];
  offset_swoa = offset_swoa + half_spinor_words;

  wilson_p->ab[2] = offset_swoa;
  half_spinor_words = HALF_SPINOR_SIZE * wilson_p->padded_subgrid_vol[2];
  offset_swoa = offset_swoa + half_spinor_words;

  wilson_p->ab[3] = offset_swoa;
  half_spinor_words = HALF_SPINOR_SIZE * wilson_p->padded_subgrid_vol[3];
  offset_swoa = offset_swoa + half_spinor_words;

#else


#if BLRTS_CX == 1
/*--------------------------------------------------------------------------*/
/* Get the pointer to an area of memory with TLB L1 setting of              */
/* write-through, store without allocate (swoa).                            */
/* Each core is assigned its own area.                                      */
/*--------------------------------------------------------------------------*/
  if(rts_get_writethrough_scratchpad_window((void **) &scratch_base, &scratch_size) != 0){
    ERR.General("wilson_init", " ", "Failed to get pointer to swoa memory\n");
  }
  offset_swoa = (IFloat *) scratch_base;

/*--------------------------------------------------------------------------*/
/* Set the pointers for the 4 forward and 4 backward spin projected half    */ 
/* spinors in the swoa memory area.                                         */
/*--------------------------------------------------------------------------*/  
  wilson_p->af[0] = offset_swoa;
  half_spinor_words = HALF_SPINOR_SIZE * wilson_p->padded_subgrid_vol[0];
  offset_swoa = offset_swoa + half_spinor_words;

  wilson_p->af[1] = offset_swoa;
  half_spinor_words = HALF_SPINOR_SIZE * wilson_p->padded_subgrid_vol[1];
  offset_swoa = offset_swoa + half_spinor_words;

  wilson_p->af[2] = offset_swoa;
  half_spinor_words = HALF_SPINOR_SIZE * wilson_p->padded_subgrid_vol[2];
  offset_swoa = offset_swoa + half_spinor_words;

  wilson_p->af[3] = offset_swoa;
  half_spinor_words = HALF_SPINOR_SIZE * wilson_p->padded_subgrid_vol[3];
  offset_swoa = offset_swoa + half_spinor_words;


  wilson_p->ab[0] = offset_swoa;
  half_spinor_words = HALF_SPINOR_SIZE * wilson_p->padded_subgrid_vol[0];
  offset_swoa = offset_swoa + half_spinor_words;

  wilson_p->ab[1] = offset_swoa;
  half_spinor_words = HALF_SPINOR_SIZE * wilson_p->padded_subgrid_vol[1];
  offset_swoa = offset_swoa + half_spinor_words;

  wilson_p->ab[2] = offset_swoa;
  half_spinor_words = HALF_SPINOR_SIZE * wilson_p->padded_subgrid_vol[2];
  offset_swoa = offset_swoa + half_spinor_words;

  wilson_p->ab[3] = offset_swoa;
  half_spinor_words = HALF_SPINOR_SIZE * wilson_p->padded_subgrid_vol[3];
  offset_swoa = offset_swoa + half_spinor_words;

#else

/*--------------------------------------------------------------------------*/
/* Reserve memory for the 4 forward and 4 backward spin projected half      */ 
/* spinors using standar smalloc.                                           */
/*--------------------------------------------------------------------------*/
  for(i=0; i<4; i++){
    half_spinor_words = HALF_SPINOR_SIZE * wilson_p->padded_subgrid_vol[i];

    wilson_p->af[i] = (IFloat *) smalloc (cname,fname,
		"af[i]", half_spinor_words*sizeof(IFloat));

    wilson_p->ab[i] = (IFloat *) smalloc (cname,fname,
		"ab[i]", half_spinor_words*sizeof(IFloat));
  }

#endif

#endif


/*--------------------------------------------------------------------------*/
/* BGL specific                                                             */
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/* Set the total number of chuncks for each direction                       */
/*--------------------------------------------------------------------------*/
  for(mu=0; mu<ND; ++mu){
    wfm_numchunk[mu]   = wilson_p->comm_numblk[mu]*wilson_p->comm_blklen[mu];
    wfm_numchunk[4+mu] = wfm_numchunk[mu];
  }

/*--------------------------------------------------------------------------*/
/* find the maximum numchunk from all 8 dirs.                               */
/*--------------------------------------------------------------------------*/
  wfm_max_numchunk = 0;
  for(dir=0; dir<8; dir++){
    if(wfm_numchunk[dir] > wfm_max_numchunk){
      wfm_max_numchunk = wfm_numchunk[dir];
    }
  }

/*--------------------------------------------------------------------------*/
/* Reserve memory for the send and receive address arrays                   */
/*--------------------------------------------------------------------------*/
  wfm_send_ad = (IFloat **) smalloc(cname,fname,"wfm_send_ad",8*wfm_max_numchunk*sizeof(int));

  wfm_recv_ad = (IFloat **) smalloc (cname,fname,
	      "wfm_recv_ad", 8*wfm_max_numchunk*sizeof(int));


/*--------------------------------------------------------------------------*/
/* Set the send address for all +directions -> x+,y+,z+,t+                  */
/* Set the recv address for all -directions -> x-,y-,z-,t-                  */
/*--------------------------------------------------------------------------*/
  for(mu=0; mu<ND; ++mu){
    wfm_s_start[mu]= wilson_p->comm_offset[mu] + wilson_p->ab[mu];
    wfm_r_start[mu]= wilson_p->ab[mu];
    wfm_numblk[mu] = wilson_p->comm_numblk[mu];
    wfm_blklen[mu] = wilson_p->comm_blklen[mu]*BLOCK*sizeof(IFloat);
    wfm_stride[mu] = (wilson_p->comm_stride[mu]-1)*sizeof(IFloat)+wfm_blklen[mu];
#if 0
    if(!UniqueID())
    printf("Node 0: wfm %d: %p %p %d %d %d\n", mu, wfm_s_start[mu],
          wfm_r_start[mu],wfm_numblk[mu],wfm_blklen[mu],wfm_stride[mu]);
#endif
  }
  for(mu=0; mu<ND; ++mu){
    wfm_s_start[4+mu]= wilson_p->af[mu];
    wfm_r_start[4+mu]= wilson_p->comm_offset[mu] + wilson_p->af[mu];
    wfm_numblk[4+mu] = wilson_p->comm_numblk[mu];
    wfm_blklen[4+mu] = wilson_p->comm_blklen[mu]*BLOCK*sizeof(IFloat);
    wfm_stride[4+mu] = (wilson_p->comm_stride[mu]-1)*sizeof(IFloat)+wfm_blklen[mu];
  }

  for(mu=0; mu<ND; ++mu)
    {
      k = 0;
      send_ad = wilson_p->comm_offset[mu] + wilson_p->ab[mu];
      receive_ad = 0 + wilson_p->ab[mu];
      for(i=0; i<wilson_p->comm_numblk[mu]; ++i) { 

	for(j=0; j<wilson_p->comm_blklen[mu]; ++j) {
	  wfm_send_ad[mu+8*k] = send_ad;
	  wfm_recv_ad[mu+8*k] = receive_ad;
//    if(!UniqueID()) printf("Node 0: wfm_ad[%d][%d]: %p %p \n",mu,k,send_ad,receive_ad);
	  k = k+1;
	  //sends to +dir (2*mu), receives from - dir (2*mu+1)
	  //getMinusData(receive_ad, send_ad, BLOCK, mu);
	  send_ad    = send_ad    + BLOCK;
	  receive_ad = receive_ad + BLOCK;
	}	 
	send_ad = send_ad + wilson_p->comm_stride[mu] - 1;
	receive_ad = receive_ad + wilson_p->comm_stride[mu] - 1;
      }
    }

/*--------------------------------------------------------------------------*/
/* Set the send address for all -directions -> x-,y-,z-,t-                  */
/* Set the recv address for all +directions -> x+,y+,z+,t+                  */
/*--------------------------------------------------------------------------*/

  for(mu=0; mu<ND; ++mu)
    {
      k = 0;
      receive_ad = wilson_p->comm_offset[mu] + wilson_p->af[mu];
      send_ad = 0 + wilson_p->af[mu];
      for(i=0; i<wilson_p->comm_numblk[mu]; ++i) { 

	for(j=0; j<wilson_p->comm_blklen[mu]; ++j) {
	  wfm_send_ad[4+mu+8*k] = send_ad;
	  wfm_recv_ad[4+mu+8*k] = receive_ad;

	  k = k+1;
	  //sends to -dir (2*mu+1), receives from + dir (2*mu)
	  //getPlusData(receive_ad, send_ad, BLOCK, mu);
	  send_ad    = send_ad    + BLOCK;
	  receive_ad = receive_ad + BLOCK;
    if(!UniqueID()) printf("Node 0: wfm_ad[%d][%d]: %p %p \n",mu+4,k,send_ad,receive_ad);
	}
	send_ad = send_ad + wilson_p->comm_stride[mu] - 1;
	receive_ad = receive_ad + wilson_p->comm_stride[mu] - 1;
      }
   }

   wilson_initted=1;
}


CPS_END_NAMESPACE
