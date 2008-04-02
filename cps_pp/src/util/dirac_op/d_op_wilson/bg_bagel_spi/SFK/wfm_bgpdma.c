#define _WFM_BGPDMA_C_
/****************************************************************************
All routines relevant to the Wilson fermion kernel and  
BlueGene/P dma communications are defined in this file.
                                                   P. Vranas 7/7/06
Corrected, parallelized in the time direction and converted to VN mode
and userspace, fixed timing and added communication routines for single 
precision and get{Plus,Minus}Data routines. 
                                                   S. Krieg  07/2007
*****************************************************************************/

/*----------------------------------------------------------------------------
In this routine the BG/P dma is initialized so that it can do the Wilson 
fermion kernel persistent communications. 

This routine must be called after the Wilson structure has been initialized
by the wilson_init.c routine. The half spinor addresses (af, ab) are 
initialized so that each core has its own addresses.

The needed block-strided-moves (bsm) for each direction are in the input 
argument structure wilson_p. 

All bsm are "expanded" into descriptors (one for each block transfer). 

The direction numbering used is:
x+, x-, y+, y-, z+, z-, t+, t- ==> 0, 1, 2, 3, 4, 5, 6, 7

There are:
4 forward  communication directions (send to + and  receive from -) and
4 backward communication directions (send to - and  receive from +) 

From the 4 forward  comms 3 are performed using the torus netwok 
in order to communicate outside the chip. The other forward comm is
between cores within the chip (local).

From the 4 backward comms 3 are performed using the torus netwok 
in order to communicate outside the chip. The other backward comm is
between cores within the chip (local).

The forward comms are performed at a different time than the
backward comms (they overlap with different parts of the code). 
Therefore it makes sense that one counter is assigned for all forward 
comms, torus and local, and another counter for all backward comms,
torus and local.

Therefore the notation fwd/bwd reffers to these two comm phases and
also indicates that both the physics and hardware comms are along the
+/- directions respectively.

In order to have all 4 fwd(bwd) directions busy one needs one fifo 
per direction for a total of 8 fifos. The descriptors in these fifos
should communicate roughly the same amount of data so that one 
descriptor does not occupy the dma with a large transfer along one 
direction and as a result the other directions are idle (such a 
situation would effectively serialize the transfers and would 
result to large transfer times)

Therefore the notation tors/local refers to the dimension that the 
hardware communication takes place. By convention:
hardware dimension = x, y, z ==> 0, 1, 2 ==> torus
hardware dimension =       t ==>       3 ==> local

The physics lattice convention is:
physics dimension = x, y, z, t ==> 0, 1, 2, 3

Each core must have its own set of fifos therefore each core is assigned
a separate group of counters and fifos. For BGP there are 4 cores and
8 groups of fifos per core.

For each of the 4 groups the assignements are:
Fifo 0:  fwd torus comms x+ \
Fifo 2:  fwd torus comms y+  \
                              ==> counter 0
Fifo 4:  fwd torus comms z+  /
Fifo 6:  fwd local comms t+ /


Fifo 1:  bwd torus comms x- \
Fifo 3:  bwd torus comms y-  \
                              ==> counter 1
Fifo 5:  bwd torus comms z-  /
Fifo 7:  bwd local comms t- /

Therefore:

The even numbered fifos correspond to fwd (+) direction and all use
counter 0.

The odd numbered fifos correspond to bwd (-) direction and all use
counter 1.

These assignements are followed both for the injection and reception
counters and fifos.
----------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h> //shmem
#include <unistd.h>   //shmem
#include <fcntl.h>    //shmem
#include <sys/types.h>//shmem
#include <unistd.h>   //shmem
#include <spi/bgp_SPI.h>
#include <errno.h>
#include <memory.h>
#include "param_init.h"
#include "wilson.h"
#include "wfm_bgpdma.h"
#include "Layout.h"
#include "GlobBarrier.h"
#include "NodeBarrier.h"
#include "Mem.h" 
#include "Helpers.h"
#include "Errors.h"
extern int errno;


/*----------------------------------------------------------------------------
  # Defines
  ----------------------------------------------------------------------------*/

#define DMA_PACKET_HINT_NULL  (0x00) // looks official but it is not
#define SFK_MAX_D 4
#define NumInjFifos 10 // MUST be >=8
#define AuxInjFifoTorus 8  // the aux inj fifo, default is 9
#define AuxInjFifoLocal 9  // the aux inj fifo, default is 9
#define BGP_WFM_NUM_DEFAULT_DESCRIPTORS 128
#define NumInjFifos_s 8
/* shared memory stuff */
#define MAX_SHARED_SIZE 4096 /* hardcoded for the moment */
#define SHM_FILE "/shm-file"
/* a verbose print */
#define Vprintf if (!bgp_wfm_hw_my_rank) printf

/*----------------------------------------------------------------------------
  Globals
  ----------------------------------------------------------------------------*/
int bgp_wfm_initialized = 0;
int bgp_wfm_hw_my_rank = -1;
int bgp_wfm_hw_my_coords[4] = {-1, -1, -1, -1};
int bgp_wfm_hw_num_nodes[4] = {-1, -1, -1, -1};
int bgp_wfm_dim_map[4] = {0, 1, 2, 3};
/* shared memory stuff */
char *bgp_wfm_shmptr = NULL;
/* The memory fifos */
static unsigned char *IFifo[NumInjFifos];
static unsigned char *IFifo_s[NumInjFifos_s];
/* injection counter group id */
static int inj_c_grp;
/* reception counter group id */
static int rec_c_grp[NUM_CORES];
/* injection fifo group id */
static int inj_f_grp;
/* injection counter group structure */
static DMA_CounterGroup_t inj_counter_group;
/* reception counter group structure */
static DMA_CounterGroup_t rec_counter_group;
/* injection fifo group structure */
static DMA_InjFifoGroup_t inj_fifo_group;
/* injection counters for forward/backward, for both torus and local */
static int inj_counter_id_fwd;
static int inj_counter_id_bwd;
/* reception counters for forward/backward, for both torus and local */
static int rec_counter_id_fwd[NUM_CORES];
static int rec_counter_id_bwd[NUM_CORES];
/* Bytes to communicate allong all 4 fwd dims. It is the same for bwd */
static int comm_bytes;
static int comm_bytes_s;
/* the nearest neighbours */
static int nn_x[8];
static int nn_y[8];
static int nn_z[8];
static int nn_t[8];


/*----------------------------------------------------------------------------
  Functions
  ----------------------------------------------------------------------------*/

void wfm_bgpdma_init(Wilson *wilson_p)
{
    char *fname = " wfm_bgpdma_init";
    int rc, Pid, i, mu, hdim, blk;
    int iPid;
    _BGP_Personality_t mem_pers;
    _BGP_Personality_t *pers=&mem_pers;
    int dim_map[4];
    int my_x, my_y, my_z, my_t;
    int nodes_x, nodes_y, nodes_z, nodes_t;
    int max_x, max_y, max_z, max_t;
    int inj_fifo_num_descriptors[8];
    int inj_fifo_size_in_bytes[NumInjFifos];
    int subgroups[2];
    int num_inj_fifos;
    int num_inj_fifos_s;
    int inj_fifo_ids[NumInjFifos+NumInjFifos_s];
    unsigned short int inj_priorities[NumInjFifos+NumInjFifos_s];
    unsigned short int inj_locals[NumInjFifos+NumInjFifos_s];
    unsigned char ts_inj_maps[NumInjFifos+NumInjFifos_s];
    uint32_t tail;
    int blk_byte_size;
    uint32_t send_off, recv_off;
    extern uint32_t comm_off[4];
    /* Descriptor for preparation and immediate storage into fifo */
    DMA_InjDescriptor_t inj_desc;
    /* this is initialized in wfm_bgpdma_init. contains no hints, actually... */
    unsigned char hints[6];
    int fd; // file id thing for shared mem

    /*--------------------------------------------------------------------------
      Set the communication dimension map:
      hdim = dim_map[pdim] where pdim is the physics lattice dimension 
      and hdim is the hardware dimention that it corresponds to. 
      This must be set by the physics system.
      For now it is set here with hdim = pdim.
      The dimension numbering is:
      x, y, z, t ==> 0, 1, 2, 3
      Hardware fimensions 0, 1, 2, are assigned to torus comms and 3 is 
      assigned to local comms.
      --------------------------------------------------------------------------*/
    dim_map[0] = bgp_wfm_dim_map[0];
    dim_map[1] = bgp_wfm_dim_map[1];
    dim_map[2] = bgp_wfm_dim_map[2];
    dim_map[3] = bgp_wfm_dim_map[3];

    /*--------------------------------------------------------------------------
      Get the processor id
      --------------------------------------------------------------------------*/
    Pid = GET_PID;

    /*--------------------------------------------------------------------------
      Get the personality
      --------------------------------------------------------------------------*/
    Kernel_GetPersonality(pers, sizeof(_BGP_Personality_t));

    /*--------------------------------------------------------------------------
      Set this node coordinates and direction sizes
      --------------------------------------------------------------------------*/
    my_x        = pers->Network_Config.Xcoord;
    my_y        = pers->Network_Config.Ycoord;
    my_z        = pers->Network_Config.Zcoord;
    my_t        = Pid;
    nodes_x     = pers->Network_Config.Xnodes;
    nodes_y     = pers->Network_Config.Ynodes;
    nodes_z     = pers->Network_Config.Znodes;
    nodes_t     = NUM_CORES;
    max_x       = (nodes_x - 1);
    max_y       = (nodes_y - 1);
    max_z       = (nodes_z - 1);
    max_t       = (nodes_t - 1);
    /* set the globals */
    bgp_wfm_hw_my_rank          = my_t + 
	nodes_t * ( my_z + nodes_z * ( my_y + nodes_y * my_x ) ); 
    bgp_wfm_hw_my_coords[DIR_X] = my_x;
    bgp_wfm_hw_my_coords[DIR_Y] = my_y;
    bgp_wfm_hw_my_coords[DIR_Z] = my_z;
    bgp_wfm_hw_my_coords[DIR_T] = my_t;
    bgp_wfm_hw_num_nodes[DIR_X] = nodes_x;
    bgp_wfm_hw_num_nodes[DIR_Y] = nodes_y;
    bgp_wfm_hw_num_nodes[DIR_Z] = nodes_z;
    bgp_wfm_hw_num_nodes[DIR_T] = nodes_t;
  
    /*--------------------------------------------------------------------------
      Set the x,y,z,t coordinates of each nearest neighbor along the 
      8 directions (0,1,2,3,4,5,6,7 = +x,-x,+y,-y,+z,-z,+t,-t)
      --------------------------------------------------------------------------*/
    nn_x[0] = (my_x + nodes_x + 1) % nodes_x;
    nn_x[1] = (my_x + nodes_x - 1) % nodes_x;
    nn_x[2] =  my_x;
    nn_x[3] =  my_x;
    nn_x[4] =  my_x;
    nn_x[5] =  my_x;
    nn_x[6] =  my_x;
    nn_x[7] =  my_x;

    nn_y[0] =  my_y;
    nn_y[1] =  my_y;
    nn_y[2] = (my_y + nodes_y + 1) % nodes_y;
    nn_y[3] = (my_y + nodes_y - 1) % nodes_y;
    nn_y[4] =  my_y;
    nn_y[5] =  my_y;
    nn_y[6] =  my_y;
    nn_y[7] =  my_y;

    nn_z[0] =  my_z;
    nn_z[1] =  my_z;
    nn_z[2] =  my_z;
    nn_z[3] =  my_z;
    nn_z[4] = (my_z + nodes_z + 1) % nodes_z;
    nn_z[5] = (my_z + nodes_z - 1) % nodes_z;
    nn_z[6] =  my_z;
    nn_z[7] =  my_z;

    nn_t[0] =  my_t;
    nn_t[1] =  my_t;
    nn_t[2] =  my_t;
    nn_t[3] =  my_t;
    nn_t[4] =  my_t;
    nn_t[5] =  my_t;
    nn_t[6] = (my_t + nodes_t + 1) % nodes_t;
    nn_t[7] = (my_t + nodes_t - 1) % nodes_t;

    GLOB_BARRIER;
    Vprintf("%s: DMA comm layout set up.\n", fname);
    fflush(stdout);
    GLOB_BARRIER;

    /*--------------------------------------------------------------------------
      Set the hint bits for nearest neighbor comms along the 6 torus directions
      --------------------------------------------------------------------------*/
    hints[0] = DMA_PACKET_HINT_NULL;
    hints[1] = DMA_PACKET_HINT_NULL;
    hints[2] = DMA_PACKET_HINT_NULL;
    hints[3] = DMA_PACKET_HINT_NULL;
    hints[4] = DMA_PACKET_HINT_NULL;
    hints[5] = DMA_PACKET_HINT_NULL;


    /*--------------------------------------------------------------------------
      Various initializations
      REMEMBER: DMA is a compute node global thing. even in VN mode each
      counter/fifo id, either group or idividual thing can be 
      allocated only once.
      --------------------------------------------------------------------------*/
    /* Injection counter group id */
    inj_c_grp = Pid;
    /* Reception counter. For local copies need neighbours counter group id */
    for (iPid=0; iPid<NUM_CORES; iPid++){
	rec_c_grp[iPid] = iPid;
    }
    inj_f_grp = Pid;

    /* Forward/Backward counter subgroups*/
    subgroups[0] = 0;
    subgroups[1] = 1;

    /* Set injection counter ids */
    inj_counter_id_fwd = 0;
    inj_counter_id_bwd = 1;

    /* Set reception counter ids. For local copies need neighbours counter id */
    for (iPid=0; iPid<NUM_CORES; iPid++){
	rec_counter_id_fwd[iPid] = 0;
	rec_counter_id_bwd[iPid] = 1;
    }
  

    /*--------------------------------------------------------------------------
      Allocate injection counters
      --------------------------------------------------------------------------*/
    rc = DMA_CounterGroupAllocate(DMA_Type_Injection,
				  inj_c_grp,           //  grp,
				  2,                        //  num_subgroups (fwd, bwd)
				  subgroups,                //  *subgroups, (fwd=0, bwd=1)
				  0,                        //  target,
				  NULL,                     //  Kernel_InterruptHandler_t handler,
				  NULL,                     //  *handler_parm,
				  0,                        //  Kernel_InterruptGroup_t   interruptGroup,
				  &inj_counter_group); //  DMA_CounterGroup_t *cg_ptr

    if ( rc ==0) { 
	// BGPDMA_VPRINT(" wfm_bgpdma: SUCCESS for DMA_CounterGroupAllocate DMA_Type_Injection\n"); 
    }
    else { 
	printf(" wfm_bgpdma: FAILURE for DMA_CounterGroupAllocate DMA_Type_Injection\n");
	KillJob(rc);
    }

    /*--------------------------------------------------------------------------
      Allocate reception counters
      --------------------------------------------------------------------------*/
    rc = DMA_CounterGroupAllocate(DMA_Type_Reception,
				  rec_c_grp[Pid],           //  grp,
				  2,                        //  num_subgroups, (fwd, bwd)
				  subgroups,                //  *subgroups, (fwd=0, bwd=1)
				  0,                        //  target,
				  NULL,                     //  Kernel_InterruptHandler_t handler,
				  NULL,                     //  *handler_parm,
				  0,                        //  Kernel_InterruptGroup_t   interruptGroup,
				  &rec_counter_group); //  DMA_CounterGroup_t *cg_ptr

    if ( rc ==0) { 
	// BGPDMA_VPRINT(" wfm_bgpdma: SUCCESS for DMA_CounterGroupAllocate DMA_Type_Reception\n"); 
    }
    else { 
	printf(" wfm_bgpdma: FAILURE for DMA_CounterGroupAllocate DMA_Type_Reception \n");
	KillJob(rc);
    }

    GLOB_BARRIER;
    Vprintf("%s: Counters allocated.\n", fname);
    fflush(stdout);
    GLOB_BARRIER;

    // Allocate the 8 injection fifos
    // 0 -> torus x+
    // 1 -> torus x-
    // 2 -> torus y+
    // 3 -> torus y-
    // 4 -> torus z+
    // 5 -> torus z-
    // 6 -> local t+
    // 7 -> local t-
    //--------------------------------------------------------------------
    num_inj_fifos = NumInjFifos;
    num_inj_fifos_s = NumInjFifos_s; // the fifos for the single prec kernel

    inj_fifo_ids[0] = 0;
    inj_fifo_ids[1] = 1;
    inj_fifo_ids[2] = 2;
    inj_fifo_ids[3] = 3;
    inj_fifo_ids[4] = 4;
    inj_fifo_ids[5] = 5;
    inj_fifo_ids[6] = 6;
    inj_fifo_ids[7] = 7;

    inj_priorities[0] = 0;
    inj_priorities[1] = 0;
    inj_priorities[2] = 0;
    inj_priorities[3] = 0;
    inj_priorities[4] = 0;
    inj_priorities[5] = 0;
    inj_priorities[6] = 0;
    inj_priorities[7] = 0;

    inj_locals[0] = 0;
    inj_locals[1] = 0;
    inj_locals[2] = 0;
    inj_locals[3] = 0;
    inj_locals[4] = 0;
    inj_locals[5] = 0;
    inj_locals[6] = 1;
    inj_locals[7] = 1;

    // torus: can inject in any non-priority ts fifo
    ts_inj_maps[0] = 0x80;   
    ts_inj_maps[1] = 0x40;   
    ts_inj_maps[2] = 0x20;   
    ts_inj_maps[3] = 0x10;   
    ts_inj_maps[4] = 0x08;   
    ts_inj_maps[5] = 0x04;   

    // local
    ts_inj_maps[6] = 0x00;   
    ts_inj_maps[7] = 0x00;   

    // additional fifos for comms
    for (i=8; i<num_inj_fifos; i++){
	inj_fifo_ids[i] = i;
	inj_priorities[i] = 0;
	inj_locals[i] = i%2;
	if (i%2){
	    ts_inj_maps[i] = 0x00;
	}else{
	    ts_inj_maps[i] = ts_inj_maps[i%6];
	}
    }
    for (i=num_inj_fifos; i<num_inj_fifos+num_inj_fifos_s; i++){
	inj_fifo_ids[i] = i;
	inj_priorities[i] = 0;
	inj_locals[i] = inj_locals[i-num_inj_fifos];
	ts_inj_maps[i] = ts_inj_maps[i-num_inj_fifos];
    }

    rc = DMA_InjFifoGroupAllocate(inj_f_grp,                //  injection fifo grp,
				  num_inj_fifos+num_inj_fifos_s, //  number of injection fifos
				  inj_fifo_ids,             //  injection fifo ids
				  inj_priorities,           //  injection fifo priorities
				  inj_locals,               //  local injection fifos
				  ts_inj_maps,              //  torus injection maps
				  &inj_fifo_group);    //  return structure
    if ( rc == 0) {
	// BGPDMA_VPRINT(" wfm_bgpdma: SUCCESS for DMA_InjFifoGroupAllocate\n");
    }
    else {
	printf(" wfm_bgpdma: Injection Fifo  Allocate FAILURE rc = %d\n",rc);
	KillJob(rc);
    }

    GLOB_BARRIER;
    Vprintf("%s: Injection Fifo Group allocated.\n", fname);
    fflush(stdout);
    GLOB_BARRIER;


    /*--------------------------------------------------------------------------
      Calculate the number of descriptors to be stored in the injection fifo
      (the number is the same for + and - directions).
      --------------------------------------------------------------------------*/

    /* Loop over the lattice directions */
    for(mu=0; mu<SFK_MAX_D; mu++) {

	/* Set the hardware dimension hdim */
	hdim = dim_map[mu];

	/* Add the number of descriptors needed to comunicate this phase. The 
	   result is the same for the - direction. */
	inj_fifo_num_descriptors[(2*mu)]   = wilson_p->comm_numblk[mu];
	inj_fifo_num_descriptors[(2*mu)+1] = wilson_p->comm_numblk[mu];
    
    }


    /*--------------------------------------------------------------------------
      Calculate the injection fifos size in bytes. Each descriptor is 32 Bytes 
      and an additional num Bytes must be allocated for internal dma purposes.
      --------------------------------------------------------------------------*/
    for(i=0; i<8; i++){
	inj_fifo_size_in_bytes[i] = inj_fifo_num_descriptors[i]
	    * DMA_FIFO_DESCRIPTOR_SIZE_IN_BYTES // this is 32 
	    + 2 * DMA_MIN_INJ_FIFO_SIZE_IN_BYTES;// 2times 4*256 just to be safe
    }

    /* The additional fifos for the generic comms */
    for(i=8; i<NumInjFifos; i++){
	inj_fifo_size_in_bytes[i] = BGP_WFM_NUM_DEFAULT_DESCRIPTORS
	    * DMA_FIFO_DESCRIPTOR_SIZE_IN_BYTES // this is 32
	    + 2 * DMA_MIN_INJ_FIFO_SIZE_IN_BYTES;// 2times 4*256 just to be safe
    }
    /* The fifos for single precision persitent comms */
    for(i=NumInjFifos; i<NumInjFifos+NumInjFifos_s; i++){
	inj_fifo_size_in_bytes[i] = inj_fifo_size_in_bytes[i-NumInjFifos];
    }
  

    /*--------------------------------------------------------------------------
      Calculate the total number of bytes that will be communicated
      along all fwd dimensions. The total includes both torus and local 
      transfers. The total is the same for fwd and bwd.
      --------------------------------------------------------------------------*/
    comm_bytes = 0;

    /* Loop over the 4 lattice directions */
    for(mu=0; mu<SFK_MAX_D; mu++) {

	/* Set the hardware dimension hdim */
	hdim = dim_map[mu];

	/* If the grid does not end in the + direction of this dimension add 
	   the bytes needed to comunicate this phase. The result is the same 
	   for the - direction. */
	comm_bytes = comm_bytes +
	    (  wilson_p->comm_numblk[mu] 
	       * wilson_p->comm_blklen[mu]
	       * BLOCK
	       * IFLOAT_IN_BYTES );
    }
    /* single prec */
    comm_bytes_s = comm_bytes/2;

    /*--------------------------------------------------------------------------
      Allocate memory for the 8 injection fifos that will hold the descriptors.
      --------------------------------------------------------------------------*/
  
    for (i=0; i<NumInjFifos; i++){
	IFifo[i] = 
	    (unsigned char *) alloc_L1_aligned(inj_fifo_size_in_bytes[i]);
    }
    for (i=NumInjFifos; i<NumInjFifos+NumInjFifos_s; i++){
	IFifo_s[i-NumInjFifos] = 
	    (unsigned char *) alloc_L1_aligned(inj_fifo_size_in_bytes[i]);
    }

    GLOB_BARRIER;
    Vprintf("%s: Memory allocation for njection fifos successfull.\n", fname);
    fflush(stdout);
    GLOB_BARRIER;

    /*--------------------------------------------------------------------------
      Initialize the 8 injection fifos
      --------------------------------------------------------------------------*/
    rc = 0;

    for (i=0; i<NumInjFifos; i++){
	rc = DMA_InjFifoInitById(&inj_fifo_group,         // fifo group structure
				 i,                       // fifo id
				 (void *) IFifo[i],     // fifo start address
				 (void *) IFifo[i],     // fifo header address
				 (void *) &IFifo[i][inj_fifo_size_in_bytes[i]]); 
	if (rc) {
	    printf(" wfm_bgpdma: FAILURE for DMA_InjFifoInitById, "
		   "Ififo[%i], rc=%i \n", i, rc);
	    KillJob(rc);
	}
	// fifo end address (+32 to ensure it is a multiple of 32)
    }
    for (i=NumInjFifos; i<NumInjFifos+NumInjFifos_s; i++){
	rc = DMA_InjFifoInitById(&inj_fifo_group,         // fifo group structure
				 i,                       // fifo id
				 (void *) IFifo_s[i-NumInjFifos],     // fifo start address
				 (void *) IFifo_s[i-NumInjFifos],     // fifo header address
				 (void *) &IFifo_s[i-NumInjFifos][inj_fifo_size_in_bytes[i]]); 
	if (rc) {
	    printf(" wfm_bgpdma: FAILURE for DMA_InjFifoInitById, "
		   "Ififo[%i], rc=%i \n", i, rc);
	    KillJob(rc);
	}
	// fifo end address (+32 to ensure it is a multiple of 32)
    }

    /* activate AUX fifos */
    for (i=8; i<NumInjFifos; i++){
	DMA_InjFifoSetActivateById(&inj_fifo_group, // fifo group structure
				   i);                // fifo id
    }

    /*--------------------------------------------------------------------------
      Deactivate all fifos so that the dma does not start sending once the 
      descriptors are stored into the injection memory fifos.
      --------------------------------------------------------------------------*/
    for(i=0; i<8; i++){
	DMA_InjFifoSetDeactivateById( &inj_fifo_group, // fifo group structure
				      i);                   // fifo id
    }
    for(i=NumInjFifos; i<NumInjFifos+NumInjFifos_s; i++){
	DMA_InjFifoSetDeactivateById( &inj_fifo_group, // fifo group structure
				      i);                   // fifo id
    }

    GLOB_BARRIER;
    Vprintf("%s: Injection fifos initialized.\n", fname);
    fflush(stdout);
    GLOB_BARRIER;

    /*--------------------------------------------------------------------------
      Generate and store the descriptors that communicate the 
      "backward-projected-spinors" in the forward direction (+).
      Store the descriptors in the  injection fifos.
      This does not start the transfer because the inj fifos are deactivated.
      --------------------------------------------------------------------------*/


    /* Loop over the 4 directions */

    for (iPid=0; iPid<NUM_CORES; iPid++){
	if (Pid==iPid){

	    for(mu=0; mu<SFK_MAX_D; mu++) {
	      
		/* Set the hardware dimension hdim */
		hdim = dim_map[mu];
	      
		/* calculate number of bytes to send */
		blk_byte_size = wilson_p->comm_blklen[mu] * BLOCK * IFLOAT_IN_BYTES;

		/* initialize send address */
		send_off = comm_off[mu] + 8 * wilson_p->comm_offset[mu];

		/* initialize receive address */
		recv_off = comm_off[mu] + 0;

		/* Generate one describtor per consecutive data block */
		for(blk = 0; blk < wilson_p->comm_numblk[mu]; blk = blk+1) {
      
		    /* For torus */
		    if(hdim < 3) {
			rc = DMA_TorusDirectPutDescriptor(&inj_desc,
							  nn_x[hdim*2],
							  nn_y[hdim*2],
							  nn_z[hdim*2],
							  hints[hdim*2],
							  DMA_PACKET_VC_BN,
							  inj_c_grp,
							  inj_counter_id_fwd,
							  send_off,
							  rec_c_grp[Pid],
							  rec_counter_id_fwd[Pid],
							  recv_off,
							  blk_byte_size);

			if( rc == 0) {
			    // BGPDMA_VPRINT(" wfm_bgpdma: SUCCESS for DMA_TorusDirectPutDescriptor torus\n");
			    // printf("send_ad = 0x%08x, recv_ad = 0x%08x\n", send_ad, recv_ad);
			} else {
			    printf(" wfm_bgpdma: FAILURE for DMA_TorusDirectPutDescriptor torus rc=%d\n",rc);  
			    KillJob(rc);
			}

			rc = DMA_InjFifoInjectDescriptorById(&inj_fifo_group, (2*hdim), &inj_desc);
			if ( rc != 1) {
			    printf(" FAILURE for torus DMA_InjFifoInjectDescriptorById rc=%d\n",rc);
			    KillJob(rc);
			}
		    }
      
		    /* For local */
		    if(hdim == 3) {
			rc = DMA_LocalDirectPutDescriptor(
			    &inj_desc,
			    inj_c_grp,
			    inj_counter_id_fwd,
			    send_off,
			    rec_c_grp[(Pid+1)%NUM_CORES],
			    rec_counter_id_fwd[(Pid+1)%NUM_CORES],
			    recv_off,
			    blk_byte_size
			    );

			if( rc == 0) { 
			    // BGPDMA_VPRINT(" wfm_bgpdma: SUCCESS for DMA_TorusDirectPutDescriptor torus\n");
			    // printf("send_ad = 0x%08x, recv_ad = 0x%08x\n", send_ad, recv_ad);
			} else {
			    printf(" wfm_bgpdma: FAILURE for DMA_TorusDirectPutDescriptor torus rc=%d\n",rc);  
			    KillJob(rc);
			}
	
			rc = DMA_InjFifoInjectDescriptorById(&inj_fifo_group, (2*hdim), &inj_desc);
			if ( rc != 1) {
			    printf(" FAILURE for local DMA_InjFifoInjectDescriptorById rc=%d\n",rc);
			    KillJob(rc);
			}

		    }


		    /* calculate the next send address */
		    send_off = send_off + 8 * ( 
			wilson_p->comm_blklen[mu] * BLOCK 
			+ wilson_p->comm_stride[mu] 
			- 1);

		    /* calculate the next receive address */
		    recv_off = recv_off + 8 * (
			wilson_p->comm_blklen[mu] * BLOCK 
			+ wilson_p->comm_stride[mu] 
			- 1);
		}
	    }
	}
	NODE_BARRIER;
	_bgp_msync();
    }
    fflush(stdout);
    GLOB_BARRIER;


    /*--------------------------------------------------------------------------
      Generate and inject the descriptors that communicate the 
      "forward-projected-spinors" in the backward direction (-).
      Store the descriptors in the IFifo_bwd_torus[Pid], IFifo_bwd_local[Pid]
      injection fifos.
      This does not start the transfer because the inj fifos are deactivated.
      --------------------------------------------------------------------------*/

    for (iPid=0; iPid<NUM_CORES; iPid++){
	if (Pid==iPid){

	    /* Loop over the 4 directions */
	    // bytes_total = 0;
	    for(mu = 0; mu < SFK_MAX_D; mu++) {
     
		/* Set the hardware dimension hdim */
		hdim = dim_map[mu];

		/* calculate number of bytes to send */
		blk_byte_size = wilson_p->comm_blklen[mu] * BLOCK * IFLOAT_IN_BYTES;

		/* initialize send address */
		send_off = comm_off[mu];

		/* initialize receive address */
		recv_off = comm_off[mu] + 8 * wilson_p->comm_offset[mu];


		/* Generate one describtor per consecutive data block */
		for(blk = 0; blk < wilson_p->comm_numblk[mu]; blk = blk+1) {

		    /* For torus */
		    if(hdim < 3) {
			rc = DMA_TorusDirectPutDescriptor(&inj_desc,
							  nn_x[1+hdim*2],
							  nn_y[1+hdim*2],
							  nn_z[1+hdim*2],
							  hints[1+hdim*2],
							  DMA_PACKET_VC_BN,
							  inj_c_grp,
							  inj_counter_id_bwd,
							  send_off,
							  rec_c_grp[Pid],
							  rec_counter_id_bwd[Pid],
							  recv_off,
							  blk_byte_size);

			if( rc == 0) { 
			    // BGPDMA_VPRINT(" wfm_bgpdma: SUCCESS for DMA_TorusDirectPutDescriptor torus\n");
			    // printf("send_ad = 0x%08x, recv_ad = 0x%08x\n", send_ad, recv_ad);
			} else {
			    printf(" wfm_bgpdma: FAILURE for DMA_TorusDirectPutDescriptor torus rc=%d\n",rc);  
			    KillJob(rc);
			}

			rc = DMA_InjFifoInjectDescriptorById(&inj_fifo_group, ((2*hdim)+1), &inj_desc);
			if ( rc != 1) {
			    printf(" FAILURE for torus DMA_InjFifoInjectDescriptorById rc=%d\n",rc);
			    KillJob(rc);
			}

		    }
	      
 		    /* For local */
		    if(hdim == 3) {
			rc = DMA_LocalDirectPutDescriptor(
			    &inj_desc,
			    inj_c_grp,
			    inj_counter_id_bwd,
			    send_off,
			    rec_c_grp[(Pid+NUM_CORES-1)%NUM_CORES],
			    rec_counter_id_bwd[(Pid+NUM_CORES-1)%NUM_CORES],
			    recv_off,
			    blk_byte_size
			    );

			if( rc == 0) {
			    // BGPDMA_VPRINT(" wfm_bgpdma: SUCCESS for DMA_TorusDirectPutDescriptor torus\n");
			    // printf("send_ad = 0x%08x, recv_ad = 0x%08x\n", send_ad, recv_ad);
			} else {
			    printf(" wfm_bgpdma: FAILURE for DMA_TorusDirectPutDescriptor torus rc=%d\n",rc);  
			    KillJob(rc);
			}
	
			rc = DMA_InjFifoInjectDescriptorById(&inj_fifo_group, ((2*hdim)+1), &inj_desc);
			if ( rc != 1) {
			    printf(" FAILURE for local DMA_InjFifoInjectDescriptorById rc=%d\n",rc);
			    KillJob(rc);
			}

		    }
      
      
		    /* calculate the next send address */
		    send_off = send_off + 8 * ( 
			wilson_p->comm_blklen[mu] * BLOCK 
			+ wilson_p->comm_stride[mu] 
			- 1);

		    /* calculate the next receive address */
		    recv_off = recv_off + 8 * (
			wilson_p->comm_blklen[mu] * BLOCK 
			+ wilson_p->comm_stride[mu] 
			- 1);
		}

	    }

	    fflush(stdout);
	    GLOB_BARRIER;

	}
	NODE_BARRIER;
	_bgp_msync();
    }

    GLOB_BARRIER;
    Vprintf("%s: Double precision persistent comms set up, descriptors injected.\n", fname);
    fflush(stdout);
    GLOB_BARRIER;



    /*--------------------------------------------------------------------------
      Generate and store the descriptors that communicate the 
      "backward-projected-spinors" in the forward direction (+).
      Store the descriptors in the  injection fifos.
      This does not start the transfer because the inj fifos are deactivated.
      --------------------------------------------------------------------------*/


    /* Loop over the 4 directions */

    for (iPid=0; iPid<NUM_CORES; iPid++){
	if (Pid==iPid){

	    for(mu=0; mu<SFK_MAX_D; mu++) {
	      
		/* Set the hardware dimension hdim */
		hdim = dim_map[mu];
	      
		/* calculate number of bytes to send */
		blk_byte_size = wilson_p->comm_blklen[mu] * BLOCK * 4;

		/* initialize send address */
		send_off = comm_off[mu] + 4 * wilson_p->comm_offset[mu];

		/* initialize receive address */
		recv_off = comm_off[mu] + 0;

		/* Generate one describtor per consecutive data block */
		for(blk = 0; blk < wilson_p->comm_numblk[mu]; blk = blk+1) {
      
		    /* For torus */
		    if(hdim < 3) {
			rc = DMA_TorusDirectPutDescriptor(&inj_desc,
							  nn_x[hdim*2],
							  nn_y[hdim*2],
							  nn_z[hdim*2],
							  hints[hdim*2],
							  DMA_PACKET_VC_BN,
							  inj_c_grp,
							  inj_counter_id_fwd,
							  send_off,
							  rec_c_grp[Pid],
							  rec_counter_id_fwd[Pid],
							  recv_off,
							  blk_byte_size);

			if( rc == 0) {
			    // BGPDMA_VPRINT(" wfm_bgpdma: SUCCESS for DMA_TorusDirectPutDescriptor torus\n");
			    // printf("send_ad = 0x%08x, recv_ad = 0x%08x\n", send_ad, recv_ad);
			} else {
			    printf(" wfm_bgpdma: FAILURE for DMA_TorusDirectPutDescriptor torus rc=%d\n",rc);  
			    KillJob(rc);
			}

			rc = DMA_InjFifoInjectDescriptorById(&inj_fifo_group, NumInjFifos + (2*hdim), &inj_desc);
			if ( rc != 1) {
			    printf(" FAILURE for torus DMA_InjFifoInjectDescriptorById rc=%d\n",rc);
			    KillJob(rc);
			}
		    }
      
		    /* For local */
		    if(hdim == 3) {
			rc = DMA_LocalDirectPutDescriptor(
			    &inj_desc,
			    inj_c_grp,
			    inj_counter_id_fwd,
			    send_off,
			    rec_c_grp[(Pid+1)%NUM_CORES],
			    rec_counter_id_fwd[(Pid+1)%NUM_CORES],
			    recv_off,
			    blk_byte_size
			    );

			if( rc == 0) { 
			    // BGPDMA_VPRINT(" wfm_bgpdma: SUCCESS for DMA_TorusDirectPutDescriptor torus\n");
			    // printf("send_ad = 0x%08x, recv_ad = 0x%08x\n", send_ad, recv_ad);
			} else {
			    printf(" wfm_bgpdma: FAILURE for DMA_TorusDirectPutDescriptor torus rc=%d\n",rc);  
			    KillJob(rc);
			}
	
			rc = DMA_InjFifoInjectDescriptorById(&inj_fifo_group, NumInjFifos + (2*hdim), &inj_desc);
			if ( rc != 1) {
			    printf(" FAILURE for local DMA_InjFifoInjectDescriptorById rc=%d\n",rc);
			    KillJob(rc);
			}

		    }


		    /* calculate the next send address */
		    send_off = send_off + 4 * ( 
			wilson_p->comm_blklen[mu] * BLOCK 
			+ wilson_p->comm_stride[mu] 
			- 1);

		    /* calculate the next receive address */
		    recv_off = recv_off + 4 * (
			wilson_p->comm_blklen[mu] * BLOCK 
			+ wilson_p->comm_stride[mu] 
			- 1);
		}
	    }
	}
	NODE_BARRIER;
	_bgp_msync();
    }
    fflush(stdout);
    GLOB_BARRIER;


    /*--------------------------------------------------------------------------
      Generate and inject the descriptors that communicate the 
      "forward-projected-spinors" in the backward direction (-).
      Store the descriptors in the IFifo_bwd_torus[Pid], IFifo_bwd_local[Pid]
      injection fifos.
      This does not start the transfer because the inj fifos are deactivated.
      --------------------------------------------------------------------------*/

    for (iPid=0; iPid<NUM_CORES; iPid++){
	if (Pid==iPid){

	    /* Loop over the 4 directions */
	    // bytes_total = 0;
	    for(mu = 0; mu < SFK_MAX_D; mu++) {
     
		/* Set the hardware dimension hdim */
		hdim = dim_map[mu];

		/* calculate number of bytes to send */
		blk_byte_size = wilson_p->comm_blklen[mu] * BLOCK * 4;

		/* initialize send address */
		send_off = comm_off[mu];

		/* initialize receive address */
		recv_off = comm_off[mu] + 4 * wilson_p->comm_offset[mu];


		/* Generate one describtor per consecutive data block */
		for(blk = 0; blk < wilson_p->comm_numblk[mu]; blk = blk+1) {

		    /* For torus */
		    if(hdim < 3) {
			rc = DMA_TorusDirectPutDescriptor(&inj_desc,
							  nn_x[1+hdim*2],
							  nn_y[1+hdim*2],
							  nn_z[1+hdim*2],
							  hints[1+hdim*2],
							  DMA_PACKET_VC_BN,
							  inj_c_grp,
							  inj_counter_id_bwd,
							  send_off,
							  rec_c_grp[Pid],
							  rec_counter_id_bwd[Pid],
							  recv_off,
							  blk_byte_size);

			if( rc == 0) { 
			    // BGPDMA_VPRINT(" wfm_bgpdma: SUCCESS for DMA_TorusDirectPutDescriptor torus\n");
			    // printf("send_ad = 0x%08x, recv_ad = 0x%08x\n", send_ad, recv_ad);
			} else {
			    printf(" wfm_bgpdma: FAILURE for DMA_TorusDirectPutDescriptor torus rc=%d\n",rc);  
			    KillJob(rc);
			}

			rc = DMA_InjFifoInjectDescriptorById(&inj_fifo_group, NumInjFifos + ((2*hdim)+1), &inj_desc);
			if ( rc != 1) {
			    printf(" FAILURE for torus DMA_InjFifoInjectDescriptorById rc=%d\n",rc);
			    KillJob(rc);
			}

		    }
	      
 		    /* For local */
		    if(hdim == 3) {
			rc = DMA_LocalDirectPutDescriptor(
			    &inj_desc,
			    inj_c_grp,
			    inj_counter_id_bwd,
			    send_off,
			    rec_c_grp[(Pid+NUM_CORES-1)%NUM_CORES],
			    rec_counter_id_bwd[(Pid+NUM_CORES-1)%NUM_CORES],
			    recv_off,
			    blk_byte_size
			    );

			if( rc == 0) {
			    // BGPDMA_VPRINT(" wfm_bgpdma: SUCCESS for DMA_TorusDirectPutDescriptor torus\n");
			    // printf("send_ad = 0x%08x, recv_ad = 0x%08x\n", send_ad, recv_ad);
			} else {
			    printf(" wfm_bgpdma: FAILURE for DMA_TorusDirectPutDescriptor torus rc=%d\n",rc);  
			    KillJob(rc);
			}
	
			rc = DMA_InjFifoInjectDescriptorById(&inj_fifo_group, NumInjFifos + ((2*hdim)+1), &inj_desc);
			if ( rc != 1) {
			    printf(" FAILURE for local DMA_InjFifoInjectDescriptorById rc=%d\n",rc);
			    KillJob(rc);
			}

		    }
      
      
		    /* calculate the next send address */
		    send_off = send_off + 4 * ( 
			wilson_p->comm_blklen[mu] * BLOCK 
			+ wilson_p->comm_stride[mu] 
			- 1);

		    /* calculate the next receive address */
		    recv_off = recv_off + 4 * (
			wilson_p->comm_blklen[mu] * BLOCK 
			+ wilson_p->comm_stride[mu] 
			- 1);
		}

	    }

	    fflush(stdout);
	    GLOB_BARRIER;

	}
	NODE_BARRIER;
	_bgp_msync();
    }

    GLOB_BARRIER;
    Vprintf("%s: Single precision persistent comms set up, descriptors injected.\n", fname);
    fflush(stdout);
    GLOB_BARRIER;

    for (i=0; i<8; i++) {
	tail = (uint32_t) DMA_FifoGetTail(&(inj_fifo_group.fifos[i].dma_fifo));
	DMA_InjFifoSetHead(&(inj_fifo_group.fifos[i]), (void *) tail);
    }
    /* single prec. fifos */
    for (i = NumInjFifos; i< NumInjFifos + NumInjFifos_s; i++) {
	tail = (uint32_t) DMA_FifoGetTail(&(inj_fifo_group.fifos[i].dma_fifo));
	DMA_InjFifoSetHead(&(inj_fifo_group.fifos[i]), (void *) tail);
    }
    
    for(i=0; i<2; i++){
	DMA_CounterSetEnableById(&inj_counter_group, i);
	DMA_CounterSetEnableById(&rec_counter_group, i);
    }
    
    for(i=0; i<8; i++){
	DMA_InjFifoSetActivateById(&inj_fifo_group, // fifo group structure
				   i);                // fifo id
    }
    /* single prec. fifos*/
    for (i = NumInjFifos; i< NumInjFifos + NumInjFifos_s; i++) {
	DMA_InjFifoSetActivateById(&inj_fifo_group, // fifo group structure
				   i);                // fifo id
    }

    GLOB_BARRIER;
    Vprintf("%s: Fifo heads set, counters enabled and fifos activated.\n", fname);
    fflush(stdout);
    GLOB_BARRIER;

    /* initialize the shared memory window. */

    switch( Pid )
    {
	case 0:
	    // Map the shared memory using the SHM_FILE handle
	    fd = shm_open( SHM_FILE, O_CREAT | O_RDWR, 0600 );

	    if ( fd == -1 )
	    {
		printf(" wfm_bgpdma: Initial shm_open failed with errno=%d\n", errno);
		KillJob(ERR_ALLOC_SHMEM_FAIL);
	    }

	    // Set up the size of the shared memory
	    rc = ftruncate( fd, MAX_SHARED_SIZE );

	    if ( rc == -1 )
	    {
		printf(" wfm_bgpdma: ftruncate failed with errno=%d\n", errno);
		KillJob(ERR_ALLOC_SHMEM_FAIL);

	    }

	    // Get a pointer to the shared memory
	    bgp_wfm_shmptr = mmap( NULL, MAX_SHARED_SIZE, PROT_READ | PROT_WRITE, 
			   MAP_SHARED, fd, 0);

	    //printf(" wfm_bgpdma: Pid 0: bgp_wfm_shmptr=0x%08x\n", (unsigned int) bgp_wfm_shmptr);
	
	    CHIP_BARRIER;
	    break;

	case 1:

	    CHIP_BARRIER;

	    fd = shm_open( SHM_FILE, O_RDWR, 0600 );
 
	    bgp_wfm_shmptr = mmap( NULL, MAX_SHARED_SIZE, PROT_READ | PROT_WRITE, 
			    MAP_SHARED, fd, 0);
	    //printf(" wfm_bgpdma: Pid 1: bgp_wfm_shmptr=0x%08x\n", (unsigned int) bgp_wfm_shmptr);

	    break;

	case 2:

	    CHIP_BARRIER;

	    fd = shm_open( SHM_FILE, O_RDWR, 0600 );

	    bgp_wfm_shmptr = mmap( NULL, MAX_SHARED_SIZE, PROT_READ | PROT_WRITE, 
			    MAP_SHARED, fd, 0);
	    //printf(" wfm_bgpdma: Pid 2: bgp_wfm_shmptr=0x%08x\n", (unsigned int) bgp_wfm_shmptr);

	    break;

	case 3:

	    CHIP_BARRIER;

	    fd = shm_open( SHM_FILE, O_RDWR, 0600 );

	    bgp_wfm_shmptr = mmap( NULL, MAX_SHARED_SIZE, PROT_READ | PROT_WRITE, 
			    MAP_SHARED, fd, 0);
	    //printf(" wfm_bgpdma: Pid 3: bgp_wfm_shmptr=0x%08x\n", (unsigned int) bgp_wfm_shmptr);

	    break;

	default:
	    printf(" wfm_bgpdma: ERROR: Pid not in [0123].\n");
	    KillJob(ERR_INTERNAL_INCONSISTENCY);

	    break;
    }

    CHIP_BARRIER;

    if ( bgp_wfm_shmptr == MAP_FAILED )
    {
	printf(" wfm_bgpdma: Initial mmap failed with errno=%d\n", errno);
	KillJob(ERR_ALLOC_SHMEM_FAIL);
    }    

    GLOB_BARRIER;
    Vprintf("%s: Shared memory allocated. Pointer: %x\n", fname, (unsigned int) bgp_wfm_shmptr);
    fflush(stdout);
    GLOB_BARRIER;

    /* set the initialized flag */
    bgp_wfm_initialized=1;
    GLOB_BARRIER;

    Vprintf("%s: Initialization complete.\n", fname);
    fflush(stdout);
    GLOB_BARRIER;

    return;
}



/*----------------------------------------------------------------------------
  In this routine the BG/P dma is started in order to service the descriptors 
  for the backward-projected-spinor transfer along the forward direction.

  Calling this subroutine sets the fwd memory fifo head to 
  WFM_BGPDMA_FWD_FIFO_HEAD which starts the dma.

  -

  I will have to use correct virtual addresses to run in cna mode. Therefore
  the prototype was changed to have the wilson_p as parameter and I will use
  ab[0] as smallest and &ab[0][comm_off_stride] as largest va.
  sfk 07/2007
  ----------------------------------------------------------------------------*/
void wfm_bgpdma_fwd_start(Wilson *wilson_p)
{
    char *fname = "wfm_bgpdma_fwd_start";
    int Pid = GET_PID;
    int rc;
    extern uint32_t comm_off_stride;
    
    // setting reception counters
    rc = DMA_CounterSetValueBaseMaxById(&rec_counter_group, 
					rec_counter_id_fwd[Pid],
					comm_bytes, 
					wilson_p->ab[0], 
					&wilson_p->ab[0][comm_off_stride]);
    if ( rc != 0) {
        printf(" %s: FAILURE in  DMA_CounterSetValueBaseMaxById rc=%d, "
	       "errno=%i\n",fname, rc, errno);
        KillJob(rc);
    }

    // setting injection counters
    rc = DMA_CounterSetValueBaseById(&inj_counter_group, 
				     inj_counter_id_fwd,
				     comm_bytes,
				     wilson_p->ab[0]);
    if ( rc != 0) {
        printf(" %s: FAILURE in  DMA_CounterSetValueBaseById rc=%d, "
	       "errno=%i\n",fname, rc, errno);
        KillJob(rc);
    }

    GLOB_BARRIER;

    if (SFK_MAX_D>0) DMA_InjFifoSetHead(&(inj_fifo_group.fifos[0]), (void *) IFifo[0]);
    if (SFK_MAX_D>1) DMA_InjFifoSetHead(&(inj_fifo_group.fifos[2]), (void *) IFifo[2]);
    if (SFK_MAX_D>2) DMA_InjFifoSetHead(&(inj_fifo_group.fifos[4]), (void *) IFifo[4]);
    if (SFK_MAX_D>3) DMA_InjFifoSetHead(&(inj_fifo_group.fifos[6]), (void *) IFifo[6]);

    return;

}


void wfm_bgpdma_fwd_start_single(Wilson *wilson_p)
{
    char *fname = "wfm_bgpdma_fwd_start_single";
    int Pid = GET_PID;
    int rc;
    extern uint32_t comm_off_stride;
    
    // setting reception counters
    rc = DMA_CounterSetValueBaseMaxById(&rec_counter_group, 
					rec_counter_id_fwd[Pid],
					comm_bytes_s, 
					wilson_p->ab[0], 
					&wilson_p->ab[0][comm_off_stride]);
    if ( rc != 0) {
        printf(" %s: FAILURE in  DMA_CounterSetValueBaseMaxById rc=%d, "
	       "errno=%i\n",fname, rc, errno);
        KillJob(rc);
    }

    // setting injection counters
    rc = DMA_CounterSetValueBaseById(&inj_counter_group, 
				     inj_counter_id_fwd,
				     comm_bytes_s,
				     wilson_p->ab[0]);
    if ( rc != 0) {
        printf(" %s: FAILURE in  DMA_CounterSetValueBaseById rc=%d, "
	       "errno=%i\n",fname, rc, errno);
        KillJob(rc);
    }

    GLOB_BARRIER;

    if (SFK_MAX_D>0) DMA_InjFifoSetHead(&(inj_fifo_group.fifos[NumInjFifos + 0]), (void *) IFifo_s[0]);
    if (SFK_MAX_D>1) DMA_InjFifoSetHead(&(inj_fifo_group.fifos[NumInjFifos + 2]), (void *) IFifo_s[2]);
    if (SFK_MAX_D>2) DMA_InjFifoSetHead(&(inj_fifo_group.fifos[NumInjFifos + 4]), (void *) IFifo_s[4]);
    if (SFK_MAX_D>3) DMA_InjFifoSetHead(&(inj_fifo_group.fifos[NumInjFifos + 6]), (void *) IFifo_s[6]);

    return;
}



/*----------------------------------------------------------------------------
  In this routine the BG/P dma is started in order to service the descriptors 
  for the forward-projected-spinor transfer along the backward direction.


  Calling this subroutine sets the bwd memory fifo head to 
  WFM_BGPDMA_BWD_FIFO_HEAD which starts the dma.
  -

  I will have to use correct virtual addresses to run in cna mode. Therefore
  the prototype was changed to have the wilson_p as parameter and I will use
  af[0] as smallest and &af[0][comm_off_stride] as largest va.
  sfk 07/2007
  ----------------------------------------------------------------------------*/
void wfm_bgpdma_bwd_start(Wilson *wilson_p)
{
    char *fname = "wfm_bgpdma_bwd_start";
    int Pid = GET_PID;
    extern uint32_t comm_off_stride;
    int rc = 0;

    rc = DMA_CounterSetValueBaseMaxById(&rec_counter_group, 
					rec_counter_id_bwd[Pid],
					comm_bytes, 
					wilson_p->af[0], 
					&wilson_p->af[0][comm_off_stride]);

    if ( rc != 0) {
        printf(" %s: FAILURE in  DMA_CounterSetValueBaseMaxById rc=%d, "
	       "errno=%i\n", fname, rc, errno);
        KillJob(rc);
    }

    rc = DMA_CounterSetValueBaseById(&inj_counter_group, 
				     inj_counter_id_bwd,
				     comm_bytes,
				     wilson_p->af[0]);

    if ( rc != 0) {
	printf(" %s: FAILURE in  DMA_CounterSetValueBaseById rc=%d, "
	       "errno=%i\n", fname, rc, errno);
	KillJob(rc);
    }

    GLOB_BARRIER;

    if (SFK_MAX_D>0) DMA_InjFifoSetHead(&(inj_fifo_group.fifos[1]), (void *) IFifo[1]);
    if (SFK_MAX_D>1) DMA_InjFifoSetHead(&(inj_fifo_group.fifos[3]), (void *) IFifo[3]);
    if (SFK_MAX_D>2) DMA_InjFifoSetHead(&(inj_fifo_group.fifos[5]), (void *) IFifo[5]);
    if (SFK_MAX_D>3) DMA_InjFifoSetHead(&(inj_fifo_group.fifos[7]), (void *) IFifo[7]);
    
    return;
}

void wfm_bgpdma_bwd_start_single(Wilson *wilson_p)
{
    char *fname = "wfm_bgpdma_bwd_start_single";
    int Pid = GET_PID;
    extern uint32_t comm_off_stride;
    int rc = 0;

    rc = DMA_CounterSetValueBaseMaxById(&rec_counter_group, 
					rec_counter_id_bwd[Pid],
					comm_bytes_s, 
					wilson_p->af[0], 
					&wilson_p->af[0][comm_off_stride]);

    if ( rc != 0) {
        printf(" %s: FAILURE in  DMA_CounterSetValueBaseMaxById rc=%d, "
	       "errno=%i\n", fname, rc, errno);
        KillJob(rc);
    }

    rc = DMA_CounterSetValueBaseById(&inj_counter_group, 
				     inj_counter_id_bwd,
				     comm_bytes_s,
				     wilson_p->af[0]);

    if ( rc != 0) {
	printf(" %s: FAILURE in  DMA_CounterSetValueBaseById rc=%d, "
	       "errno=%i\n", fname, rc, errno);
	KillJob(rc);
    }

    GLOB_BARRIER;

    if (SFK_MAX_D>0) DMA_InjFifoSetHead(&(inj_fifo_group.fifos[NumInjFifos + 1]), (void *) IFifo_s[1]);
    if (SFK_MAX_D>1) DMA_InjFifoSetHead(&(inj_fifo_group.fifos[NumInjFifos + 3]), (void *) IFifo_s[3]);
    if (SFK_MAX_D>2) DMA_InjFifoSetHead(&(inj_fifo_group.fifos[NumInjFifos + 5]), (void *) IFifo_s[5]);
    if (SFK_MAX_D>3) DMA_InjFifoSetHead(&(inj_fifo_group.fifos[NumInjFifos + 7]), (void *) IFifo_s[7]);
    
    return;
}


/*----------------------------------------------------------------------------
  This routine waits until the counter for the forward dma communication, 
  that transfers the backward-projected-spinors, (with counter id  
  WFM_BGPDMA_FWD_CNT_ID) reaches zero. This indicates that the transfer is 
  done and then the routine returns. 

  This routine is "blocking" until the dma completes the communications. 
  ----------------------------------------------------------------------------*/
void wfm_bgpdma_fwd_wait(void)
{
    int ctr_rec;
    int ctr_inj;
    int c=0;
    int Pid = GET_PID;
    char *fname = "wfm_bgpdma_fwd_wait";

    while (1) {
        ctr_rec = DMA_CounterGetValueById(&rec_counter_group, rec_counter_id_fwd[Pid]);
        if (ctr_rec == 0) break;
	if ((c++ >= 100) ) {
	    printf(" %s: pid %i: fwd: ctr_rec = %d\n", fname, Pid, ctr_rec);
	    ctr_inj = DMA_CounterGetValueById(&inj_counter_group, inj_counter_id_fwd);
	    printf(" %s: pid %i: fwd: ctr_inj = %d\n", fname, Pid, ctr_inj);
	    fflush(stdout);
	    KillJob(ERR_COMM_TIMEOUT);
	}
	ppc450_Delay(50);
    }
    
    return;
}


/*----------------------------------------------------------------------------
  This routine waits until the counter for the backward dma communication, 
  that transfers the forward-projected-spinors, (with counter id  
  WFM_BGPDMA_BWD_CNT_ID) reaches zero. This indicates that the transfer is 
  done and then the routine returns. 

  This routine is "blocking" until the dma completes the communications. 
  ----------------------------------------------------------------------------*/
void wfm_bgpdma_bwd_wait(void)
{
    int ctr_rec;
    int ctr_inj;
    int c=0;
    int Pid = GET_PID;
    char *fname = "wfm_bgpdma_bwd_wait";

    while (1) {
        ctr_rec = DMA_CounterGetValueById(&rec_counter_group, rec_counter_id_bwd[Pid]);
        if (ctr_rec == 0) break;
	if ((c++ >= 100) ) {
	    printf(" %s: pid %i: bwd: ctr_rec = %d\n", fname, Pid, ctr_rec);
	    ctr_inj = DMA_CounterGetValueById(&inj_counter_group, inj_counter_id_bwd);
	    printf(" %s: pid %i: bwd: ctr_inj = %d\n", fname, Pid, ctr_inj);
	    fflush(stdout);
	    KillJob(ERR_COMM_TIMEOUT);
	}
	ppc450_Delay(50);
    }

    return;
}


/*----------------------------------------------------------------------------
  Generic comm routine to send data to another node.
  Limitations:
              - To reduce the number of counters used, the FWD counters 
	        will be used here as well. -> blocking comms
		                                            05.07.07 sfk
  ----------------------------------------------------------------------------*/

void wfm_bgpdma_comm_start(void *recv_a, void *send_a, int num_bytes, int dir)
{
    char *fname = "wfm_bgpdma_comm";
    int Pid = GET_PID;
    /* int blk; */
    /* int numblk=1; */
    int blk_byte_size;
    int hdim;
    uint32_t send_off, recv_off;
    DMA_InjDescriptor_t inj_desc; /* Descriptor for prep&immediate inj into fifo */
    int rc;
    int mu = dir/2;
    int pn = dir%2;
    char *recv_addr = (char *) recv_a;
    char *send_addr = (char *) send_a;
 

    printf ("%s: mu = %i, pn = %i\n", fname, mu, pn);
    fflush(stdout);

    // setting reception counters
    if (pn){
	rc = DMA_CounterSetValueBaseMaxById(&rec_counter_group, 
					    rec_counter_id_bwd[Pid],
					    num_bytes, 
					    recv_addr, 
					    &recv_addr[num_bytes]);
    }else{
	rc = DMA_CounterSetValueBaseMaxById(&rec_counter_group, 
					    rec_counter_id_fwd[Pid],
					    num_bytes, 
					    recv_addr, 
					    &recv_addr[num_bytes]);
    }
    
    if ( rc != 0) {
	printf(" %s: FAILURE in  DMA_CounterSetValueBaseMaxById rc=%d, "
	       "errno=%i\n",fname, rc, errno);
	KillJob(rc);
    } 
    
    // setting injection counters
    if (pn){
	rc = DMA_CounterSetValueBaseById(&inj_counter_group, 
					 inj_counter_id_bwd,
					 num_bytes,
					 send_addr);
    }else{
	rc = DMA_CounterSetValueBaseById(&inj_counter_group, 
					 inj_counter_id_fwd,
					 num_bytes,
					 send_addr);
    }	
    if ( rc != 0) {
	printf(" %s: FAILURE in  DMA_CounterSetValueBaseById rc=%d, "
	       "errno=%i\n",fname, rc, errno);
	KillJob(rc);
    }

    GLOB_BARRIER;

    /* initialize send address */
    send_off = 0;
    
    /* initialize receive address */
    recv_off = 0;
    
    /* Set the hardware dimension hdim */
    hdim = bgp_wfm_dim_map[mu];
    
    /* calculate number of bytes to send */
    blk_byte_size = num_bytes;


    /* For torus */
    if (hdim<3){
	if (pn){
	    rc = DMA_TorusDirectPutDescriptor(&inj_desc,
					      nn_x[hdim*2+1], //odd: neg dir
					      nn_y[hdim*2+1],
					      nn_z[hdim*2+1],
					      DMA_PACKET_HINT_NULL,
					      DMA_PACKET_VC_BN,
					      inj_c_grp,
					      inj_counter_id_bwd,
					      send_off,
					      rec_c_grp[Pid],
					      rec_counter_id_bwd[Pid],
					      recv_off,
					      blk_byte_size);
	}else{
	    rc = DMA_TorusDirectPutDescriptor(&inj_desc,
					      nn_x[hdim*2],
					      nn_y[hdim*2],
					      nn_z[hdim*2],
					      DMA_PACKET_HINT_NULL,
					      DMA_PACKET_VC_BN,
					      inj_c_grp,
					      inj_counter_id_fwd,
					      send_off,
					      rec_c_grp[Pid],
					      rec_counter_id_fwd[Pid],
					      recv_off,
					      blk_byte_size);
	}
	if(rc){
	    printf(" wfm_bgpdma: FAILURE for DMA_TorusDirectPutDescriptor torus rc=%d\n",rc);  
	    KillJob(rc);
	}
	rc = DMA_InjFifoInjectDescriptorById(&inj_fifo_group, AuxInjFifoTorus, &inj_desc);
	if ( rc != 1) {
	    printf(" FAILURE for torus DMA_InjFifoInjectDescriptorById rc=%d\n",rc);
	    KillJob(rc);
	}

    }
    /* For local */
    if(hdim == 3) {
	if (pn){
	    rc = DMA_LocalDirectPutDescriptor(&inj_desc,
					      inj_c_grp,
					      inj_counter_id_bwd,
					      send_off,
					      rec_c_grp[(Pid+NUM_CORES-1)%NUM_CORES],
					      rec_counter_id_bwd[(Pid+NUM_CORES-1)%NUM_CORES],
					      recv_off,
					      blk_byte_size);
	}else{
	    rc = DMA_LocalDirectPutDescriptor(&inj_desc,
					      inj_c_grp,
					      inj_counter_id_fwd,
					      send_off,
					      rec_c_grp[(Pid+1)%NUM_CORES],
					      rec_counter_id_fwd[(Pid+1)%NUM_CORES],
					      recv_off,
					      blk_byte_size);
	}
	    
	if(rc){
	    printf(" wfm_bgpdma: FAILURE for DMA_TorusDirectPutDescriptor torus rc=%d\n",rc);  
	    KillJob(rc);
	}

	rc = DMA_InjFifoInjectDescriptorById(&inj_fifo_group, AuxInjFifoLocal, &inj_desc);
	if ( rc != 1) {
	    printf(" FAILURE for torus DMA_InjFifoInjectDescriptorById rc=%d\n",rc);
	    KillJob(rc);
	}
    }

    return;
}

void wfm_bgpdma_comm_wait(int dir)
{
    int pn = dir%2;
    if (pn){
	wfm_bgpdma_bwd_wait();
    }else{
	wfm_bgpdma_fwd_wait();
    }
    
    return;
}
