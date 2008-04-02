/****************************************************************************
This header file has all necessary declarations and external variable 
definitions needed by the BlueGene/P dma communications for the
Wilson fermion kernel.                       
                                            P.Vranas  10/2006
Corrected, parallelized in the time direction and converted to VN mode
and userspace, fixed timing and added communication routines for single 
precision and get{Plus,Minus}Data routines. 
                                            S. Krieg  07/2007
*****************************************************************************/
#ifndef _WFM_BGPDMA_H_
#define _WFM_BGPDMA_H_

#include "wilson.h"
/*----------------------------------------------------------------------------
 Globals
----------------------------------------------------------------------------*/
#ifndef _WFM_BGPDMA_C_
extern int bgp_wfm_initialized;
extern int bgp_wfm_hw_my_rank;
extern int bgp_wfm_hw_my_coords[4];
extern int bgp_wfm_hw_num_nodes[4];
extern int bgp_wfm_dim_map[4];
extern char *bgp_wfm_shmptr;
#endif

/*----------------------------------------------------------------------------
 Function declarations
----------------------------------------------------------------------------*/

/* 
Initializes the bgpdma for the Wilson comms
*/
void wfm_bgpdma_init(Wilson *wilson_p);

/* 
Start the forward, backward dma comms
*/
void wfm_bgpdma_fwd_start(Wilson *wilson_p);
void wfm_bgpdma_bwd_start(Wilson *wilson_p);

void wfm_bgpdma_fwd_start_single(Wilson *wilson_p);
void wfm_bgpdma_bwd_start_single(Wilson *wilson_p);

/* 
Wait for  the forward, backward dma comms to complete
*/
void wfm_bgpdma_fwd_wait(void);
void wfm_bgpdma_bwd_wait(void);

/* 
Generic comms along direction dir 
WARINING: CANNOT SEND MORE THAN ONE DIR AT THE TIME
          because counters are set, not incremented.
	  Subject to change.
*/
void wfm_bgpdma_comm_start(void *recv_a, void *send_a, int num_bytes, 
		     int dir);
void wfm_bgpdma_comm_wait(int dir);
#endif
