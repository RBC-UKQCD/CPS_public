/*
 *
 *  Extern "C" routine to augment the set of optimised implementations
 *  for the Columbia physics codes.
 *
 *  Calls PAB's assembler routines to give an implementation
 *  of the wilson dslash
 *
 *  Replicates funcionality of wfm_dslash.asm in qcdsp optimised code
 *
 */
#include <util/wilson.h>
#include "prec.h"
#include <string.h>

extern "C" {

void QCDOC_ChDecom_m      (void *psi,void *len,void *tab_bwd);
void QCDOC_ChDecom_p      (void *psi,void *len,void *tab_bwd);
void QCDOC_ChDecom_m_hsu3 (void *psi,void *Ucb ,void *len,void *tab_fwd);
void QCDOC_ChDecom_p_hsu3 (void *psi,void *Ucb ,void *len,void *tab_fwd);

void QCDOC_ChRecon_m_su3  (void *chi,void *Un  ,void *chib,void *len);
void QCDOC_ChRecon_p_su3  (void *chi,void *Un  ,void *chib,void *len);
void QCDOC_ChRecon_m_add  (void *chi,void * chif,void *len);
void QCDOC_ChRecon_p_add  (void *chi,void * chif,void *len);

void wfm_abi_save(void);
void wfm_abi_restore(void);


void wilson_dslash(Float *chi, 
		   Float *u, 
		   Float *psi, 
		   int cb,
		   int dag,
		   Wilson *wilson_p)
{

  /*
   * To a first approximation, we simply
   * remap the arguments into a form acceptable
   * to the assembler, then call it
   */

  fpoint *gauge[2];
  fpoint *gauge_par;
  fpoint *gauge_notpar;
  fpoint *ChiF;
  fpoint *ChiB;
  asmint *tabledec;
  asmint *tabledecsu3;
  asmint length;
  int mu;

  wfm_abi_save();

  length = wilson_p->vol[0];

  /*2 spinor args - just refer to wilson_p*/
  ChiF = (fpoint *)wilson_p->af;
  ChiB = (fpoint *)wilson_p->ab;

  /*Gauge args*/
  gauge[0] = (fpoint *)u;
  gauge[1] = (fpoint *)u + GAUGE_SIZE*wilson_p->vol[0];
  gauge_par = gauge[cb] ;
  gauge_notpar = gauge[1-cb];
  
  /*
   * Decompose and shift back from x+mu to x 
   */
  tabledec = (asmint *)wilson_p->shift_table[cb][1];
  tabledecsu3 = (asmint *)wilson_p->shift_table[cb][0];


  if (dag) { 

    QCDOC_ChDecom_p(psi,
                    &length,
                    tabledec
		    );

  }else{

    QCDOC_ChDecom_m(psi,
                    &length,
                    tabledec
		    );

  }

  /*
   * Hit the PEC streams early
   */

  /*
   * Send backwards
   */
  wfm_comm_backward_start(wilson_p);

  /*
   * Now do fused decompose Udagger(x) Chi(x) and shift forwards
   */

  if ( dag ) {
    QCDOC_ChDecom_m_hsu3(psi,
                         gauge_par,
                         &length,
                         tabledecsu3
                        );

  }  else {
    QCDOC_ChDecom_p_hsu3(psi,
                         gauge_par,
                         &length,
                         tabledecsu3
			 );
  }


  /*
   * Send forwards
   */
  wfm_comm_forward_start(wilson_p);

  /*
   * Receive backwards
   */
  wfm_comm_backward_complete(wilson_p);


  /*
   * Scatter the backward face from buffer
   */
  wfm_scatter_face(wilson_p,1,cb);

  /*
   * U multiply and Reconstruct the received 2 spinors
   */

  if ( dag ) {

    QCDOC_ChRecon_p_su3( chi,
                         gauge_notpar,
                         ChiB,
			 &length
                        );

  } else {

    QCDOC_ChRecon_m_su3( chi,
                         gauge_notpar,
                         ChiB,
                         &length
                        );
  }

  /*Receive forwards*/
  wfm_comm_forward_complete(wilson_p);
  /*Scatter the forward face from buffer*/
  wfm_scatter_face(wilson_p,0,cb);


  /*
   * Reconstruct the final 2 spinors
   */

  if ( dag ) { 

    QCDOC_ChRecon_m_add(chi,
			ChiF,
			&length
			);
  } else {

    QCDOC_ChRecon_p_add(chi,
			ChiF,
                        &length
			);
  }

  /*
   * Take a bow, exit stage right, and come back for 
   * an encore since the matrix is preconditioned...
   */
  wfm_abi_restore();
  return;
}

}







