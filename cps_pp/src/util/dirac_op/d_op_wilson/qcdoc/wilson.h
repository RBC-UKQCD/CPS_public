/****************************************************************************/
/* 10/16/97                                                                 */
/*                                                                          */
/* wilsonpab.h                                                              */
/*                                                                          */
/* C header file for the fermion wilson library wilson.lib                  */
/* 10/1/2002 PAB... modify for QCDOC optimised matrix multiply              */
/*                                                                          */
/****************************************************************************/

#ifndef INCLUDED_WILSON_PAB_H
#define INCLUDED_WILSON_PAB_H

#include <qcdocos.h> /*Needed for SCUDirArg*/
typedef double Float ;
typedef double IFloat ;

/*--------------------------------------------------------------------------*/
/* Definitions                                                              */
/*--------------------------------------------------------------------------*/
#define ND                 4      /* Space time dimension                   */
#define SPINOR_SIZE        24     /* Dirac spinor size                      */
#define HALF_SPINOR_SIZE   12     /* Half of the Dirac spinor size          */
#define BLOCK   HALF_SPINOR_SIZE  /* Size of spin projected 2 comp. spinor  */
#define PAD_HALF_SPINOR_SIZE 16   /* We pad to a divisor of the PEC reg size*/
                                  /* to avoid compulsory reg misses         */
#define COLUMN_SPINOR_SIZE  6     /* Size of one spinor component           */
#define GAUGE_SIZE         72     /* Gauge field size per site              */

/*--------------------------------------------------------------------------*/
/* Type definitions                                                         */
/*--------------------------------------------------------------------------*/
/* The Wilson structure typedef */
typedef struct{
   unsigned long *ptr;
   unsigned long *shift_table[2][2];/*pointer to an array with addressing offsets */
   unsigned long *face_table[2][2][ND];/*pointer to an array with addressing offsets */

   int   vol[2];
   IFloat *spinor_tmp;     /* temp spinor needed by mdagm                 */
   IFloat *af;     /* point. array to 4 interleaved fwd proj half spinors  */
   IFloat *ab;     /* point. array to 4 interleaved bwd proj half spinors  */
   IFloat *send_f[ND];
   IFloat *recv_f[ND];
   IFloat *send_b[ND];
   IFloat *recv_b[ND];

   int local_latt[ND];
   int nbound[4];
   int allbound;
   int local_comm[4];
   SCUDirArgMulti *comm_f;
   SCUDirArgMulti *comm_b;
} Wilson;



/*--------------------------------------------------------------------------*/
/* Function declarations                                                    */
/*--------------------------------------------------------------------------*/
#ifdef __cplusplus
extern "C" {
#endif

/*==========================================================================*/
/* This routine performs all initializations needed before wilson library   */
/* funcs are called. It sets the addressing related arrays and reserves     */
/* memory for the needed temporary buffers. It only needs to be called only */
/* once at the begining of the program (or after a wilson_end call)         */
/* before any number of calls to wilson functions are made.                 */
/*==========================================================================*/
void wilson_init(Wilson *wilson_p);     /* pointer to Wilson struct. */


/*==========================================================================*/
/* This routine frees any memory reserved by wilson_init                    */
/*==========================================================================*/
void wilson_end(Wilson *wilson_p);     /* pointer to Wilson struct. */


/*Comms routines*/

void wfm_scatter_face(Wilson *wilson_p, int pm, int cb);
void wfm_comm_init(Wilson *wp);
void wfm_comm_forward_start(Wilson *wp);
void wfm_comm_backward_start(Wilson *wp);
void wfm_comm_forward_complete(Wilson *wp);
void wfm_comm_backward_complete(Wilson *wp);

/*==========================================================================*/
/* The Wilson fermion matrix * vector routines                              */
/* They are external "C"                                                    */
/*==========================================================================*/

void wilson_mdagm(IFloat  *chi,        /* chi = MdagM(u) psi          */
		  IFloat  *u,          /* Gauge field                 */
		  IFloat  *psi,        /* chi = MdagM(u) psi          */
		  IFloat  *mp_sq_p,    /* pointer to Sum |M psi|^2    */
		  IFloat  Kappa,       /* Wilson's kappa parameter    */
		  Wilson *wilson_p);  /* pointer to a Wilson struct. */

void wilson_dslash(IFloat *chi, 
		   IFloat *u, 
		   IFloat *psi, 
		   int cb,
		   int dag,
		   Wilson *wilson_p);

void wilson_m(IFloat *chi, 
	      IFloat *u, 
	      IFloat *psi, 
	      IFloat kappa,
	      Wilson *wilson_p);

void wilson_mdag(IFloat *chi, 
		 IFloat *u, 
		 IFloat *psi, 
		 IFloat kappa,
		 Wilson *wilson_p);
#ifdef __cplusplus
}
#endif

#endif
