#include <config.h>
CPS_START_NAMESPACE
/*--------------------------------------------------------------------*/
/*! \file
  \brief  Routine used internally in the DiracOpWilson class.
*/
/*--------------------------------------------------------------------*/

/***************************************************************************/
/*                                                                         */
/* wfm_spproj:                                                             */
/*                                                                         */
/* Spin project psi with (1- sign*gamma) into the 4 forward half spinors   */
/* af[4]. It performs the spin projection by equating the top 2 spin       */
/* components  of [(1 - sign*gamma_mu)/2 psi] to af_mu. The bottom two     */
/* components can be reconstructed  using the trick routines.              */
/*                                                                         */
/* wfm_spproj(float *af0,                 ; af for mu = 0                  */
/*	      float *af1,                 ; af for mu = 1                  */
/*	      float *af2,                 ; af for mu = 2                  */
/*   	      float *af3,                 ; af for mu = 3                  */
/* 	      float *psi,		  ; 4 comp. spinor                 */
/*	      float sign, 		  ; the +/- sign                   */
/*	      Wilson *wilson_p,           ; Wilson struct.                 */
/*	      int cb);		          ; checkerboard 0/1 = even/odd    */
/*                                                                         */
/* WARNING:                                                                */
/*                                                                         */
/* This set of routines will work only if the node sublattices have        */
/* even number of sites in each direction.                                 */
/*                                                                         */
/***************************************************************************/

CPS_END_NAMESPACE
#include <util/data_types.h>
#include <util/vector.h>
#include <util/wilson.h>
#include <util/error.h>
#include <comms/scu.h>
CPS_START_NAMESPACE

#define PSI(r,c,s)      *(psi_p +(r+2*(c+3*s)))
#define AF0(r,c,s)      *(af0 + af0_ptr_b +(r+2*(c+3*s)))
#define AF1(r,c,s)      *(af1 + af1_ptr_b +(r+2*(c+3*s)))
#define AF2(r,c,s)      *(af2 + af2_ptr_b +(r+2*(c+3*s)))
#define AF3(r,c,s)      *(af3 + af3_ptr_b +(r+2*(c+3*s)))

extern void wilson_dslash_spproj(double *out0, double *out1, double *out2, double *out3, 
				 double *inf, double sign);


void wfm_spproj(IFloat *af0, IFloat *af1, IFloat *af2, IFloat *af3,
		IFloat *psi, 
		IFloat sign, 
		Wilson *wilson_p, 
		int cb)
{
  char *fname = "wfm_spproj";
  int iyzt;
  int ix;
  int ip;
  int c;
  int ixb;
  int sx;
  int bpad;
  int shift;
  int *ptr;
  int af0_ptr_a, af1_ptr_a, af2_ptr_a, af3_ptr_a;
  int af0_ptr_b, af1_ptr_b, af2_ptr_b, af3_ptr_b;
  IFloat *psi_p;

  /*************************************************************************/
  /* Initializations                                                       */
  /*************************************************************************/

  /* Construct the base pointer for the array of pointers that holds the   */
  /* addresses of the af spinors.                                          */
  bpad = 0;
  shift = 0;
  ptr = wilson_p->ptr + (shift + 2 * (bpad + 2 * cb)) * wilson_p->offset;

  /*************************************************************************/
  /* Loop over sites                                                       */
  /*************************************************************************/
  psi_p = psi;
  ip = 0;
  for(iyzt=0; iyzt <= wilson_p->yztmax; iyzt++){

    af0_ptr_a = ptr[ip++];
    af1_ptr_a = ptr[ip++];
    af2_ptr_a = ptr[ip++];
    af3_ptr_a = ptr[ip++];
    sx        = ptr[ip++];

    for(ix=0; ix <= sx; ix++){

      ixb = ix * BLOCK;
      af0_ptr_b = af0_ptr_a + ixb;
      af1_ptr_b = af1_ptr_a + ixb;
      af2_ptr_b = af2_ptr_a + ixb;
      af3_ptr_b = af3_ptr_a + ixb;

      
      /*********************************************************************/
      /* Spin project to [1 - sign * gamma_mu]                             */
      /*********************************************************************/    
      wilson_dslash_spproj((af0 + af0_ptr_b),
			   (af1 + af1_ptr_b),
			   (af2 + af2_ptr_b),
			   (af3 + af3_ptr_b),
			   psi_p,
			   -sign);

      psi_p = psi_p + 24;

    }
  }

}

CPS_END_NAMESPACE
