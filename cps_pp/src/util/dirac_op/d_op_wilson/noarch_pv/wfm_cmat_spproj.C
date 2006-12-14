#include <config.h>
CPS_START_NAMESPACE
/*--------------------------------------------------------------------*/
/*! \file
  \brief  Routine used internally in the DiracOpWilson class.
*/
/*--------------------------------------------------------------------*/

/***************************************************************************/
/*                                                                         */
/* wfm_cmat_spproj:                                                        */
/*                                                                         */
/* Spin project psi with (1+sign*gamma) and do the color multiplication    */
/* of the resulting half spinors with the complex conjugate gauge field.   */
/* Put the result into the 4 backward half spinors ab[4].                  */
/*                                                                         */
/* ab       --> 2 component spproj'ed padded spinors for the backward      */
/*              direction                                                  */
/* u        --> gauge fields (not padded)                                  */
/* psi      --> 4 component spinors (not padded)                           */
/* sign     --> sign in spin project                                       */
/* wilson_p --> pointer to Wilson struct                                   */
/* cb       --> checkerboard                                               */
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
#define AB0(r,c,s)      *(ab0 + ab0_ptr_b +(r+2*(c+3*s)))
#define AB1(r,c,s)      *(ab1 + ab1_ptr_b +(r+2*(c+3*s)))
#define AB2(r,c,s)      *(ab2 + ab2_ptr_b +(r+2*(c+3*s)))
#define AB3(r,c,s)      *(ab3 + ab3_ptr_b +(r+2*(c+3*s)))

#define TMP0(r,c,s)     *(tmp0_spproj +(r+2*(c+3*s)))
#define TMP1(r,c,s)     *(tmp1_spproj +(r+2*(c+3*s)))
#define TMP2(r,c,s)     *(tmp2_spproj +(r+2*(c+3*s)))
#define TMP3(r,c,s)     *(tmp3_spproj +(r+2*(c+3*s)))

#define U(r,row,col,d)  *(u_p+(r+2*(row+3*(col+3*d))))


IFloat tmp0_spproj[HALF_SPINOR_SIZE];
IFloat tmp1_spproj[HALF_SPINOR_SIZE];
IFloat tmp2_spproj[HALF_SPINOR_SIZE];
IFloat tmp3_spproj[HALF_SPINOR_SIZE];

extern void wilson_dslash_csmatdag(Float *out, Float *u, Float *in, int mu);
extern void wilson_dslash_spproj(double *out0, double *out1, double *out2, double *out3, 
				 double *inf, double sign);


void wfm_cmat_spproj(IFloat *ab0, IFloat *ab1, IFloat *ab2, IFloat *ab3,
		     IFloat *u, 
		     IFloat *psi, 
		     IFloat sign, 
		     Wilson *wilson_p, 
		     int cb)
{
  char *fname = "wfm_cmat_spproj";
  int iyzt;
  int ix;
  int ip;
  int c;
  int s;
  int ixb;
  int sx;
  int bpad;
  int shift;
  int *ptr;
  int ab0_ptr_a, ab1_ptr_a, ab2_ptr_a, ab3_ptr_a;
  int ab0_ptr_b, ab1_ptr_b, ab2_ptr_b, ab3_ptr_b;
  IFloat *psi_p;
  IFloat *u_p;

  /*************************************************************************/
  /* Initializations                                                       */
  /*************************************************************************/

  /* Construct the base pointer for the array of pointers that holds the   */
  /* addresses of the ab spinors.                                          */
  bpad  = 1;
  shift = 0;
  ptr = wilson_p->ptr + (shift + 2 * (bpad + 2 * cb)) * wilson_p->offset;

  /*************************************************************************/
  /* Loop over sites                                                       */
  /*************************************************************************/
  psi_p = psi;
  u_p = u;
  ip = 0;
  for(iyzt=0; iyzt <= wilson_p->yztmax; iyzt++){

    ab0_ptr_a = ptr[ip++];
    ab1_ptr_a = ptr[ip++];
    ab2_ptr_a = ptr[ip++];
    ab3_ptr_a = ptr[ip++];
    sx        = ptr[ip++];

    for(ix=0; ix <= sx; ix++){

      ixb = ix * BLOCK;
      ab0_ptr_b = ab0_ptr_a + ixb;
      ab1_ptr_b = ab1_ptr_a + ixb;
      ab2_ptr_b = ab2_ptr_a + ixb;
      ab3_ptr_b = ab3_ptr_a + ixb;


      /*********************************************************************/
      /* Spin project to [1 + sign * gamma_mu]                             */
      /*********************************************************************/    
      wilson_dslash_spproj(tmp0_spproj,
			   tmp1_spproj,
			   tmp2_spproj,
			   tmp3_spproj,
			   psi_p,
			   sign);

      psi_p = psi_p + 24;


      /*********************************************************************/
      /* Color multiply the projected half spinors with U_mu^dag           */ 
      /*********************************************************************/
      /* mu = 0 */
      wilson_dslash_csmatdag((ab0+ab0_ptr_b), u_p, tmp0_spproj, 0);

      /* mu = 1 */
      wilson_dslash_csmatdag((ab1+ab1_ptr_b), u_p, tmp1_spproj, 1);

      /* mu = 2 */
      wilson_dslash_csmatdag((ab2+ab2_ptr_b), u_p, tmp2_spproj, 2);

      /* mu = 3 */
      wilson_dslash_csmatdag((ab3+ab3_ptr_b), u_p, tmp3_spproj, 3);

      u_p = u_p + 72;


    }
  }

}

CPS_END_NAMESPACE
