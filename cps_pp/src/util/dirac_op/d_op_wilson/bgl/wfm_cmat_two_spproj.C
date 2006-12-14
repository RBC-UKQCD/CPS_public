#include <config.h>
CPS_START_NAMESPACE
/*--------------------------------------------------------------------*/
/*! \file
  \brief  Routine used internally in the DiracOpWilson class.
*/
/*--------------------------------------------------------------------*/

/***************************************************************************/
/*                                                                         */
/* wfm_cmat_two_spproj:                                                    */
/*                                                                         */
/* Combination routine of spproj and cmat_spproj in one volume loop        */
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
#include <sys/bgl/bgl_sys.h>
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

#define STR(s) #s
#define QuadLoad(v,f)  asm volatile( "lfpdx " STR(f) ",0,%0" :: "r" (v) : "fr" STR(f) )

extern void wilson_dslash_spproj(double *out0, 
				 double *out1, 
				 double *out2, 
				 double *out3, 
				 double *inf);

extern void wilson_dslash_cmat_spproj(double *out0, 
				      double *out1, 
				      double *out2, 
				      double *out3, 
				      double *u, 
				      double *inf);


void wfm_cmat_two_spproj(IFloat *af0, IFloat *af1, IFloat *af2, IFloat *af3,
			 IFloat *ab0, IFloat *ab1, IFloat *ab2, IFloat *ab3,
			 IFloat *u, 
			 IFloat *psi, 
			 IFloat sign, 
			 Wilson *wilson_p, 
			 int cb)
{
  char *fname = "wfm_cmat_spproj";
  int iyzt;
  int ix;
  int ip_f, ip_b;
  int c;
  int s;
  int ixb;
  int sx;
  int bpad;
  int shift;
  int *ptr_f;
  int *ptr_b;
  int ab0_ptr_a, ab1_ptr_a, ab2_ptr_a, ab3_ptr_a;
  int ab0_ptr_b, ab1_ptr_b, ab2_ptr_b, ab3_ptr_b;
  int af0_ptr_a, af1_ptr_a, af2_ptr_a, af3_ptr_a;
  int af0_ptr_b, af1_ptr_b, af2_ptr_b, af3_ptr_b;
  IFloat *psi_p;
  IFloat *psi_pp;
  IFloat *u_p;
  IFloat *u_pp;
  double sign_f[2] __attribute__((aligned(BGL_QUAD_ALIGNSIZE)));
  double sign_b[2] __attribute__((aligned(BGL_QUAD_ALIGNSIZE)));
  double zero[2]   __attribute__((aligned(BGL_QUAD_ALIGNSIZE)));
  double one[2]    __attribute__((aligned(BGL_QUAD_ALIGNSIZE)));

  /*************************************************************************/
  /* Initializations                                                       */
  /*************************************************************************/
  /* Register initializations */
  sign_f[0] = -sign;
  sign_f[1] = -sign;
  QuadLoad(sign_f, 28);

  sign_b[0] = sign;
  sign_b[1] = sign;
  QuadLoad(sign_b, 29);

  zero[0] = 0.0;
  zero[1] = 0.0;
  QuadLoad(zero, 30);

  one[0] = 1.0;
  one[1] = 1.0;
  QuadLoad(one, 31);

  /* Construct the base pointer for the array of pointers that holds the   */
  /* addresses of the ab spinors.                                          */
  bpad  = 1;
  shift = 0;
  ptr_b = wilson_p->ptr + (shift + 2 * (bpad + 2 * cb)) * wilson_p->offset;


  /* Construct the base pointer for the array of pointers that holds the   */
  /* addresses of the af spinors.                                          */
  bpad = 0;
  shift = 0;
  ptr_f = wilson_p->ptr + (shift + 2 * (bpad + 2 * cb)) * wilson_p->offset;


  /*************************************************************************/
  /* Loop over sites                                                       */
  /*************************************************************************/
  psi_p = psi;
  u_p   = u;
  ip_f  = 0;
  ip_b  = 0;

  for(iyzt=0; iyzt <= wilson_p->yztmax; iyzt++){

    af0_ptr_a = ptr_f[ip_f++];
    af1_ptr_a = ptr_f[ip_f++];
    af2_ptr_a = ptr_f[ip_f++];
    af3_ptr_a = ptr_f[ip_f++];
    sx        = ptr_f[ip_f++];

    ab0_ptr_a = ptr_b[ip_b++];
    ab1_ptr_a = ptr_b[ip_b++];
    ab2_ptr_a = ptr_b[ip_b++];
    ab3_ptr_a = ptr_b[ip_b++];
    sx        = ptr_b[ip_b++];

    for(ix=0; ix <= sx; ix++){

      ixb = ix * BLOCK;

      af0_ptr_b = af0_ptr_a + ixb;
      af1_ptr_b = af1_ptr_a + ixb;
      af2_ptr_b = af2_ptr_a + ixb;
      af3_ptr_b = af3_ptr_a + ixb;

      ab0_ptr_b = ab0_ptr_a + ixb;
      ab1_ptr_b = ab1_ptr_a + ixb;
      ab2_ptr_b = ab2_ptr_a + ixb;
      ab3_ptr_b = ab3_ptr_a + ixb;

      wilson_dslash_spproj((af0 + af0_ptr_b),
			   (af1 + af1_ptr_b),
			   (af2 + af2_ptr_b),
			   (af3 + af3_ptr_b),
			   psi_p);

      wilson_dslash_cmat_spproj((ab0 + ab0_ptr_b),
				(ab1 + ab1_ptr_b),
				(ab2 + ab2_ptr_b),
				(ab3 + ab3_ptr_b),
				u_p,
				psi_p);

      psi_p = psi_p + 24;
      u_p   = u_p   + 72;

    }
  }

}

CPS_END_NAMESPACE
