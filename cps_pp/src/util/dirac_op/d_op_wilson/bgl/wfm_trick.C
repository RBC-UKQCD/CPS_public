#include <config.h>
CPS_START_NAMESPACE
/*--------------------------------------------------------------------*/
/*! \file
  \brief  Routine used internally in the DiracOpWilson class.
*/
/*--------------------------------------------------------------------*/

/***************************************************************************/
/*                                                                         */
/* wfm_trick:                                                              */
/*                                                                         */
/* Expand the backward half spinors ab to full spinors (trick) add them    */
/* and then add their sum into the full spinor chi.                        */
/*                                                                         */
/* chi      --> 4 component spinors (not padded)                           */
/* ab       --> 2 component spproj'ed padded spinors for the backward      */
/*              direction                                                  */
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
#include <sys/bgl/bgl_sys.h>
#include <comms/scu.h>
CPS_START_NAMESPACE

#define CHI(r,c,s)      *(chi_p +(r+2*(c+3*s)))
#define AB0(r,c,s)      *(ab0 + ab0_ptr_b +(r+2*(c+3*s)))
#define AB1(r,c,s)      *(ab1 + ab1_ptr_b +(r+2*(c+3*s)))
#define AB2(r,c,s)      *(ab2 + ab2_ptr_b +(r+2*(c+3*s)))
#define AB3(r,c,s)      *(ab3 + ab3_ptr_b +(r+2*(c+3*s)))

#define STR(s) #s
#define QuadLoad(v,f)  asm volatile( "lfpdx " STR(f) ",0,%0" :: "r" (v) : "fr" STR(f) )
#define BGL_QUAD_ALIGNSIZE  16


extern void wilson_dslash_trick(double *outf, 
				double *in0, double *in1, double *in2, double *in3, 
				double sign, int accum);


void wfm_trick(IFloat *chi, 
	       IFloat *ab0, IFloat *ab1, IFloat *ab2, IFloat *ab3,
	       IFloat sign, 
	       Wilson *wilson_p, 
	       int cb)
{
  char *fname = "wfm_trick";


  int iyzt;
  int ix;
  int ip;
  int c;
  int ixb;
  int sx;
  int bpad;
  int shift;
  int *ptr;
  int ab0_ptr_a, ab1_ptr_a, ab2_ptr_a, ab3_ptr_a;
  int ab0_ptr_b, ab1_ptr_b, ab2_ptr_b, ab3_ptr_b;
  IFloat *chi_p;
  double sign_a[2] __attribute__((aligned(BGL_QUAD_ALIGNSIZE)));
  double zero[2]   __attribute__((aligned(BGL_QUAD_ALIGNSIZE)));
  double one[2]    __attribute__((aligned(BGL_QUAD_ALIGNSIZE)));

  /*************************************************************************/
  /* Initializations                                                       */
  /*************************************************************************/
  /* Register initializations */
  sign_a[0] = sign;
  sign_a[1] = sign;
  QuadLoad(sign_a, 29);

  zero[0] = 0.0;
  zero[1] = 0.0;
  QuadLoad(zero, 30);

  one[0] = 1.0;
  one[1] = 1.0;
  QuadLoad(one, 31);

  /* Construct the base pointer for the array of pointers that holds the   */
  /* addresses of the ab spinors.                                          */
  bpad  = 1;
  shift = 1;
  ptr = wilson_p->ptr + (shift + 2 * (bpad + 2 * cb)) * wilson_p->offset;

  /*************************************************************************/
  /* Loop over sites                                                       */
  /*************************************************************************/
  chi_p = chi;
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
      /* Expand the half spinors in to a full spinor (trick), where        */
      /* the spin projection was done with [1+sign*gamma].                 */
      /* Aaccumulate the result into the full spinor (accum=1).            */
      /*********************************************************************/
      wilson_dslash_trick(chi_p,
			  (ab0 + ab0_ptr_b),
			  (ab1 + ab1_ptr_b),
			  (ab2 + ab2_ptr_b),
			  (ab3 + ab3_ptr_b),
			  sign, 1);

      chi_p = chi_p + 24;

    }
  }




}

CPS_END_NAMESPACE
