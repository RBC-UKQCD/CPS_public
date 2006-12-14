#include <config.h>
CPS_START_NAMESPACE
/*--------------------------------------------------------------------*/
/*! \file
  \brief  Routine used internally in the DiracOpWilson class.
*/
/*--------------------------------------------------------------------*/

/***************************************************************************/
/*                                                                         */
/* wfm_mat_trick:                                                          */
/*                                                                         */
/* Do the color multiplication for each forward half spinor af with the    */
/* gauge field. Then expand the half spinors to full spinors (trick)       */
/* add them and store their sum into the full spinor chi.                  */
/*                                                                         */
/* chi      --> 4 component spinors (not padded)                           */
/* u        --> gauge fields (not padded)                                  */
/* af       --> 2 component spproj'ed padded spinors for the backward      */
/*              direction                                                  */
/* sign     --> sign in spin project                                       */
/* wilson_p --> pointer to Wilson structure                                */
/* cb       --> chekerboard                                                */
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
#define AF0(r,c,s)      *(af0 + af0_ptr_b +(r+2*(c+3*s)))
#define AF1(r,c,s)      *(af1 + af1_ptr_b +(r+2*(c+3*s)))
#define AF2(r,c,s)      *(af2 + af2_ptr_b +(r+2*(c+3*s)))
#define AF3(r,c,s)      *(af3 + af3_ptr_b +(r+2*(c+3*s)))

#define TMP0(r,c,s)     *(tmp0_trick +(r+2*(c+3*s)))
#define TMP1(r,c,s)     *(tmp1_trick +(r+2*(c+3*s)))
#define TMP2(r,c,s)     *(tmp2_trick +(r+2*(c+3*s)))
#define TMP3(r,c,s)     *(tmp3_trick +(r+2*(c+3*s)))

#define U(r,row,col,d)  *(u_p+(r+2*(row+3*(col+3*d))))

#define STR(s) #s
#define QuadLoad(v,f)  asm volatile( "lfpdx " STR(f) ",0,%0" :: "r" (v) : "fr" STR(f) )


extern void wilson_dslash_mat_trick(double *outf, 
				    double *u, 
				    double *wfm_tmp0, 
				    double *wfm_tmp1, 
				    double *wfm_tmp2, 
				    double *wfm_tmp3, 
				    double *in0, 
				    double *in1, 
				    double *in2, 
				    double *in3,
				    double *in0p);


double wfm_tmp0[HALF_SPINOR_SIZE] __attribute__((aligned(BGL_L1_ALIGNSIZE)));
double wfm_tmp1[HALF_SPINOR_SIZE] __attribute__((aligned(BGL_L1_ALIGNSIZE)));
double wfm_tmp2[HALF_SPINOR_SIZE] __attribute__((aligned(BGL_L1_ALIGNSIZE)));
double wfm_tmp3[HALF_SPINOR_SIZE] __attribute__((aligned(BGL_L1_ALIGNSIZE)));


void wfm_mat_trick(IFloat *chi, 
		   IFloat *u,
		   IFloat *af0, IFloat *af1, IFloat *af2, IFloat *af3,
		   IFloat sign, 
		   Wilson *wilson_p, 
		   int cb)
{
  char *fname = "wfm_mat_trick";
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
  int af0_ptr_a, af1_ptr_a, af2_ptr_a, af3_ptr_a;
  int af0_ptr_b, af1_ptr_b, af2_ptr_b, af3_ptr_b, af00_ptr_b;
  IFloat *chi_p;
  IFloat *u_p;
  IFloat *u_pp;
  double sign_a[2] __attribute__((aligned(BGL_L1_ALIGNSIZE)));
  double zero[2]   __attribute__((aligned(BGL_L1_ALIGNSIZE)));
  double one[2]    __attribute__((aligned(BGL_L1_ALIGNSIZE)));

  /*************************************************************************/
  /* Initializations                                                       */
  /*************************************************************************/
  /* Register initializations */
  sign_a[0] = -sign;
  sign_a[1] = -sign;
  QuadLoad(sign_a, 29);

  zero[0] = 0.0;
  zero[1] = 0.0;
  QuadLoad(zero, 30);

  one[0] = 1.0;
  one[1] = 1.0;
  QuadLoad(one, 31);

  /* Construct the base pointer for the array of pointers that holds the   */
  /* addresses of the af spinors.                                          */
  bpad  = 0;
  shift = 1;
  ptr = wilson_p->ptr + (shift + 2 * (bpad + 2 * cb)) * wilson_p->offset;

  /*************************************************************************/
  /* Loop over sites                                                       */
  /*************************************************************************/
  chi_p = chi;
  u_p = u;
  ip = 0;
  for(iyzt=0; iyzt <= wilson_p->yztmax; iyzt++){

    af0_ptr_a = ptr[ip++];
    af1_ptr_a = ptr[ip++];
    af2_ptr_a = ptr[ip++];
    af3_ptr_a = ptr[ip++];
    sx        = ptr[ip++];


    for(ix=0; ix <= sx; ix++){

      ixb = ix * BLOCK;
      af0_ptr_b  = af0_ptr_a + ixb;
      af1_ptr_b  = af1_ptr_a + ixb;
      af2_ptr_b  = af2_ptr_a + ixb;
      af3_ptr_b  = af3_ptr_a + ixb;
      af00_ptr_b = af0_ptr_a + ixb + BLOCK;


      /*********************************************************************/
      /* Color multiply the af half spinors with U_mu and trick them       */
      /* into chi.                                                         */ 
      /*********************************************************************/
      wilson_dslash_mat_trick(chi_p,
			      u_p,
			      wfm_tmp0, 
			      wfm_tmp1, 
			      wfm_tmp2, 
			      wfm_tmp3, 
			      (af0+af0_ptr_b),
			      (af1+af1_ptr_b),
			      (af2+af2_ptr_b),
			      (af3+af3_ptr_b),
			      (af0+af00_ptr_b));

      u_p = u_p + 72;
      chi_p = chi_p + 24;

   }
  }

}

CPS_END_NAMESPACE
